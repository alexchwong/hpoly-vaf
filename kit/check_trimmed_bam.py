#!/usr/bin/env python3
"""
check_trimmed_bam.py - does the primer-trimmed BAM keep full-length reads?

Usage:
    python3 check_trimmed_bam.py --bam merged.bam --trimmed-bam merged.bam.ptrim.bam [--n 2000]
                                 [--region chr1:43800000-43830000]

It takes the first N primary reads of the trimmed BAM (or of --region), finds the same
reads in the original BAM by name, and compares them. Nothing is written or modified.

Needs either `samtools` on PATH or the Python package `pysam`. Both BAMs need an index
unless --region is omitted AND the reads sit at the very start of the file; if a fetch
fails, run `samtools index <bam>` first.

Exit status: 0 = SOFT (full reads kept), 1 = HARD or MIXED, 2 = could not decide.
"""
import argparse, re, shutil, subprocess, sys
from collections import Counter

SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x800          # unmapped, secondary, QC fail, supplementary
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def sam_lines(path, region=None):
    """Yield SAM text lines for primary mapped reads, using samtools if present, else pysam."""
    if shutil.which("samtools"):
        cmd = ["samtools", "view", "-F", str(SKIP_FLAGS), path] + ([region] if region else [])
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in p.stdout:
            yield line.rstrip("\n")
        p.stdout.close()
        p.wait()
        if p.returncode not in (0, -13, 141):       # -13/141: we stopped reading early
            sys.exit("samtools failed on %s %s:\n%s" % (path, region or "", p.stderr.read()))
    else:
        try:
            import pysam
        except ImportError:
            sys.exit("Need samtools on PATH or the pysam package.")
        with pysam.AlignmentFile(path, "rb") as f:
            it = f.fetch(region=region) if region else f.fetch(until_eof=True)
            for a in it:
                if a.flag & SKIP_FLAGS:
                    continue
                yield a.to_string()


def parse(line):
    f = line.split("\t")
    cig = CIG.findall(f[5])
    qlen = sum(int(n) for n, op in cig if op in "MIS=X")
    rlen = sum(int(n) for n, op in cig if op in "MDN=X")
    pos = int(f[3]) - 1
    return {
        "name": f[0], "rev": bool(int(f[1]) & 16), "chrom": f[2], "pos": pos, "end": pos + rlen,
        "cigar": f[5], "seq": f[9], "qlen": qlen,
        "hard": sum(int(n) for n, op in cig if op == "H"),
        "soft": sum(int(n) for n, op in cig if op == "S"),
        "tags": {t.split(":", 1)[0] for t in f[11:]},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="original flow-bearing BAM (merged.bam)")
    ap.add_argument("--trimmed-bam", required=True, help="primer-trimmed BAM (merged.bam.ptrim.bam)")
    ap.add_argument("--n", type=int, default=2000, help="number of trimmed reads to sample (default 2000)")
    ap.add_argument("--region", help="optional region, e.g. chr1:43800000-43830000")
    a = ap.parse_args()

    # 1. sample the trimmed BAM
    trimmed = []
    for line in sam_lines(a.trimmed_bam, a.region):
        r = parse(line)
        if trimmed and r["chrom"] != trimmed[0]["chrom"]:
            break                                   # keep the sample on one contig
        trimmed.append(r)
        if len(trimmed) >= a.n:
            break
    if not trimmed:
        print("No primary mapped reads found in the trimmed BAM" + (" in " + a.region if a.region else "") + ".")
        return 2

    chrom = trimmed[0]["chrom"]
    lo = max(0, min(r["pos"] for r in trimmed) - 500)
    hi = max(r["end"] for r in trimmed) + 500
    region = "%s:%d-%d" % (chrom, lo + 1, hi)

    # 2. fetch the same stretch of the original BAM
    orig = {}
    dup = 0
    for line in sam_lines(a.bam, region):
        r = parse(line)
        if r["name"] in orig:
            dup += 1
        orig[r["name"]] = r

    # 3. compare
    c = Counter()
    tag_t, tag_o = Counter(), Counter()
    examples = []
    for t in trimmed:
        c["trimmed_reads"] += 1
        c["trimmed_with_H"] += t["hard"] > 0
        c["trimmed_with_S"] += t["soft"] > 0
        c["trimmed_cigar_len_ne_seq"] += t["qlen"] != len(t["seq"])
        for tag in ("ZM", "ZP", "ZF", "ZL", "RG"):
            tag_t[tag] += tag in t["tags"]
        o = orig.get(t["name"])
        if o is None:
            c["no_original_found"] += 1
            continue
        c["matched"] += 1
        for tag in ("ZM", "ZP", "ZF", "RG"):
            tag_o[tag] += tag in o["tags"]
        c["strand_differs"] += t["rev"] != o["rev"]
        c["span_inside_original"] += (t["pos"] >= o["pos"] and t["end"] <= o["end"])
        c["span_identical"] += (t["pos"] == o["pos"] and t["end"] == o["end"])
        if t["seq"] == o["seq"]:
            c["seq_identical"] += 1
        elif len(t["seq"]) < len(o["seq"]) and t["seq"] in o["seq"]:
            c["seq_shorter_substring"] += 1
        else:
            c["seq_other_difference"] += 1
        if len(examples) < 5 and (t["pos"], t["cigar"]) != (o["pos"], o["cigar"]):
            examples.append((t["name"], o, t))

    n, m = c["trimmed_reads"], c["matched"]
    print("Sampled region            : %s" % region)
    print("Trimmed reads sampled     : %d" % n)
    print("  with hard clip (H)      : %d" % c["trimmed_with_H"])
    print("  with soft clip (S)      : %d" % c["trimmed_with_S"])
    print("  CIGAR length != SEQ     : %d" % c["trimmed_cigar_len_ne_seq"])
    print("  tags present            : " + ", ".join("%s=%d" % (k, tag_t[k]) for k in ("ZM", "ZP", "ZF", "ZL", "RG")))
    print("Matched in original BAM   : %d   (not found: %d; duplicate names in original: %d)" % (m, c["no_original_found"], dup))
    print("  original tags present   : " + ", ".join("%s=%d" % (k, tag_o[k]) for k in ("ZM", "ZP", "ZF", "RG")))
    print("  SEQ identical           : %d" % c["seq_identical"])
    print("  SEQ shorter (substring) : %d" % c["seq_shorter_substring"])
    print("  SEQ differs otherwise   : %d" % c["seq_other_difference"])
    print("  strand differs          : %d" % c["strand_differs"])
    print("  span inside original    : %d" % c["span_inside_original"])
    print("  span unchanged          : %d" % c["span_identical"])
    if examples:
        print("\nExamples where POS or CIGAR changed (original -> trimmed):")
        for name, o, t in examples:
            print("  %s\n    original pos=%d cigar=%s seqlen=%d\n    trimmed  pos=%d cigar=%s seqlen=%d"
                  % (name, o["pos"] + 1, o["cigar"], len(o["seq"]), t["pos"] + 1, t["cigar"], len(t["seq"])))

    print()
    if m == 0:
        print("VERDICT: UNDECIDED - no sampled trimmed read was found in the original BAM. "
              "Check that the two files are from the same run and that the original BAM is indexed.")
        return 2
    if c["seq_identical"] == m and c["trimmed_with_H"] == 0:
        print("VERDICT: SOFT - the trimmed BAM keeps the full read sequence; trimming is soft clipping only.")
        return 0
    if c["seq_shorter_substring"] + c["seq_identical"] == m and c["seq_shorter_substring"] > 0:
        print("VERDICT: HARD - the trimmed BAM stores shortened reads (%d of %d). "
              "Only its reference span can be used." % (c["seq_shorter_substring"], m))
        return 1
    print("VERDICT: MIXED - reads differ in ways other than simple trimming. Send this output back for review.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
