#!/usr/bin/env python3
"""Build hg19 with chrM = rCRS (NC_012920.1, 16,569 bp), matched to a BAM header.

UCSC hg19.fa carries chrM = NC_001807 (Yoruba mtDNA, 16,571 bp). Ion Torrent /
Genexus hg19 uses the rCRS for chrM. This script:

  1. downloads UCSC hg19.fa.gz (md5-checked against UCSC md5sum.txt) and the
     rCRS from NCBI (length-checked), unless local copies are given;
  2. replaces chrM with the rCRS, uppercases all bases (unless --keep-case);
  3. writes the contigs named in the BAM's @SQ lines, in the BAM's order;
  4. checks every contig's length, and its MD5 if the BAM header has M5 tags;
  5. writes the .fai index itself (samtools faidx-compatible).

Exit 0 only if every @SQ contig is present with matching length (and M5 when
given). Needs python3, plus samtools on PATH or pysam to read the BAM header.

Usage:
  python3 -m hpoly_vaf prepare-genexus-reference --bam merged.bam --out hg19_rcrs.fa
  python3 -m hpoly_vaf prepare-genexus-reference --bam merged.bam --out hg19_rcrs.fa \
      --hg19 hg19.fa.gz --mt NC_012920.1.fa      # offline, local copies
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import urllib.request

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/"
NCBI = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text")
RCRS_LEN = 16569
RCRS_MD5 = "c68f52674c9fb33aef52dcf399755519"  # SAM M5 of uppercase NC_012920.1
LINE = 60


def die(msg):
    sys.exit("ERROR: " + msg)


def fetch(url, dest):
    print(f"download {url}", file=sys.stderr)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
    os.replace(tmp, dest)


def hash_file(path, algorithm="md5"):
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def bam_header_text(bam):
    try:
        return subprocess.run(["samtools", "view", "-H", bam], check=True,
                              capture_output=True, text=True).stdout
    except FileNotFoundError:
        pass
    try:
        import pysam
    except ImportError:
        die("need samtools on PATH or the pysam module to read the BAM header")
    with pysam.AlignmentFile(bam, "rb") as b:
        return str(b.header)


def bam_sq(bam):
    hdr = bam_header_text(bam)
    sq = []
    for line in hdr.splitlines():
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            sq.append((f["SN"], int(f["LN"]), f.get("M5", "").lower()))
    if not sq:
        die(f"no @SQ lines in {bam}")
    return sq


def read_fasta_one(path):
    records = list(fasta_records(path))
    if len(records) != 1:
        die("rCRS FASTA must contain exactly one record")
    return records[0][1]


def fasta_records(path):
    op = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def add_arguments(ap):
    """Share argument definitions between module and pipeline entry points."""
    ap.add_argument("--bam", required=True, help="BAM whose @SQ lines define the target")
    ap.add_argument("--out", required=True, help="output FASTA")
    ap.add_argument("--hg19", help="local UCSC hg19.fa or hg19.fa.gz (skip download)")
    ap.add_argument("--mt", help="local rCRS FASTA, NC_012920.1 (skip download)")
    ap.add_argument("--trimmed-bam", help="also validate against the matching trimmed BAM header")
    ap.add_argument("--hg19-md5-file", help="trusted UCSC md5sum.txt for a local hg19.fa.gz")
    ap.add_argument("--cache", help="download directory (default: reference_cache beside output)")
    ap.add_argument("--keep-case", action="store_true",
                    help="keep UCSC soft-masking (lowercase); default uppercases")
    ap.add_argument("--mt-name", default="chrM", help="contig to replace (default chrM)")


def prepare_reference(a):
    """Build a separate validated reference; never modify inputs or bypass F1."""
    output = Path(a.out)
    for suffix in ("", ".fai", ".manifest.json"):
        if Path(str(output) + suffix).exists():
            die(f"refusing to overwrite {output}{suffix}")
    output.parent.mkdir(parents=True, exist_ok=True)
    a.cache = a.cache or str(output.parent / "reference_cache")
    os.makedirs(a.cache, exist_ok=True)

    sq = bam_sq(a.bam)
    if len({name for name, _, _ in sq}) != len(sq):
        die("duplicate SQ names in original BAM")
    if a.mt_name not in {name for name, _, _ in sq}:
        die(f"BAM lacks {a.mt_name}; this command is only for explicit hg19/rCRS preparation")
    if a.trimmed_bam:
        tsq = bam_sq(a.trimmed_bam)
        original_lengths = {name: length for name, length, _ in sq}
        if any(original_lengths.get(name) != length for name, length, _ in tsq):
            die("trimmed BAM contigs/lengths are incompatible with original BAM")
    else:
        tsq = []
    print(f"BAM header: {len(sq)} contigs", file=sys.stderr)

    # --- inputs -------------------------------------------------------------
    hg19 = a.hg19
    checksum_status = "local input; UCSC archive checksum not supplied"
    if not hg19:
        hg19 = os.path.join(a.cache, "hg19.fa.gz")
        if not os.path.exists(hg19):
            fetch(UCSC + "hg19.fa.gz", hg19)
        sums = os.path.join(a.cache, "hg19.md5sum.txt")
        if not os.path.exists(sums):
            fetch(UCSC + "md5sum.txt", sums)
        want = ucsc_md5(sums)
        got = hash_file(hg19)
        if got != want:
            die(f"hg19.fa.gz md5 {got} != UCSC {want}; delete it and rerun")
        print(f"hg19.fa.gz md5 OK ({got})", file=sys.stderr)
        checksum_status = "verified against UCSC md5sum.txt"
    elif a.hg19_md5_file:
        if hash_file(hg19) != ucsc_md5(a.hg19_md5_file):
            die("local hg19 archive differs from supplied UCSC checksum")
        checksum_status = "verified against supplied UCSC md5sum.txt"

    mt = a.mt
    if not mt:
        mt = os.path.join(a.cache, "NC_012920.1.fa")
        if not os.path.exists(mt):
            fetch(NCBI, mt)
    rcrs = read_fasta_one(mt).upper()
    if len(rcrs) != RCRS_LEN:
        die(f"{mt}: length {len(rcrs)} != {RCRS_LEN}; not NC_012920.1")
    if hashlib.md5(rcrs.encode()).hexdigest() != RCRS_MD5:
        die(f"{mt}: MD5 differs from NC_012920.1")
    print(f"rCRS OK: {RCRS_LEN} bp, M5 {RCRS_MD5}", file=sys.stderr)

    # --- stream hg19, replace chrM, keep only BAM contigs -------------------
    # Each kept contig goes to its own temp file so the output can be written
    # in BAM order without holding the genome in memory.
    want = {n for n, _, _ in sq}
    tmpdir = tempfile.mkdtemp(prefix=output.name + ".parts-", dir=output.parent)
    seen, problems = {}, []
    for name, seq in fasta_records(hg19):
        if name not in want:
            continue
        if name in seen:
            shutil.rmtree(tmpdir)
            die(f"duplicate FASTA contig: {name}")
        if name == a.mt_name:
            print(f"{name}: UCSC {len(seq)} bp replaced by rCRS {RCRS_LEN} bp",
                  file=sys.stderr)
            seq = rcrs
        elif not a.keep_case:
            seq = seq.upper()
        seen[name] = (len(seq), hashlib.md5(seq.upper().encode()).hexdigest())
        with open(os.path.join(tmpdir, str(len(seen) - 1)), "w") as f:
            for i in range(0, len(seq), LINE):
                f.write(seq[i:i + LINE] + "\n")
        del seq

    # --- check against BAM header -------------------------------------------
    for n, ln, m5 in sq + tsq:
        if n not in seen:
            problems.append(f"{n}: in BAM header, absent from hg19.fa")
            continue
        gl, gm = seen[n]
        if gl != ln:
            problems.append(f"{n}: length {gl} != BAM {ln}")
        elif m5 and gm != m5:
            problems.append(f"{n}: M5 {gm} != BAM {m5}")
    n_m5 = sum(1 for s_ in sq if s_[2])
    print(f"M5 tags in BAM header: {n_m5}/{len(sq)} "
          f"({'sequence identity checked' if n_m5 else 'lengths only'})",
          file=sys.stderr)
    if problems:
        shutil.rmtree(tmpdir)
        print("\n".join(problems), file=sys.stderr)
        die(f"{len(problems)} contig(s) do not match the BAM header; no output written")

    # --- concatenate in BAM order, write .fai --------------------------------
    parts = {name: str(i) for i, name in enumerate(seen)}
    off = 0
    with open(a.out, "wb") as out, open(a.out + ".fai", "w") as fai:
        for n, ln, _ in sq:
            head = f">{n}\n".encode()
            out.write(head)
            off += len(head)
            width = min(LINE, ln)
            fai.write(f"{n}\t{ln}\t{off}\t{width}\t{width + 1}\n")
            with open(os.path.join(tmpdir, parts[n]), "rb") as f:
                shutil.copyfileobj(f, out)
            off += os.path.getsize(os.path.join(tmpdir, parts[n]))
    shutil.rmtree(tmpdir)
    manifest = {
        "operation": "explicit hg19 chrM replacement with NC_012920.1",
        "bam": os.path.abspath(a.bam),
        "trimmed_bam": os.path.abspath(a.trimmed_bam) if a.trimmed_bam else None,
        "hg19": {"path": os.path.abspath(hg19), "sha256": hash_file(hg19, "sha256"),
                 "archive_checksum_status": checksum_status},
        "rcrs": {"path": os.path.abspath(mt), "accession": "NC_012920.1",
                 "length": RCRS_LEN, "sequence_md5": RCRS_MD5},
        "output": {"path": str(output.resolve()), "sha256": hash_file(a.out, "sha256"),
                   "uppercase": not a.keep_case},
        "contigs": [{"name": name, "length": length, "sequence_md5": seen[name][1],
                     "bam_m5": m5 or None} for name, length, m5 in sq],
        "original_m5_checked": n_m5,
        "trimmed_m5_checked": sum(bool(m5) for _, _, m5 in tsq),
        "identity_limit": "Only BAM contigs carrying M5 tags establish sequence identity to the alignment reference.",
    }
    Path(a.out + ".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"OK: {a.out} + .fai ({len(sq)} contigs, BAM order, "
          f"{'soft-masked' if a.keep_case else 'uppercase'}, "
          f"{a.mt_name}={RCRS_LEN} bp)")
    return 0


def ucsc_md5(path):
    """Never claim checksum verification when the manifest lacks the archive."""
    with open(path) as handle:
        values = [line.split()[0] for line in handle if len(line.split()) >= 2
                  and line.split()[-1].lstrip("*") == "hg19.fa.gz"]
    if len(values) != 1:
        die("UCSC checksum manifest must contain exactly one hg19.fa.gz entry")
    return values[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    return prepare_reference(ap.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
