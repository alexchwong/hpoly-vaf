"""Real-file preflight checks F1--F11 from extraction specification v3."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
import json
from pathlib import Path
import sys

import pysam

from .cigar import QUERY
from .islands import iter_islands, primary_mapped
from .join import join_island


class HardFailure(Exception):
    """A recorded hard check failed; no subsequent phase work is allowed."""


class CheckLog:
    def __init__(self, work):
        self.path = Path(work) / "preflight.json"
        self.entries = []

    def add(self, name, passed, observed, expected, kind="HARD"):
        entry = dict(name=name, status="PASS" if passed else "FAIL",
                     observed=observed, expected=expected, kind=kind)
        self.entries.append(entry)
        self.path.write_text(json.dumps(self.entries, indent=2) + "\n")
        print(f"{entry['status']} {name}", flush=True)
        if not passed and kind == "HARD":
            raise HardFailure(name)


def ensure_index(bam, work, label):
    """Reuse a fresh adjacent BAI/CSI, otherwise build a separate work index."""
    bam, work = Path(bam), Path(work)
    candidates = [Path(str(bam) + ext) for ext in (".bai", ".csi")]
    candidates += [bam.with_suffix(ext) for ext in (".bai", ".csi")]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_mtime >= bam.stat().st_mtime:
            return str(candidate)
    index = work / (label + ".csi")
    pysam.index("-c", "-o", str(index), str(bam))
    return str(index)


def sort_observation(bam):
    """Inspect the first 100,000 records, with unmapped tail last."""
    with pysam.AlignmentFile(bam, "rb") as handle:
        header_ok = handle.header.to_dict().get("HD", {}).get("SO") == "coordinate"
        previous = None
        violations = 0
        checked = 0
        for read in islice(handle, 100000):
            key = (read.reference_id, read.reference_start) if read.reference_id >= 0 else (len(handle.references), 0)
            violations += previous is not None and key < previous
            previous = key
            checked += 1
    return dict(header_coordinate=header_ok, records=checked, regressions=violations)


def primary_sample(bam):
    with pysam.AlignmentFile(bam, "rb") as handle:
        return list(islice((r for r in handle if primary_mapped(r)), 10000))


def clip_agrees(agreements, joined):
    """Exact integer comparison avoids floating point at the 99.9% boundary."""
    return joined > 0 and agreements * 1000 >= joined * 999


def check_headers_and_samples(bam, trimmed, reference, log):
    with pysam.FastaFile(reference) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    mismatches = []
    for label, path in (("original", bam), ("trimmed", trimmed)):
        with pysam.AlignmentFile(path, "rb") as handle:
            for name, length in zip(handle.references, handle.lengths):
                if lengths.get(name) != length:
                    mismatches.append(dict(bam=label, contig=name, bam_length=length,
                                           fasta_length=lengths.get(name)))
    log.add("F1 contigs", not mismatches, mismatches, "Every BAM SQ name/length matches FASTA .fai")
    orders = {label: sort_observation(path) for label, path in (("original", bam), ("trimmed", trimmed))}
    log.add("F2 sort order", all(v["header_coordinate"] and not v["regressions"] for v in orders.values()),
            orders, "SO:coordinate and non-decreasing coordinates in first 100000 records of both BAMs")
    with pysam.AlignmentFile(bam, "rb") as handle:
        groups = handle.header.to_dict().get("RG", [])
    samples = sorted({g["SM"] for g in groups if "SM" in g})
    missing = [{"id": g.get("ID"), "missing": [k for k in ("FO", "KS", "SM") if k not in g]}
               for g in groups if any(k not in g for k in ("FO", "KS", "SM"))]
    log.add("F3 read groups", not missing and len(samples) == 1,
            dict(groups=len(groups), samples=samples, missing=missing), "Every RG has FO, KS, SM; exactly one distinct SM")
    rg = {g["ID"]: g for g in groups if "ID" in g}
    reads = primary_sample(bam)
    invalid = 0
    for read in reads:
        group = read.get_tag("RG") if read.has_tag("RG") else None
        invalid += (not read.query_name or ":" not in read.query_name or group not in rg
                    or read.query_name.split(":", 1)[0] != group.split(".", 1)[0])
    log.add("F4 run id", bool(reads) and invalid == 0,
            dict(checked=len(reads), invalid=invalid), "First 10000 primary mapped reads: name prefix equals RG prefix")
    missing_tags = Counter(tag for read in reads for tag in ("ZM", "ZP", "ZF", "RG") if not read.has_tag(tag))
    log.add("F5 flow tags", not missing_tags, dict(checked=len(reads), missing=dict(missing_tags)),
            "Every sampled original read has ZM, ZP, ZF, RG")
    reads = primary_sample(trimmed)
    hard = sum(any(op == 5 for op, _ in (r.cigartuples or [])) for r in reads)
    unequal = sum(r.query_sequence is None or sum(n for op, n in (r.cigartuples or []) if op in QUERY)
                  != len(r.query_sequence) for r in reads)
    log.add("F6 trimmed SEQ", bool(reads) and hard == 0 and unequal == 0,
            dict(checked=len(reads), hard_clips=hard, query_length_mismatches=unequal),
            "No H and CIGAR query length equals len(SEQ) in first 10000 trimmed primary mapped reads")


def check_islands(bam, trimmed, bam_index, trimmed_index, work, log):
    with pysam.AlignmentFile(trimmed, "rb") as handle:
        islands = list(iter_islands(handle))
    with (work / "islands.bed").open("w") as handle:
        for island in islands:
            handle.write(f"{island.contig}\t{island.start}\t{island.end}\n")
    totals = Counter()
    undefined = Counter()
    mismatches, examples, over_capacity = [], [], []
    depth = 0
    tasks = ((bam, trimmed, bam_index, trimmed_index, island) for island in islands)
    with ProcessPoolExecutor(max_workers=4) as pool:
        for island, result in zip(islands, pool.map(join_island, tasks, chunksize=1)):
            for key in ("joined", "violations", "duplicates", "agreements", "empty_after_clip", "no_ar"):
                totals[key] += result[key]
            undefined.update(result["undefined"])
            mismatches.extend(result["mismatches"][:max(0, 20 - len(mismatches))])
            examples.extend(result["examples"][:max(0, 20 - len(examples))])
            depth = max(depth, result["maximum_depth"])
            if result["maximum_depth"] > 200000:
                over_capacity.append(dict(contig=island.contig, start=island.start, end=island.end,
                                          maximum_depth=result["maximum_depth"]))
    diagnostics = dict(totals=totals, undefined_cases=undefined, first_20_mismatches=mismatches,
                       first_20_join_violations=examples, maximum_depth=depth, over_capacity=over_capacity)
    (work / "join_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    log.add("F7 join", totals["violations"] == 0 and bool(islands),
            dict(islands=len(islands), joined=totals["joined"], violations=totals["violations"], examples=examples),
            "Every trimmed primary read has exactly one matching original; AR inside original span; zero violations")
    log.add("F8 unique names", totals["duplicates"] == 0,
            dict(duplicate_occurrences=totals["duplicates"]), "Zero duplicate primary names within any island in either BAM")
    log.add("F9 clip agreement", clip_agrees(totals["agreements"], totals["joined"]),
            dict(joined=totals["joined"], agreements=totals["agreements"],
                 mismatches=totals["joined"] - totals["agreements"], first_20_mismatches=mismatches,
                 empty_after_clip=totals["empty_after_clip"], undefined_cases=undefined),
            "At least 99.9% exact POS/CIGAR agreement; zero joined reads cannot establish agreement")
    # Whole-input name membership prevents double counting long reads across islands.
    with pysam.AlignmentFile(trimmed, "rb") as handle:
        names = set()
        minimum_mapq = None
        for read in handle:
            if primary_mapped(read):
                names.add(read.query_name)
                minimum_mapq = read.mapping_quality if minimum_mapq is None else min(minimum_mapq, read.mapping_quality)
    no_trimmed = 0
    with pysam.AlignmentFile(bam, "rb") as handle:
        for read in handle:
            if primary_mapped(read) and read.query_name not in names:
                no_trimmed += 1
    log.add("F10 filtered reads", True, dict(original_without_trimmed=no_trimmed, minimum_trimmed_mapq=minimum_mapq),
            "Report whole-input original primary mapped reads without trimmed records and minimum trimmed MAPQ", "REPORT")
    log.add("F11 depth", True, dict(maximum_island_depth=depth, islands_above_200000=over_capacity),
            "Report maximum depth and all islands exceeding 200000", "REPORT")


def run_preflight(bam, trimmed_bam, reference, outdir):
    """Shared preflight API; writes partial evidence and stops on a HARD failure."""
    outdir = Path(outdir)
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        print("Output directory must be absent or empty", file=sys.stderr)
        return 2
    work = outdir / "work"
    work.mkdir(parents=True)
    log = CheckLog(work)
    try:
        for path in (bam, trimmed_bam, reference, str(reference) + ".fai"):
            with open(path, "rb"):
                pass
        check_headers_and_samples(bam, trimmed_bam, reference, log)
        original_index = ensure_index(bam, work, "original")
        trimmed_index = ensure_index(trimmed_bam, work, "trimmed")
        (work / "indexes.json").write_text(json.dumps(dict(original=original_index, trimmed=trimmed_index), indent=2) + "\n")
        check_islands(bam, trimmed_bam, original_index, trimmed_index, work, log)
    except HardFailure:
        return 2
    except (OSError, ValueError, pysam.SamtoolsError) as error:
        try:
            log.add("input or processing error", False, str(error), "Readable, valid inputs and usable indexes")
        except HardFailure:
            pass
        print(str(error), file=sys.stderr)
        return 2
    return 0