"""Per-island original/trimmed join and clipping diagnostics."""

from collections import Counter

import pysam

from .cigar import MATCH, clip_to_span
from .islands import maximum_depth, primary_mapped


def join_island(task):
    """Worker with private handles; task/results are safe to send between processes."""
    bam, trimmed, bam_index, trimmed_index, island = task
    result = dict(joined=0, violations=0, duplicates=0, agreements=0,
                  empty_after_clip=0, no_ar=0, mismatches=[], examples=[],
                  undefined=Counter(), maximum_depth=0)
    with pysam.AlignmentFile(trimmed, "rb", index_filename=trimmed_index) as t:
        reads = [r for r in t.fetch(island.contig, island.start, island.end)
                 if primary_mapped(r)]
    names = Counter(r.query_name for r in reads)
    result["duplicates"] += sum(n - 1 for n in names.values())
    ar = {r.query_name: r for r in reads if names[r.query_name] == 1}
    result["maximum_depth"] = maximum_depth((r.reference_start, r.reference_end) for r in reads)
    with pysam.AlignmentFile(bam, "rb", index_filename=bam_index) as original:
        end = min(original.get_reference_length(island.contig), island.end + 200)
        originals = [r for r in original.fetch(island.contig, max(0, island.start - 200), end)
                     if primary_mapped(r)]
    original_names = Counter(r.query_name for r in originals)
    result["duplicates"] += sum(n - 1 for n in original_names.values())
    by_name = {r.query_name: r for r in originals if original_names[r.query_name] == 1}
    for read in originals:
        if read.query_name not in names and read.reference_start < island.end and read.reference_end > island.start:
            result["no_ar"] += 1
    for name, t in ar.items():
        o = by_name.get(name)
        valid = (o is not None and o.reference_name == t.reference_name
                 and o.is_reverse == t.is_reverse and o.query_sequence == t.query_sequence
                 and o.reference_start <= t.reference_start and t.reference_end <= o.reference_end)
        if not valid:
            result["violations"] += 1
            if len(result["examples"]) < 20:
                result["examples"].append(dict(name=name, reason="missing, nonunique, mismatching original or AR containment"))
            continue
        result["joined"] += 1
        try:
            pos, cigar = clip_to_span(o.reference_start, o.cigartuples, t.reference_start, t.reference_end)
        except ValueError as error:
            result["undefined"][str(error)] += 1
            pos, cigar = None, []
        if not any(op in MATCH for op, _ in cigar):
            result["empty_after_clip"] += 1
        if pos == t.reference_start and cigar == t.cigartuples:
            result["agreements"] += 1
        elif len(result["mismatches"]) < 20:
            result["mismatches"].append(dict(name=name, contig=island.contig,
                                            clipped_pos=pos, clipped_cigar=cigar,
                                            trimmed_pos=t.reference_start, trimmed_cigar=t.cigartuples))
    # Duplicated trimmed names are excluded, never overwritten into a passing join.
    result["violations"] += sum(n for n in names.values() if n > 1)
    result["undefined"] = dict(result["undefined"])
    return result