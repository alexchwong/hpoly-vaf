"""Assessable islands and interval-event depth; never uses pileup."""

from collections import Counter
from dataclasses import dataclass

SKIP_FLAGS = 0x4 | 0x100 | 0x800 | 0x200


def primary_mapped(read):
    return not read.flag & SKIP_FLAGS


@dataclass(frozen=True)
class Island:
    contig: str
    start: int
    end: int


def iter_islands(reads):
    """Stream section 6.1's merging algorithm in coordinate order."""
    current = None
    for read in reads:
        if not primary_mapped(read):
            continue
        start, end = read.reference_start, read.reference_end
        if end is None or end <= start:
            raise ValueError("mapped read without a positive reference span")
        contig = read.reference_name
        if current and contig == current.contig and start <= current.end:
            current = Island(contig, current.start, max(current.end, end))
        else:
            if current:
                yield current
            current = Island(contig, start, end)
    if current:
        yield current


def maximum_depth(spans):
    """Maximum half-open reference-span depth (adjacent endpoints do not overlap)."""
    events = Counter()
    for start, end in spans:
        events[start] += 1
        events[end] -= 1
    depth = maximum = 0
    for position in sorted(events):
        depth += events[position]
        maximum = max(maximum, depth)
    return maximum