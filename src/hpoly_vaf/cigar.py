"""Literal implementation of specification v3 section 6.3.

CIGAR tuples use pysam's (numeric operation, length) representation.
"""

MATCH = {0, 7, 8}
REFERENCE = MATCH | {2, 3}
QUERY = MATCH | {1, 4}


def clip_to_span(pos, cigar, a, b):
    """Return clipped POS/CIGAR without changing SEQ or normalizing the middle."""
    if not pos <= a <= b:
        raise ValueError("clip span must start at or after POS and have a <= b")
    ops = list(cigar)
    if any(op not in REFERENCE | {1, 4, 5} or n <= 0 for op, n in ops):
        raise ValueError("unsupported or non-positive CIGAR operation")
    i, r, lc, rc = 0, pos, 0, 0
    mid = []
    while i < len(ops) and r <= a:
        op, n = ops[i]
        if op == 5:
            i += 1
            continue
        if op in {4, 1}:
            lc += n
            i += 1
            continue
        gap = a - r
        if gap == 0 and op != 2:
            break
        if op != 2:
            if n > gap:
                lc += gap
                ops[i] = (op, n - gap)
                r = a
                break
            lc += n
            r += n
            i += 1
        else:
            r += n
            i += 1
            if n > gap:
                break
    new_pos = r
    while i < len(ops) and r < b:
        op, n = ops[i]
        if op == 5:
            i += 1
            continue
        if op in {4, 1}:
            mid.append((op, n))
            i += 1
            continue
        gap = b - r
        if op != 2:
            if n > gap:
                mid.append((op, gap))
                ops[i] = (op, n - gap)
                r = b
                break
            mid.append((op, n))
            r += n
            i += 1
        else:
            r += n
            i += 1
            if n >= gap:
                break
            mid.append((op, n))
    for op, n in ops[i:]:
        if op not in {5, 2}:
            rc += n
    return new_pos, ([(4, lc)] if lc else []) + mid + ([(4, rc)] if rc else [])