# P2 plan: preflight on the real files

Status: completed locally on 15 November 2026 after explicit Genexus hg19/rCRS
reference preparation. All F1–F9 HARD checks passed; F10/F11 recorded. The prior
UCSC-reference F1 failure is retained in REPORT.md. Next: docs/P3_PLAN.md.
The original implementation plan follows for reference.

## 1. Required operator inputs and safety boundary

Obtain absolute paths for:

1. Original flow-bearing `merged.bam`.
2. Matching primer-trimmed `merged.bam.ptrim.bam`.
3. Matching reference FASTA, with `.fai`.
4. A new empty or absent preflight output directory outside the repository.

Use a fresh child directory of
`/home/alex-wong/vscode_projects/hpoly_vaf_workfolder` for P2 output. Keep P1's
synthetic evidence intact. Mount real data read-only and place all indexes,
logs and generated outputs in the external output/work directory. Do not copy
sequencing inputs into Git or scan unrelated directories to guess their paths.

Continue using the exact P0/P1 image:
`sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43`.
Do not rebuild the evaluators or change the compiler. Python dependencies are
limited to pysam, numpy, scipy and the standard library; pytest is test-only.

Because the repository is public, inspect any real-data report before pushing:
checker output and F9 examples may contain read names, regions and file paths.
Do not automatically publish real-data logs. If that evidence cannot be
published, agree on a private evidence location before proceeding; do not
silently omit the full checker output required by the specification.

## 2. Pass the P2 input gate before writing implementation

Confirm the three input files exist and are readable. Run the unchanged
checker at
`/home/alex-wong/vscode_projects/hpoly_vaf/kit/check_trimmed_bam.py`
against the supplied original and trimmed BAMs inside the fixed image.

Required: checker prints `VERDICT: SOFT` (normally exit 0). Capture its full
output and exit status and paste the full output into
`/home/alex-wong/vscode_projects/hpoly_vaf/REPORT.md`, as required by the P2
input gate. Historical operator results and P1 synthetic clipping are not
substitutes for this real-pair check.

Missing or stale BAM indexes must never be written beside real input BAMs.
The checker has no index-path option. Where usable adjacent indexes are absent,
prepare a view of the same input files using symlinks in the external work
directory, and create local indexes alongside those links before calling the
unchanged checker. Preserve the original source paths in the local record.
Do not modify the checker or input data to obtain a passing verdict.

On unreadable inputs, a non-SOFT verdict, or another HARD failure: record the
command, expected value and observed value, then stop without implementing P2.

## 3. Implement the minimal CLI and shared preflight components

Use the specification's module layout under
`/home/alex-wong/vscode_projects/hpoly_vaf/src/hpoly_vaf/`:

| Module | P2 responsibility |
|---|---|
| `__main__.py`, `cli.py` | Both module entry points; `preflight --bam --trimmed-bam --reference --outdir` |
| `preflight.py` | Checks F1–F11, index handling, structured output and exit status |
| `islands.py` | Exact streaming island algorithm from section 6.1 |
| `join.py` | Exact per-island join and diagnostics from section 6.2 |
| `cigar.py` | Exact `clip_to_span` algorithm from section 6.3 |

These modules will expose reusable functions so the later `run` command calls
the same preflight implementation. Do not add a pretend `run` command, package
installer, classifier, variant output, or later-phase functionality in P2.

Support both `python3 -m hpoly_vaf` and `python3 -m hpoly_vaf.cli` with
`PYTHONPATH` pointing to the repository's absolute `src` directory. Use argparse
and existing allowed libraries, not a new CLI framework.

Create `<outdir>/work/preflight.json`, one entry per check with exactly the
required fields: `name`, `status` (`PASS` or `FAIL`), `observed`, `expected`.
Represent gate kind separately from status. Return 0 only when every HARD
check passes; use nonzero failure status (2, consistent with preflight failure
in the later `run` command). REPORT failures must not change the HARD decision.
Do not implement a bypass switch or continue after a HARD gate failure.

Find BAM indexes beside each input using the specified `.bai`/`.csi` rule.
If missing or older than its BAM, rebuild under `<outdir>/work/`, using distinct
names for the two inputs and `index_filename` when reopening them. Never write
beside the input BAMs. Validate the FASTA `.fai`; do not silently generate one
beside a read-only input.

## 4. Implement and evaluate every output check

| Check | Kind | Required behavior |
|---|---|---|
| F1 contigs | HARD | Every SQ name and length in both BAMs matches an entry in the FASTA index |
| F2 sort order | HARD | Both headers say `SO:coordinate`; contig/position order is non-decreasing over the first 100,000 records |
| F3 read groups | HARD | Every original RG contains FO, KS and SM; exactly one distinct SM |
| F4 run ID | HARD | First 10,000 original primary mapped reads have a colon; name prefix equals RG ID prefix before its first dot |
| F5 flow tags | HARD | Those same original reads all carry ZM, ZP, ZF and RG |
| F6 trimmed SEQ | HARD | First 10,000 trimmed primary mapped reads have no hard clip and CIGAR query length equals SEQ length |
| F7 join | HARD | Every trimmed primary read has exactly one original match by name, contig, strand and SEQ; AR is contained in original reference span; zero violations over all islands |
| F8 unique names | HARD | No repeated primary read name within any island in either BAM |
| F9 clip agreement | HARD | Exact clipped POS and CIGAR match for at least 99.9% of joined reads; report counts and first 20 mismatches |
| F10 filtered reads | REPORT | Count original primary mapped reads without a trimmed record and minimum trimmed MAPQ |
| F11 depth | REPORT | Maximum island depth and all islands with any position above 200,000 |

Follow sections 6.1–6.3 literally:

- Mapped records only; exclude secondary (`0x100`), supplementary (`0x800`)
  and QC-fail (`0x200`) records for islands and joins. Do not invent duplicate
  or MAPQ filters.
- Merge same-contig spans when the next start is less than or equal to the
  current end. Write `<outdir>/work/islands.bed` in input order.
- Fetch original reads with the specified 200 bp island extension. Detect
  duplicate names before building dictionaries so overwrite cannot hide F8
  failures. Collect results in island order; workers use separate BAM handles.
- Preserve SEQ and implement clipping with the specified left/middle/right
  phases and no post-processing of `mid`. Count empty-after-clip reads.
- Use interval boundary events for depth, not pileup.
- Distinguish F10's whole-input count from the per-island `NO_AR` diagnostic;
  avoid counting an original read twice merely because it spans two islands.

If F6 fails, stop: F9 cannot be evaluated, and the operator must decide. Do not
manufacture a passing F9 entry or treat skipped checks as PASS. Undefined
cases must be excluded and documented under "Undefined cases" in the report,
per the specification; ambiguities affecting gates must not weaken a gate.

## 5. Tests before real-data acceptance

Add pytest tests under
`/home/alex-wong/vscode_projects/hpoly_vaf/tests/`, creating small BAM/FASTA
fixtures in temporary directories with pysam. Do not use real patient data.

Required coverage:

- Both CLI entry points, argument errors, empty-output enforcement, JSON shape,
  HARD failure exit status and REPORT-only behavior.
- F1 mismatch; F2 wrong header, decreasing positions and contig regression;
  missing F3 fields and multiple samples; F4 prefix errors; every missing F5
  tag; F6 hard clipping and query/SEQ length disagreement where representable.
- Island overlap, adjacency, gaps, contig transitions and excluded flags.
- F7 missing originals, mismatching strand/contig/SEQ and AR containment;
  F8 duplicates in either BAM without dictionary overwrite.
- F9 exact equality, below/equal/above 99.9% boundary and first-20 diagnostic
  cap; explicit undefined-case handling for zero joined reads.
- F10 no-trimmed counts, low MAPQ and no double counting across islands;
  F11 boundary depths 200,000 and 200,001 using interval-event fixtures.
- Missing and stale indexes rebuilt only in work; fresh indexes reused;
  original/trimmed index naming collisions avoided. Read-only inputs unchanged.
- The four exact section 7 clipping examples, plus boundary insertions,
  deletions, skips, hard clips and empty-after-clip behavior. Preserve the
  specified algorithm even in unusual cases; do not substitute a library trim.

Run pytest in the fixed image with bytecode and pytest caches directed to the
external work directory. Run both CLI `--help` entry points there too. Inspect
all new files and the final diff before acceptance.

## 6. Real-data acceptance and stop boundary

Once unit tests pass, run preflight on the supplied BAM pair and reference in
the same fixed image. Real data remain read-only. Preserve the command,
environment ID, exit status, JSON output and mismatch diagnostics locally.
Validate all F1–F9 HARD results and record F10/F11 REPORT observations. Append
results to the report without rewriting prior entries. A HARD failure stops
the phase immediately and is not repaired by loosening the check.

P2 is complete only when the input gate and every HARD output check pass on
the real files, REPORT observations are recorded, and tests pass. Review
real-data evidence for public-release suitability before any GitHub push.
Stop after P2; do not start P3 raw classification, evaluator execution on real
data, or adoption of the pruned evaluator.