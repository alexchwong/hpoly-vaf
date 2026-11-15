# P3 plan: raw alignment classification and counting

Planned only. Do not start until the latest real P2 output has passed all HARD
checks. Retain the P0/P1 Docker image and explicit prepared Genexus reference;
do not silently substitute the original UCSC chrM sequence.

## Inputs and environment

- Repository: `/home/alex-wong/vscode_projects/hpoly_vaf`.
- Work root: `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder`.
- P2 evidence: `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/p2_preflight_rcrs/work/`.
- Prepared FASTA: `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/data/reference/hg19_genexus_rcrs.fa`.
- Image: `sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43`.

Use the original and trimmed BAMs already relocated under the external data
directory. Keep input mounts read-only, use new output directories for each
validation run, and keep real read classes/logs private. The reference manifest
records zero BAM M5 tags: reference sequence identity to the original alignment
reference is not established, despite explicit recipe and length compatibility.

## Implementation sequence

1. Extend `cli.py` with `run --stop-after raw`, input options from specification
   section 3.1 and `--threads` (default 4). Call the same shared preflight code
   first; do not reuse the failed original-reference output or bypass checks.
2. Reuse islands, original/trimmed joining and literal clipping from P2. Keep
   original spans for overlap selection, count NO_AR and EMPTY_AFTER_CLIP, and
   never overwrite duplicate names. Workers open separate file handles and
   results are collected in island order.
3. Build the section 6.4 clipped-read object: uppercase SEQ, match blocks,
   insertion positions, deletion/skip intervals and aligned span. Implement
   `query_index` by binary search; do not use pileup or decode flow tags here.
4. In `runs.py`, implement section 6.5 maximal reference runs and flank runs,
   clamped island windows, A/C/G/T boundaries and optional BED-overlap selection.
   Implement section 6.6 insertion/deletion seeding before classification.
5. In `classify.py`, implement every step of section 6.7 in order: AR containment,
   skips, flank anchors, anchor base matches, flank-base stripping, homopolymer
   core and signed length change. Return exactly OUTSIDE_AR, REF, ALT(change),
   OTHER or UNINFORMATIVE. Undefined cases are excluded and recorded.
6. In `rawcount.py` and `groups.py`, implement section 6.8 counts, strand splits,
   integer 3% qualification, deterministic group fill, unsupplied counts,
   group-capacity and 200000-depth rules. Preserve the distinction between
   OUTSIDE_AR, NO_AR, OTHER and UNINFORMATIVE.
7. Write islands.bed, raw_counts.tsv, read_classes/<contig>.tsv and the raw VCF
   with FLOW_STATUS=NOT_RUN. Use specified coordinates and REF/ALT padding;
   sort in FASTA contig order, then position. Do not create P4 evaluator inputs,
   run the engine, or pretend flow evaluation has occurred.

## Mandatory tests

Run all existing tests plus every section 7 case: 14 exact classifier examples,
four clipping examples, threshold boundaries, group fill, candidate construction,
allele mapping and the same-read arithmetic example. Implement pure helpers
needed by these required tests without initiating later-phase execution.

Also test binary-search boundaries, reference ambiguity, contig ends, island
boundaries, seeding at inclusive insertion limits, skip overlap, no flank anchor,
strand counts, deterministic fill ties, depth/capacity limits, missing values,
unique read-class keys and malformed/undefined input exclusion. Do not replace
an algorithm with a library routine that only appears equivalent.

## Real-data output gates

| Gate | Kind | Acceptance |
|---|---|---|
| R1 | HARD | `pytest tests/` passes, including every section 7 case |
| R2 | HARD | Every raw-count row satisfies DP = RO + AO + RAW_OTHER + RAW_UNINFORMATIVE |
| R3 | HARD | No duplicate read name for the same run in read_classes |
| R4 | HARD | pysam parses the VCF, REF matches FASTA and records are sorted |
| R5 | HARD | Fresh runs with threads 1 and 4 produce byte-identical raw_counts.tsv |
| R6 | REPORT | Record RO, AO and RAWAF at chr1:43818247 TC>T, chr11:533987 T>TG, chr7:148506170 C>CA; reference export percentages 5.997%, 13.144%, 3.299% are not acceptance thresholds |
| R7 | REPORT | Count runs having more than four qualifying ALTs, by homopolymer length |

Preserve commands, image/reference hashes, logs and gate evidence. Append to
REPORT.md without rewriting past failures. Stop at a HARD failure; do not relax
checks. Review all real-data evidence before public release. Stop after P3;
P4 evaluator-input preparation requires both P1 and P3 outputs to pass.