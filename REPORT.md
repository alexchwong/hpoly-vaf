# Gate report

Append to this file. Do not rewrite earlier entries.

## Kit authoring environment, 14 and 15 November 2026

Ubuntu 24.04, g++ 13.3, bamtools v2.5.3, Armadillo 15.6 headers, Boost 1.90 headers.
LAPACK came from the `scipy-openblas32` wheel through `kit/gen/lapack_shim.c`. Docker, conda and system `libopenblas` routes were not executed.

| Phase | Gate | Result |
|---|---|---|
| P0 | five binaries built from a clean work directory | PASS |
| P1 | G1 single allele `FAO=1200 FDP=6000 FRO=4800` | PASS |
| P1 | G2 three alleles `FAO=1200,0,0 FDP=6000 FRO=4800` | PASS |
| P1 | G3 rerun identical | PASS |
| P1 | G4 default sampling `FDP=2000` | PASS |
| P1 | G5 trace VCF body equals unmodified | PASS |
| P1 | G6 `groups.tsv` | PASS |
| P1 | G7 `reads.tsv` classes | PASS |
| P1 | G8 pruned VCF body equals unmodified | PASS |
| P1 | G9 pruned trace equals trace | PASS |
| P1 | G10 soft-clipped reads give the same trace | PASS |

## Operator, 15 November 2026

`kit/check_trimmed_bam.py` on an example BAM pair: `VERDICT: SOFT`. Full output not recorded here; paste it when P2 starts.

## Local machine

Rerun P0 and P1 and record results below. Phases P2 onward have not been run anywhere.

## Undefined cases
