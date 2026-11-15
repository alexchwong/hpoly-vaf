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

## Local setup and P0, 15 November 2026

Public GitHub repository: `alexchwong/hpoly-vaf`, branch `main`.
Initial commit: `8ad9cdb` (`kit as delivered`), authored by Alex Wong with
GitHub-linked no-reply email `80015046+alexchwong@users.noreply.github.com`.
Original delivery duplicates, local secrets/configuration, sequencing inputs,
build products and run outputs are ignored before the initial commit.
The kit remains byte-identical to the delivered archive, including its original
whitespace; no kit files were edited.

Work directory: `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder`.
Free space before environment setup: 1,369,147,961,344 bytes (requirement: at least 10 GB).

### Environment

Docker build exited 0 using the unchanged `kit/Dockerfile`.
Image tag: `hpoly-vaf-env:latest`.
Image ID: `sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43`.
Ubuntu 24.04; g++ 13.3.0; CMake 3.28.3; Python packages:
pysam 0.22.0, numpy 1.26.4, scipy 1.11.4.
System OpenBLAS was used without the compatibility shim.
Keep this image and toolchain unchanged for subsequent phases.

The repository was mounted read-only at `/repo` and the external work directory
read-write at `/work`. Containers ran with the host user's UID and GID.

### P0 input gate (HARD): PASS

Command inside the image:

```bash
for cmd in git cmake g++ make; do command -v "$cmd"; done
python3 -c "import pysam, numpy, scipy"
git ls-remote https://github.com/iontorrent/TS.git HEAD
```

Expected: all tools present, Python imports exit 0, GitHub reachable.
Observed: all tools present under `/usr/bin`; imports exit 0; GitHub returned
`5fe78dc2f22873845dbbecc7fc2ac8fc33ca42fb` for HEAD. PASS.

### P0 build and output gate (HARD): PASS

Build command inside the image: `/repo/kit/build_tvc.sh /work`.
Observed build exit status: 0. Default build parallelism: 4.

| Command/check inside container | Expected | Observed |
|---|---|---|
| `/work/bin/tvc 2>&1 \| head -1` | `tvc 5.18-6 (5fe78dc2) - Torrent Variant Caller` | Exact match; PASS |
| `ls /work/bin`, plus executable checks | `hpoly-eval`, `hpoly-eval-pruned`, `tvc`, `tvc-trace`, `tvcutils`, all executable | All five present and executable; PASS |
| `cmp /work/bin/hpoly-eval /work/bin/tvc-trace` | Identical | Exit 0; PASS |

Additional checks: TVC checkout is at the pinned commit above with no tracked
source changes; no binary has an unresolved shared library inside the image;
every kit file matches the delivered ZIP. All PASS.

| Binary | SHA-256 |
|---|---|
| `hpoly-eval` | `8a46e55f767b1824b38ca3dd1a7d3ff0e0297e9e5c5af19d635aede36bbae0ac` |
| `hpoly-eval-pruned` | `72f37a6772214c7ef50d91f2d568796a7b3cf28b3571773eded64170b1b50073` |
| `tvc` | `00bef98037aec53ed8e2fa09cb62723e7c7b95375bb848d06895e43a7654fb9d` |
| `tvc-trace` | `8a46e55f767b1824b38ca3dd1a7d3ff0e0297e9e5c5af19d635aede36bbae0ac` |
| `tvcutils` | `b84e0171a06caa6d4fabbc8642ff4ab02e63fb8ee946c1828938e9f57685c5fd` |

Logs are local, outside Git, under
`/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/`:
`docker-build.log`, `docker-build.exit`, `p0-input-gate.log`,
`p0-build.log`, `p0-build.exit`, and `p0-output-gate.log`.

Stopped after successful P0 at the operator's request. P1 has not run locally;
the synthetic output directory does not exist. The historical P1 PASS entries
above belong to the kit author, not this machine. See `docs/P1_PLAN.md`.

## Local P1, 15 November 2026

### Input gate (HARD): PASS

Reconfirmed the exact P0 version string, all five executable binaries, and
`cmp /work/bin/hpoly-eval /work/bin/tvc-trace` (exit 0). All five binary SHA-256
values match the P0 report. All kit files match the delivered archive. The
synthetic output path did not exist before this first local run.
Input-check log: `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-input-gate.log`;
recorded exit status: 0.

### Execution

Used the same P0 Docker image without rebuilding or changing the compiler:
`sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43`.

Exact host command (Bash; logs captured outside the repository):

```bash
docker run --rm --name hpoly-vaf-p1 \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v /home/alex-wong/vscode_projects/hpoly_vaf:/repo:ro \
  -v /home/alex-wong/vscode_projects/hpoly_vaf_workfolder:/work \
  -w /repo \
  sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43 \
  /repo/kit/run_synth_gates.sh /work \
  > /home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-gates.log 2>&1
rc=$?
printf '%s\n' "$rc" > /home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-gates.exit
```

### Output gate (HARD): PASS

Observed script exit status: **0**. Complete console output:

```text
PASS G1_single_allele
PASS G2_three_alleles
PASS G3_rerun_identical
PASS G4_default_downsample
PASS G5_trace_vcf_equals_unmodified
PASS G6_groups_tsv
PASS G7_reads_tsv_classes
PASS G8_pruned_vcf_equals_unmodified
PASS G9_pruned_trace_equals_trace
PASS G10_soft_clipped_reads_same_trace
```

| Gate | Expected | Observed |
|---|---|---|
| G1 | `FAO=1200 FDP=6000 FRO=4800` | Exact match; PASS |
| G2 | `FAO=1200,0,0 FDP=6000 FRO=4800` | Exact match; PASS |
| G3 | Repeat VCF body identical | Identical; PASS |
| G4 | Default sampling `FDP=2000` | Exact match; PASS |
| G5 | Trace VCF body equals unmodified | Identical; PASS |
| G6 | `chrT 300 CAA CA,C,CAAA 6000 0 4800 1200,0,0` | Exact tab-separated match; PASS |
| G7 | 4,800 class 0; 1,200 class 1 | Exact match; PASS |
| G8 | Pruned VCF body equals unmodified | Identical; PASS |
| G9 | Pruned read and group traces equal traced build | Byte-identical; PASS |
| G10 | Soft-clipped read trace equals unclipped | Byte-identical; PASS |

Post-run verification independently checked all ten gate identifiers and exit
status, VCF body equality, the group row, read class counts, pruned traces and
soft-clipped read trace. `hpoly-eval` remains byte-identical to `tvc-trace`.
No kit or TVC source files were edited.

Retained local evidence:

- `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-gates.log`
- `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-gates.exit`
- `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/synth/` (inputs, indexes,
  per-run engine logs, VCF outputs and traces; not published)

| Evidence | SHA-256 |
|---|---|
| `logs/p1-gates.log` | `efd9576e50700b8fcb7b8d04a95c1a63df5631c935cd496860b5e53a68eaf0f6` |
| `synth/t3/reads.tsv` | `28b5bcac1aba08be6c2e84fdbfc2a29a68addb39db1a0ab1ac653f94bb91802d` |
| `synth/t3/groups.tsv` | `15ee5c49992fbbe1201c1c0509a17e53bd94b9a57e1fbac7b4ded37e14b5aa37` |

Stopped after successful P1 at the operator's request. P2 has not started:
real input paths have not been supplied, the real trimmed-BAM checker has not
been run locally, and no Python implementation was added. The detailed P2
plan is in `docs/P2_PLAN.md`. Synthetic pruning equivalence does not authorize
adoption of the pruned evaluator; that still requires P7 real-data gates.
