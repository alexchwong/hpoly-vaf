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

## Local P2 input gate, 15 November 2026 (local evidence; not for automatic publication)

Moved the supplied BAM directory and compressed reference to the external workfolder/data directory. Same-filesystem inode and size checks passed; details are in logs/p2-relocation.json. Retained hg19.fa.gz, decompressed hg19.fa (gzip exit 0), and indexed the derived FASTA. Both supplied BAM indexes are fresh; samtools quickcheck passed.

Image: `sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43`. Data mounted read-only at `/data`.

Command inside the image:

```bash
python3 /repo/kit/check_trimmed_bam.py --bam /data/Horizon_Myeloid_11-11-2026_Myeloid_NGS_12-11-2026_Result_758_1794706198_DF/merged.bam --trimmed-bam /data/Horizon_Myeloid_11-11-2026_Myeloid_NGS_12-11-2026_Result_758_1794706198_DF/merged.bam.ptrim.bam
```

Expected: readable BAM pair and FASTA, VERDICT: SOFT. Observed: all readable; checker exit 0 and SOFT. HARD input gate PASS. Full checker output:

```text
Sampled region            : chr1:36931134-36932371
Trimmed reads sampled     : 2000
  with hard clip (H)      : 0
  with soft clip (S)      : 2000
  CIGAR length != SEQ     : 0
  tags present            : ZM=0, ZP=0, ZF=2000, ZL=2000, RG=2000
Matched in original BAM   : 2000   (not found: 0; duplicate names in original: 0)
  original tags present   : ZM=2000, ZP=2000, ZF=2000, RG=2000
  SEQ identical           : 2000
  SEQ shorter (substring) : 0
  SEQ differs otherwise   : 0
  strand differs          : 0
  span inside original    : 2000
  span unchanged          : 0

Examples where POS or CIGAR changed (original -> trimmed):
  YJDCO:10517:01698
    original pos=36931619 cigar=179M5S seqlen=184
    trimmed  pos=36931634 cigar=15S164M5S seqlen=184
  YJDCO:10489:05904
    original pos=36931630 cigar=110M2I16M1D39M1D4M1D27M1I6M1I19M4S seqlen=229
    trimmed  pos=36931634 cigar=4S106M2I16M1D39M1D4M1D27M1I6M1I19M4S seqlen=229
  YJDCO:10389:05870
    original pos=36931630 cigar=166M1D32M1I27M1I9M32S seqlen=268
    trimmed  pos=36931634 cigar=4S162M1D32M1I27M1I9M32S seqlen=268
  YJDCO:10380:05620
    original pos=36931630 cigar=98M1I15M1I32M seqlen=147
    trimmed  pos=36931634 cigar=4S94M1I15M1I32M seqlen=147
  YJDCO:10350:01391
    original pos=36931619 cigar=177M1D75M seqlen=252
    trimmed  pos=36931634 cigar=15S162M1D75M seqlen=252

VERDICT: SOFT - the trimmed BAM keeps the full read sequence; trimming is soft clipping only.
```

## Local P2 output gate attempt, 15 November 2026: STOPPED at F1 (HARD FAIL)

Implemented the preflight module entry points, F1–F11 checks, islands, joining, literal clipping algorithm and synthetic tests. In the unchanged P0/P1 Docker image, `python3 -m pytest /repo/tests -q -o cache_dir=/work/test-cache --basetemp=/work/p2-test-tmp` passed all 47 tests; both module help entry points passed. An initial synthetic test expected normalized adjacent soft clips; its expectation was corrected to preserve the specification algorithm, then the full suite passed.

Real-data command inside the fixed image (repository and data mounted read-only):

```bash
PYTHONPATH=/repo/src python3 -m hpoly_vaf preflight --bam /data/Horizon_Myeloid_11-11-2026_Myeloid_NGS_12-11-2026_Result_758_1794706198_DF/merged.bam --trimmed-bam /data/Horizon_Myeloid_11-11-2026_Myeloid_NGS_12-11-2026_Result_758_1794706198_DF/merged.bam.ptrim.bam --reference /data/reference/hg19.fa --outdir /work/p2_preflight
```

Expected F1: every SQ name and length in both BAMs appears in FASTA .fai with the same length. Observed: chrM length 16569 in both original and trimmed BAM, but 16571 in the supplied FASTA. Exit status 2. Full structured output:

```json
[
  {
    "name": "F1 contigs",
    "status": "FAIL",
    "observed": [
      {
        "bam": "original",
        "contig": "chrM",
        "bam_length": 16569,
        "fasta_length": 16571
      },
      {
        "bam": "trimmed",
        "contig": "chrM",
        "bam_length": 16569,
        "fasta_length": 16571
      }
    ],
    "expected": "Every BAM SQ name/length matches FASTA .fai",
    "kind": "HARD"
  }
]
```

Stopped immediately at F1 per specification §0 rule 3. F2–F11 have not been evaluated on real data. P2 is not complete; P3 is blocked. No automatic retries, reference/header edits, or gate relaxations were performed. Kit and TVC code unchanged. All changes and real-data evidence remain local and have not been committed or pushed.

Required next action: obtain the exact matching reference FASTA and index for these BAM headers (including chrM length 16569), or supply a BAM pair aligned to the intended reference. Header-length compatibility alone does not establish sequence identity; prefer the original alignment reference. Preserve this failed attempt and use a fresh outdir for the next attempt after operator resolution.

Undefined cases in the implementation: empty original/trimmed samples fail the corresponding sample check rather than pass vacuously; zero joined reads cannot establish F9 and fail; unsupported CIGAR operations are excluded from agreement and recorded as undefined diagnostics. These policies were tested synthetically but not reached in the failed real-data attempt.

## Explicit Genexus reference integration and P2 retry, 15 November 2026

Operator authorized reference preparation from the supplied script and notes. Renamed make_reference.py to src/hpoly_vaf/genexus_reference.py and exposed prepare-genexus-reference. Integration passed 54 synthetic tests and a real reference build. Temporary hg19_chrM_fix.md deleted only after confirmation; retained explanation in docs/GENEXUS_REFERENCE.md. Kit unchanged; F1 unchanged.

Reference preparation used the fixed P0/P1 image, original and trimmed BAM headers, relocated hg19.fa.gz, UCSC md5sum.txt and NCBI NC_012920.1 FASTA. Output /work/data/reference/hg19_genexus_rcrs.fa has 25 uppercase contigs in BAM order, chrM length 16569 and fixed MD5 c68f52674c9fb33aef52dcf399755519; UCSC archive checksum verified. Provenance: data/reference/hg19_genexus_rcrs.fa.manifest.json in the external workfolder. Both BAMs have zero M5 tags: no claim of alignment-reference sequence identity is made. Original hg19 reference and prior failed attempt preserved.

P2 input gate rerun: same unchanged kit/check_trimmed_bam.py command and real BAM pair as the previous entry; fixed image, data read-only. Expected SOFT; observed exit 0, SOFT. Full output:

```text
Sampled region            : chr1:36931134-36932371
Trimmed reads sampled     : 2000
  with hard clip (H)      : 0
  with soft clip (S)      : 2000
  CIGAR length != SEQ     : 0
  tags present            : ZM=0, ZP=0, ZF=2000, ZL=2000, RG=2000
Matched in original BAM   : 2000   (not found: 0; duplicate names in original: 0)
  original tags present   : ZM=2000, ZP=2000, ZF=2000, RG=2000
  SEQ identical           : 2000
  SEQ shorter (substring) : 0
  SEQ differs otherwise   : 0
  strand differs          : 0
  span inside original    : 2000
  span unchanged          : 0

Examples where POS or CIGAR changed (original -> trimmed):
  YJDCO:10517:01698
    original pos=36931619 cigar=179M5S seqlen=184
    trimmed  pos=36931634 cigar=15S164M5S seqlen=184
  YJDCO:10489:05904
    original pos=36931630 cigar=110M2I16M1D39M1D4M1D27M1I6M1I19M4S seqlen=229
    trimmed  pos=36931634 cigar=4S106M2I16M1D39M1D4M1D27M1I6M1I19M4S seqlen=229
  YJDCO:10389:05870
    original pos=36931630 cigar=166M1D32M1I27M1I9M32S seqlen=268
    trimmed  pos=36931634 cigar=4S162M1D32M1I27M1I9M32S seqlen=268
  YJDCO:10380:05620
    original pos=36931630 cigar=98M1I15M1I32M seqlen=147
    trimmed  pos=36931634 cigar=4S94M1I15M1I32M seqlen=147
  YJDCO:10350:01391
    original pos=36931619 cigar=177M1D75M seqlen=252
    trimmed  pos=36931634 cigar=15S162M1D75M seqlen=252

VERDICT: SOFT - the trimmed BAM keeps the full read sequence; trimming is soft clipping only.
```

### P2 retry output gate: PASS

Command in fixed image: `PYTHONPATH=/repo/src python3 -m hpoly_vaf preflight --bam /data/Horizon_Myeloid_11-11-2026_Myeloid_NGS_12-11-2026_Result_758_1794706198_DF/merged.bam --trimmed-bam /data/Horizon_Myeloid_11-11-2026_Myeloid_NGS_12-11-2026_Result_758_1794706198_DF/merged.bam.ptrim.bam --reference /data/reference/hg19_genexus_rcrs.fa --outdir /work/p2_preflight_rcrs`. Exit 0. Full checks retained locally:

```json
[
  {
    "name": "F1 contigs",
    "status": "PASS",
    "observed": [],
    "expected": "Every BAM SQ name/length matches FASTA .fai",
    "kind": "HARD"
  },
  {
    "name": "F2 sort order",
    "status": "PASS",
    "observed": {
      "original": {
        "header_coordinate": true,
        "records": 100000,
        "regressions": 0
      },
      "trimmed": {
        "header_coordinate": true,
        "records": 100000,
        "regressions": 0
      }
    },
    "expected": "SO:coordinate and non-decreasing coordinates in first 100000 records of both BAMs",
    "kind": "HARD"
  },
  {
    "name": "F3 read groups",
    "status": "PASS",
    "observed": {
      "groups": 94,
      "samples": [
        "Horizon_Myeloid_11-11-2026"
      ],
      "missing": []
    },
    "expected": "Every RG has FO, KS, SM; exactly one distinct SM",
    "kind": "HARD"
  },
  {
    "name": "F4 run id",
    "status": "PASS",
    "observed": {
      "checked": 10000,
      "invalid": 0
    },
    "expected": "First 10000 primary mapped reads: name prefix equals RG prefix",
    "kind": "HARD"
  },
  {
    "name": "F5 flow tags",
    "status": "PASS",
    "observed": {
      "checked": 10000,
      "missing": {}
    },
    "expected": "Every sampled original read has ZM, ZP, ZF, RG",
    "kind": "HARD"
  },
  {
    "name": "F6 trimmed SEQ",
    "status": "PASS",
    "observed": {
      "checked": 10000,
      "hard_clips": 0,
      "query_length_mismatches": 0
    },
    "expected": "No H and CIGAR query length equals len(SEQ) in first 10000 trimmed primary mapped reads",
    "kind": "HARD"
  },
  {
    "name": "F7 join",
    "status": "PASS",
    "observed": {
      "islands": 415,
      "joined": 2764057,
      "violations": 0,
      "examples": []
    },
    "expected": "Every trimmed primary read has exactly one matching original; AR inside original span; zero violations",
    "kind": "HARD"
  },
  {
    "name": "F8 unique names",
    "status": "PASS",
    "observed": {
      "duplicate_occurrences": 0
    },
    "expected": "Zero duplicate primary names within any island in either BAM",
    "kind": "HARD"
  },
  {
    "name": "F9 clip agreement",
    "status": "PASS",
    "observed": {
      "joined": 2764057,
      "agreements": 2764057,
      "mismatches": 0,
      "first_20_mismatches": [],
      "empty_after_clip": 0,
      "undefined_cases": {}
    },
    "expected": "At least 99.9% exact POS/CIGAR agreement; zero joined reads cannot establish agreement",
    "kind": "HARD"
  },
  {
    "name": "F10 filtered reads",
    "status": "PASS",
    "observed": {
      "original_without_trimmed": 141633,
      "minimum_trimmed_mapq": 1
    },
    "expected": "Report whole-input original primary mapped reads without trimmed records and minimum trimmed MAPQ",
    "kind": "REPORT"
  },
  {
    "name": "F11 depth",
    "status": "PASS",
    "observed": {
      "maximum_island_depth": 24785,
      "islands_above_200000": []
    },
    "expected": "Report maximum depth and all islands exceeding 200000",
    "kind": "REPORT"
  }
]
```

All F1–F9 HARD checks passed. F10/F11 REPORT results recorded above. 415 islands, 2764057 joined reads, zero violations/duplicates/mismatches, 100% clip agreement, 141633 original reads without trimmed records, minimum trimmed MAPQ 1, maximum depth 24785, no islands above 200000. No undefined cases observed.

Evidence: external workfolder logs/p2-rcrs-preflight.log and .exit; p2_preflight_rcrs/work/preflight.json, islands.bed, indexes.json and join_diagnostics.json. Integration tests: 54 passed. P2 complete; P3 planned in docs/P3_PLAN.md, not started. Changes and identifiers remain local; no commit/push performed.
