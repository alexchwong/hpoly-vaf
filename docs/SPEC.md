# Homopolymer indel VAF by extracted TVC flow evaluation

Implementation specification, extraction series v3, 15 November 2026
Audience: a coding model working on a local Linux machine that holds the data
Supersedes extraction specs v1 and v2 and spec v6. This document is self-contained.

Changes from v2: trimmed-BAM hypothesis confirmed (§4.3); repository and environment setup written out (§3.3 to §3.6); pruning replaced by tested slices and a feature map (§8); kit gates extended to ten.

---

## 0. Rules for the implementer

1. Work through the phases in §5 in order. Phases P2 and P3 may be done before P0 and P1; nothing else may be reordered.
2. Each phase has an **input gate** and an **output gate**. Do not start a phase until its input gate holds. Do not declare a phase done until every check in its output gate passes.
3. Gates are of two kinds. A **HARD** gate that fails means: stop, change nothing else, and write the failure to `REPORT.md` with the command, the expected value and the observed value. A **REPORT** gate is recorded in `REPORT.md` and does not block.
4. Never weaken, skip or re-word a gate. Never add a command-line switch that bypasses one.
5. Do not add dependencies beyond those in §3.4.
6. Do not edit any file under the TVC source tree except where a phase explicitly says to.
7. Where this spec gives an algorithm, implement that algorithm. Do not substitute a library routine that "does the same thing", and do not use pileup.
8. If the spec does not cover a case, choose the behaviour that excludes the read or group from counts, record it in `REPORT.md` under "Undefined cases", and continue.

---

## 1. Goal and deliverables

Ion Torrent alignments and flow measurements can give different allele fractions for the same homopolymer indel. The tool measures how much the allele fraction changes when the **same reads** are classified first by alignment and then by TVC's flow-signal model.

The flow classification is TVC's own code. The Torrent Variant Caller evaluator is built from public source, given a small additive trace writer, and later pruned. It is never rewritten.

| Deliverable | Description |
|---|---|
| `hpoly-vaf` | Python command-line tool: `preflight`, `run`, `verify` |
| `hpoly-eval` | C++ evaluator built from TVC source. It is a copy of `tvc-trace` until phase P7 replaces it with `hpoly-eval-pruned` |
| `homopolymer_indels.vcf.gz` + `.tbi` | Output of `run` |
| `run_manifest.json` | Output of `run` |
| `REPORT.md` | Gate results for every phase |

Out of scope: comparison with the Genexus VCF, non-homopolymer variants, runs of length 1, molecular-tag libraries, clinical use.

---

## 2. Terms

| Term | Meaning |
|---|---|
| Original BAM | `merged.bam`. Flow-bearing, coordinate-sorted. Sole source of sequence, alignment and flow evidence |
| Trimmed BAM | `merged.bam.ptrim.bam`. Used only for each read's reference span |
| AR | Assessable region of a read: half-open reference interval `[ref_start, ref_end)` of its primary record in the trimmed BAM |
| Run | Maximal stretch of one base A/C/G/T in the reference, length ≥ 2. 0-based half-open `[s, e)`, base `b`, length `L = e − s` |
| Flank runs | The maximal single-base stretches immediately left `[ls, s)` and right `[e, re)` of a run. Length may be 1 |
| Window | `[ls, re)` |
| Allele key | `(contig, s, e, b, change)` where `change` is the signed length change, never 0 |
| Group | One run with REF and its submitted ALT alleles |
| Island | A maximal interval covered by trimmed-BAM primary reads with no gap |
| Stack | The reads the evaluator used at a group |

All internal coordinates are 0-based half-open. VCF `POS` is 1-based.

---

## 3. Inputs, outputs, environment

### 3.1 Inputs to `hpoly-vaf run`

| Option | Required | Meaning |
|---|---|---|
| `--bam` | yes | Original BAM. A `.bai` or `.csi` is looked for beside it |
| `--trimmed-bam` | yes | Trimmed BAM, same rule for its index |
| `--reference` | yes | FASTA with `.fai` |
| `--outdir` | yes | Created if absent. Must be empty or absent |
| `--eval` | yes | Path to the `hpoly-eval` executable |
| `--targets` | no | BED. If given, only runs overlapping a BED interval are reported |
| `--tvc-params` | no | TVC parameter JSON, passed to the evaluator unchanged |
| `--threads` | no | Worker processes for the raw side. Default 4 |

An index is stale if its modification time is older than its BAM. A stale or missing index is rebuilt into `<outdir>/work/` by copying nothing: open the BAM with `index_filename` pointing at the rebuilt index. Never write beside the input BAMs.

The reference is read uppercase. Any base other than A, C, G, T breaks a run.

### 3.2 Output directory layout

```
<outdir>/
  homopolymer_indels.vcf.gz
  homopolymer_indels.vcf.gz.tbi
  run_manifest.json
  work/
    preflight.json
    islands.bed
    raw_counts.tsv          one row per run and observed change
    read_classes/           raw classes per group, one file per contig
    eval.bam, eval.bam.bai  evaluator input
    candidates.vcf          evaluator input
    trace/reads.tsv         evaluator output
    trace/groups.tsv        evaluator output
    tvc_out/                evaluator native output
```

### 3.3 The repository zip

Unzip `hpoly-vaf.zip`. It is the starting state of the repository:

```
hpoly-vaf/
  README.md
  REPORT.md                     gate results; append to it, never rewrite earlier entries
  docs/SPEC.md                  this document
  kit/
    build_tvc.sh                builds all evaluator binaries (P0)
    run_synth_gates.sh          gates G1 to G10 (P1)
    check_trimmed_bam.py        quick test of the trimmed BAM (run before P2)
    make_synth.py               synthetic flow BAM generator used by the gates
    CMakeLists.txt              build of unmodified tvc and tvcutils
    hpoly_trace.patch           additive trace writer
    pruned/CMakeLists.txt       build of the pruned evaluator
    pruned/0*.patch             tested pruning slices, applied in name order after the trace patch
    gen/                        version stubs; lapack_shim.c (not needed on a normal machine)
    Dockerfile  environment.yml environment definitions (§3.4)
  src/hpoly_vaf/                empty; phases P2 to P6 fill it
  tests/                        empty; §7 fills it
```

First actions, in order:

1. `cd hpoly-vaf && git init && git add -A && git commit -m "kit as delivered"`.
2. Create the environment (§3.4).
3. Choose a work directory `W` outside the repository with at least 10 GB free. The TVC clone is large.
4. Run P0 and P1.

Files under `kit/` are not edited by the implementer. If a kit script fails, record the failure in `REPORT.md` and stop.

### 3.4 Environment

Two supported routes. **Use Docker when it is available.** The byte-identical gates in P5 and P7 depend on one fixed compiler, and the Docker image pins the toolchain the kit was written on (Ubuntu 24.04, g++ 13). Use conda only where Docker cannot be run, for example on a shared cluster without root.

Neither route has been executed by the spec author; the kit was built on Ubuntu 24.04 with the same packages installed another way. P0 and P1 are the test of the environment.

**Route A, Docker.**

```bash
docker build -t hpoly-vaf-env -f kit/Dockerfile kit
docker run --rm -it \
  -v "$PWD":/repo -v /path/to/data:/data:ro -v /path/to/work:/work \
  -w /repo hpoly-vaf-env
# inside the container:
kit/build_tvc.sh /work
kit/run_synth_gates.sh /work
```

Mount the data read-only. Everything the tool writes goes under `/work` or the `--outdir`.

**Route B, conda.**

```bash
conda env create -f kit/environment.yml
conda activate hpoly-vaf
kit/build_tvc.sh /path/to/work
kit/run_synth_gates.sh /path/to/work
```

The environment supplies the compiler (`gxx_linux-64` 13), `cmake`, `make`, `git`, `zlib`, `openblas`, `samtools`, and Python with `pysam`, `numpy`, `scipy`, `pytest`. The kit's CMake files add `$CONDA_PREFIX/include` and `$CONDA_PREFIX/lib` automatically when a conda environment is active. If the link step reports `cannot find -lopenblas`, the environment is not active.

**Route C, native Debian or Ubuntu** (equivalent to route A without the container):

```bash
sudo apt-get install build-essential cmake git zlib1g-dev libopenblas-dev \
     python3-pysam python3-numpy python3-scipy python3-pytest samtools tabix
```

Rules for all routes:

- Once P1 has passed in an environment, do not change the compiler or rebuild in a different environment. If the environment must change, rerun P0, P1 and P5.
- Python packages allowed in `src/`: `pysam`, `numpy`, `scipy`, and the standard library. `pytest` for tests only.
- `build_tvc.sh` is safe to rerun. It reuses existing clones and restores the TVC tree to its pinned state at the start and end.
- Environment variable `JOBS` sets build parallelism (default 4).

### 3.5 Installing the Python tool

`src/hpoly_vaf` is run in place; no packaging is required:

```bash
PYTHONPATH=src python3 -m hpoly_vaf.cli preflight ...
```

Create `src/hpoly_vaf/__main__.py` and `cli.py` so that `python3 -m hpoly_vaf` and `python3 -m hpoly_vaf.cli` both work. Wherever this spec writes `hpoly-vaf <sub-command>`, it means `PYTHONPATH=src python3 -m hpoly_vaf <sub-command>`.

### 3.6 Repository layout after implementation

```
src/hpoly_vaf/
  __main__.py  cli.py  preflight.py  islands.py  runs.py  cigar.py
  classify.py  rawcount.py  groups.py  evalio.py  join.py  vcfout.py  manifest.py
tests/          pytest unit tests (§7)
reference/      frozen evaluator output from P5
```

---

## 4. Facts

### 4.1 Tested on synthetic data, 14 and 15 November 2026

These were tested by building TVC and running it on synthetic flow data (Ubuntu 24.04, g++ 13.3). They are not assumptions. Do not re-derive them; P0 and P1 re-confirm them on your machine.

| # | Fact |
|---|---|
| T1 | Unmodified `tvc` 5.18.6 builds from commit `5fe78dc2f22873845dbbecc7fc2ac8fc33ca42fb` with the kit's `CMakeLists.txt`. |
| T2 | `tvc` runs with only a FASTA, a flow BAM and a candidate file. No target BED, parameter file, error-motif file or other resource is needed. |
| T3 | The candidate file is a plain VCF whose INFO carries `OID`, `OPOS`, `OREF`, `OALT`, `OMAPALT`. A hand-written file in that form is accepted (format in §6.9). |
| T4 | Candidates at one position are reported as one record. REF is the longest REF and ALTs are rewritten against it. Example: candidates `CA>C`, `CAA>C`, `C>CA` come back as `CAA` with ALTs `CA,C,CAAA`. |
| T5 | `FDP = FRO + ΣFAO` in every tested group, with one and with three ALTs. |
| T6 | Without `--downsample-to-coverage`, a 6,000-read site is evaluated on 2,000 reads. With `--downsample-to-coverage 200000` all 6,000 are used. 200000 is the largest value `tvc` accepts. |
| T7 | Two runs give identical VCF bodies. Only the `##fileUTCtime` header line differs. One thread and four threads give identical bodies. |
| T8 | Soft-clipping 15 bases from both ends of every read (changing only POS and CIGAR) leaves every count and every per-read responsibility unchanged. The flow model reads the full `SEQ` and tags regardless of clipping. |
| T9 | The trace patch in the kit changes nothing in the VCF body and writes per-read classes that sum to FRO and FAO. |
| T10 | The rich diagnostic (`--do-json-diagnostic on`) is about 240 KB per read. It is unusable at panel scale and is not used. |
| T12 | With consensus, VCF unification, `tvcutils`, indel assembly and JSON diagnostics removed (about 7,000 source lines no longer compiled), the evaluator gives a byte-identical VCF body and trace. |
| T11 | Synthetic flow signal that is inconsistent with the read's `ZP` phasing values produces erratic counts. The kit's generator therefore writes `ZP = 0,0,0`. |

### 4.2 Read from source, not yet tested on real data

Paths are under `Analysis/VariantCaller/` in the TVC tree.

| # | Fact | Source |
|---|---|---|
| S1 | A read's class is its most responsible hypothesis. Responsibility index 0 is the outlier, 1 is REF, 2+ are ALTs in output order. A read whose evaluation failed is an outlier. | `EnsembleEval/StackEngine.cpp` |
| S2 | Outliers are excluded from FDP. `FXX = outliers / (outliers + FDP)`. | `Filter/VariantAssist.cpp` |
| S3 | Classes come from a model fitted over all reads at the site. | `EnsembleEval/PosteriorInference.cpp` |
| S4 | The trimmed BAM contains only reads TVC did not filter, with `ZM`, `ZP`, `ZS` removed and POS/CIGAR soft-clipped to the amplicon insert. `SEQ` is not shortened. | `BAMWalkerEngine.cpp`, `TargetsManager.cpp` |
| S5 | Every read group needs `FO`. The text before the first `.` in the read-group `ID` must equal the text before the first `:` in the read names of that group. | `Bookkeeping/InputStructures.cpp`, `Reads/ExtendedReadInfo.cpp` |
| S6 | Reads need `ZM`, `ZP`, `ZF`. The read prefix is read-group `KS` followed by `ZK`, `ZT`, `ZE` when present. | `Reads/ExtendedReadInfo.cpp` |

P2 tests S4 to S6 on the real files.

### 4.3 Confirmed on real data

| # | Fact | Evidence |
|---|---|---|
| C1 | The trimmed BAM keeps the full read sequence. Trimming is soft clipping only. | `kit/check_trimmed_bam.py` run by the operator on an example BAM pair, 15 November 2026: verdict `SOFT` |

C1 closes the main risk behind check F6. F6 stays in preflight because it is cheap and applies to every new BAM pair.

---

## 5. Phases and gates

`W` is the kit work directory. `D` is the directory holding the real data. Commands are run from `repo/`.

### P0. Build unmodified TVC

**Input gate (HARD).** `git`, `cmake`, `g++`, `make` are on PATH. `python3 -c "import pysam, numpy, scipy"` exits 0. The machine can reach `github.com`.

**Work.** `kit/build_tvc.sh W`. Change nothing in the kit. If the link step cannot find `openblas`, install `libopenblas-dev` and rerun.

**Output gate (HARD).**

| Check | Expected |
|---|---|
| `W/bin/tvc 2>&1 \| head -1` | `tvc 5.18-6 (5fe78dc2) - Torrent Variant Caller` |
| `ls W/bin` | `hpoly-eval`, `hpoly-eval-pruned`, `tvc`, `tvc-trace`, `tvcutils`, all executable |
| `cmp W/bin/hpoly-eval W/bin/tvc-trace` | identical |

### P1. Confirm the engine on synthetic data

**Input gate (HARD).** P0 output gate.

**Work.** `kit/run_synth_gates.sh W`.

**Output gate (HARD).** Exit status 0 and ten `PASS` lines:

| Gate | Meaning | Expected |
|---|---|---|
| G1 | one ALT | `FAO=1200 FDP=6000 FRO=4800` |
| G2 | three ALTs | `FAO=1200,0,0 FDP=6000 FRO=4800` |
| G3 | rerun | identical VCF body |
| G4 | default sampling | `FDP=2000` |
| G5 | trace build | VCF body identical to unmodified |
| G6 | `groups.tsv` | `chrT 300 CAA CA,C,CAAA 6000 0 4800 1200,0,0` |
| G7 | `reads.tsv` classes | 4800 of class 0, 1200 of class 1 |
| G8 | pruned build | VCF body identical to unmodified |
| G9 | pruned build | trace identical to `tvc-trace` |
| G10 | soft-clipped reads | `reads.tsv` identical to unclipped |

### P2. Preflight on the real files

**Input gate (HARD).** Original BAM, trimmed BAM and FASTA exist and are readable. `python3 kit/check_trimmed_bam.py --bam <original> --trimmed-bam <trimmed>` prints `VERDICT: SOFT`. Paste its full output into `REPORT.md`.

**Work.** Implement `hpoly-vaf preflight --bam --trimmed-bam --reference --outdir`. It writes `work/preflight.json` with one entry per check: `name`, `status` (`PASS` or `FAIL`), `observed`, `expected`. Exit status is 0 only if every HARD check passes. `run` calls the same code first.

**Output gate.**

| Check | Kind | Rule |
|---|---|---|
| F1 contigs | HARD | Every `@SQ` name and length in both BAMs appears in the `.fai` with the same length |
| F2 sort order | HARD | Both BAMs are coordinate-sorted (header `SO:coordinate`, and positions non-decreasing over the first 100,000 records) |
| F3 read groups | HARD | Every `@RG` in the original BAM has `FO`, `KS` and `SM`. Exactly one distinct `SM` value |
| F4 run id | HARD | For the first 10,000 primary mapped reads of the original BAM: the name contains `:`, and the text before it equals the text before the first `.` of that read's `RG` id |
| F5 flow tags | HARD | The same 10,000 reads all carry `ZM`, `ZP`, `ZF` and `RG` |
| F6 trimmed SEQ | HARD | For the first 10,000 primary mapped reads of the trimmed BAM: the CIGAR contains no `H`, and its query length equals `len(SEQ)` |
| F7 join | HARD | Over every island (§6.1 to §6.2): each trimmed primary read has exactly one original primary read with the same name, contig, strand and `SEQ`, and `AR` lies inside the original's reference span. Zero violations |
| F8 unique names | HARD | No name occurs twice among primary records within an island, in either BAM |
| F9 clip agreement | HARD | For at least 99.9% of joined reads, `clip_to_span` (§6.3) applied to the original record gives exactly the trimmed record's POS and CIGAR. Report the count and the first 20 mismatches |
| F10 filtered reads | REPORT | Number of original primary mapped reads with no trimmed record, and the minimum MAPQ in the trimmed BAM |
| F11 depth | REPORT | Maximum island depth. If any position exceeds 200,000 reads, list those islands |

F6 is expected to pass (fact C1). If it fails on some BAM pair, that pair stores shortened reads. Stop: F9 cannot be evaluated and the operator must decide.

### P3. Raw side

**Input gate (HARD).** P2 output gate.

**Work.** Implement §6.1 to §6.8 and the unit tests in §7. `hpoly-vaf run --stop-after raw` writes `work/islands.bed`, `work/raw_counts.tsv`, `work/read_classes/` and a VCF containing only raw fields with `FLOW_STATUS=NOT_RUN`.

**Output gate.**

| Check | Kind | Rule |
|---|---|---|
| R1 unit tests | HARD | `pytest tests/` passes, including every case in §7 |
| R2 accounting | HARD | For every row of `raw_counts.tsv`: `DP = RO + AO + RAW_OTHER + RAW_UNINFORMATIVE` |
| R3 single count | HARD | No read name appears twice for the same run in `read_classes/` |
| R4 VCF validity | HARD | The VCF is parsed by `pysam.VariantFile`, every REF equals the FASTA at that position, and records are sorted |
| R5 determinism | HARD | Two runs with `--threads 1` and `--threads 4` give byte-identical `raw_counts.tsv` |
| R6 known sites | REPORT | For `chr1:43818247 TC>T`, `chr11:533987 T>TG`, `chr7:148506170 C>CA`: report `RO`, `AO`, `RAWAF`. Reference values from the Genexus export are 5.997%, 13.144% and 3.299% |
| R7 capacity | REPORT | Number of runs with more than four qualifying ALTs, by run length |

### P4. Evaluator inputs

**Input gate (HARD).** P1 and P3 output gates.

**Work.** Implement §6.9 and §6.10: write `work/eval.bam` (sorted, indexed) and `work/candidates.vcf`.

**Output gate (HARD).**

| Check | Rule |
|---|---|
| E1 tags intact | For 10,000 sampled reads, `SEQ`, `QUAL` and every tag in `eval.bam` are byte-identical to the original record |
| E2 record count | Records in `eval.bam` = joined reads counted in F7 |
| E3 candidates accepted | `W/bin/tvcutils prepare_hotspots -v work/candidates.vcf -r REF -o /tmp/chk.vcf` exits 0 and prints `Ignored 0 out of N lines` where N is the number of candidate lines. The tool merges alleles at one position into one multi-ALT line, so do not compare lines; check that the set of `OID` values in `/tmp/chk.vcf` equals the set in `candidates.vcf` |
| E4 sorted | `eval.bam` is coordinate-sorted and its index loads |

### P5. Run the evaluator and freeze the reference

**Input gate (HARD).** P4 output gate.

**Work.** Implement §6.11. Run `hpoly-eval` and unmodified `tvc` with the same arguments on the real inputs. Save `work/trace/` and both `sv.vcf` files. Copy `work/trace/reads.tsv`, `work/trace/groups.tsv` and the VCF body to `repo/reference/` and record their SHA-256 in `REPORT.md`. These three files are the **reference output**.

**Output gate (HARD).**

| Check | Rule |
|---|---|
| V1 trace is passive | VCF body from `hpoly-eval` equals VCF body from unmodified `tvc` (compare with `##` lines removed) |
| V2 identity | In every row of `groups.tsv`: stack size = outliers + REF count + sum of ALT counts |
| V3 trace vs VCF | For every group, REF count = `FRO` and each ALT count = the matching `FAO` in `sv.vcf` |
| V4 determinism | A second `hpoly-eval` run gives byte-identical `reads.tsv` and `groups.tsv` |
| V5 one row per group | Every submitted group appears in exactly one `groups.tsv` row. A group that is missing or split is listed; more than 1% missing or split fails the gate |
| V6 sensitivity | Copy `eval.bam`, and for one read at one group replace `ZM` with the `ZM` of a read of the opposite class at the same group. Rerun. `reads.tsv` must differ from the reference |

### P6. Join and final output

**Input gate (HARD).** P5 output gate.

**Work.** Implement §6.12 to §6.15. `hpoly-vaf run` now produces the final VCF and manifest. Implement `hpoly-vaf verify --run-dir OUT --tvc W/bin/tvc`, which reruns `work/eval.bam` and `work/candidates.vcf` through unmodified `tvc` and checks V1 and V3.

**Output gate (HARD).**

| Check | Rule |
|---|---|
| J1 flow accounting | Every row with `FLOW_STATUS=OK` has `FDP = FRO + ΣFAO` over its group and `FDP > 0` |
| J2 paired accounting | `PAIRED_FLOW_AO = PAIRED_RAW_AO − PAIRED_ALT_LOST + PAIRED_ALT_GAINED` on every row with `PAIRED_DP > 0` |
| J3 raw rows preserved | The set of (CHROM, POS, REF, ALT) equals that of the P3 VCF |
| J4 VCF validity | As R4, plus the `.tbi` loads and a region fetch returns records |
| J5 verify | `hpoly-vaf verify` exits 0 |
| J6 status coverage | Every row has `FLOW_STATUS` in `OK`, `NOT_EVALUABLE`, `ERROR`, with `FLOW_REASON` set unless `OK` |

### P7. Adopt the pruned evaluator

**Input gate (HARD).** P5 and P6 output gates. `W/bin/hpoly-eval-pruned` exists and P1 gates G8 and G9 passed.

**Work.** Run `hpoly-eval-pruned` on the real evaluator inputs with the §6.11 command into a fresh directory. Do not write new pruning code in this phase; §8 describes what the kit already removed and the rules for any further slice.

**Output gate (HARD).**

| Check | Rule |
|---|---|
| X1 synthetic | G8 and G9 pass |
| X2 reference | On the real inputs, `reads.tsv`, `groups.tsv` and the VCF body from `hpoly-eval-pruned` are byte-identical to `repo/reference/` |
| X3 adoption | Only if X1 and X2 pass: `cp W/bin/hpoly-eval-pruned W/bin/hpoly-eval`, then rerun the P6 output gate |

If X2 fails, keep `hpoly-eval` as the copy of `tvc-trace`, record the first differing line of each file in `REPORT.md`, and stop. The project is complete at P6; P7 changes only which binary is used.

---

## 6. Algorithms

### 6.1 Islands

Stream the trimmed BAM once. Use only records that are mapped and have none of flags `0x100`, `0x800`, `0x200`. Keep a current island `(contig, start, end)`. A read with `ref_start ≤ end` on the same contig extends `end` to `max(end, ref_end)`; otherwise close the island and start a new one. Write `work/islands.bed`. Every trimmed read lies wholly inside one island.

### 6.2 Join

For each island, in order:

1. Fetch trimmed primary reads overlapping the island into `ar[name] = (ref_start, ref_end, is_reverse, seq)`.
2. Fetch original primary reads overlapping the island extended by 200 bp on each side. Apply the same flag filter.
3. An original read with an entry in `ar` is **joined** if contig, strand and `SEQ` match and `AR` is inside its reference span. Otherwise it is a violation for F7.
4. An original read with no entry in `ar` is counted as `NO_AR` if it overlaps the island, and ignored otherwise.

Workers process islands in parallel. Each worker opens its own file handles. Results are collected in island order.

### 6.3 `clip_to_span(pos, cigar, a, b)`

Input: 0-based `pos`, CIGAR as a list of `(op, length)`, and `AR = [a, b)`. Output: new `pos` and CIGAR. `SEQ` is never changed.

Walk the CIGAR with a reference cursor `r = pos`. Maintain a left soft-clip count `lc`, a list `mid`, and a right soft-clip count `rc`.

- **Left phase**, while ops remain and `r ≤ a`:
  - `H`: skip.
  - `S` or `I`: add its length to `lc`; next op.
  - Let `gap = a − r`. If `gap == 0` and the op is not `D`: stop the left phase.
  - `M`, `=`, `X` or `N`: if length > `gap`, add `gap` to `lc`, reduce the op's length by `gap`, set `r = a`, stop the left phase. Otherwise add the whole length to `lc`, advance `r`, next op.
  - `D`: advance `r` by its full length and drop it. If its length was > `gap`, stop the left phase; otherwise next op.
- New `pos = r`.
- **Middle phase**, while ops remain and `r < b`:
  - `H`: skip.
  - `S` or `I`: append to `mid` unchanged; next op.
  - Let `gap = b − r`.
  - `M`, `=`, `X` or `N`: if length > `gap`, append `(op, gap)`, reduce the op by `gap`, set `r = b`, stop. Otherwise append it, advance `r`, next op.
  - `D`: if length ≥ `gap`, advance `r`, drop it, stop. Otherwise append it, advance `r`, next op.
- **Right phase**: every remaining op that is not `H` or `D` adds its length to `rc`.
- Result CIGAR: `(S, lc)` if `lc > 0`, then `mid`, then `(S, rc)` if `rc > 0`. Do not post-process `mid`.

This reproduces TVC's own trimming routine. Check F9 confirms it against the trimmed BAM.

A joined read whose clipped CIGAR has no `M`, `=` or `X` is dropped and counted as `EMPTY_AFTER_CLIP`.

### 6.4 Clipped read object

For each joined read build once:

- `name`, `is_reverse`, `seq` (as stored in the BAM, uppercase), clipped `pos` and CIGAR.
- `blocks`: list of `(ref_start, query_start, length)` for each `M`, `=`, `X` op of the clipped CIGAR.
- `ins`: list of reference positions at which an `I` op sits (the reference position of the next reference base).
- `dels`: list of `(ref_start, ref_end)` for `D` ops.
- `skips`: list of `(ref_start, ref_end)` for `N` ops.
- `span = (blocks[0].ref_start, last block ref_end)`.

`query_index(p)` returns the query index of reference position `p` if `p` falls in a block, else `None`. Use binary search over `blocks`.

Raw classification uses only these objects. Flow tags are not decoded on the raw side.

### 6.5 Runs

For each island, read the reference from `island_start − 100` to `island_end + 100`, clamped to the contig. Scan for maximal single-base stretches. A run is **in scope** when:

- its length is ≥ 2;
- it lies wholly inside the island;
- `ls ≥ 1` positions exist on each side, meaning `s > 0`, `e < contig length`, and the bases at `s − 1` and `e` are A, C, G or T;
- if `--targets` was given, `[s, e)` overlaps a BED interval.

Store `s, e, b, ls, re` for each.

### 6.6 Seeding

A run is **seeded** if at least one clipped read has:

- an insertion position `x` with `ls ≤ x ≤ re`, or
- a deletion `[d0, d1)` with `d0 < re` and `d1 > ls`.

Only seeded runs are classified and counted.

### 6.7 Classification of one read at one run

Return exactly one of `OUTSIDE_AR`, `REF`, `ALT(change)`, `OTHER`, `UNINFORMATIVE`.

1. If `span` does not contain `[s − 1, e + 1)`: return `OUTSIDE_AR`.
2. If any skip interval overlaps `[ls, re)`: return `UNINFORMATIVE`.
3. **Left anchor.** For `p` from `s − 1` down to `ls`: if `query_index(p)` is not `None`, set `pl = p`, `ql = query_index(p)` and stop. If none is found: return `UNINFORMATIVE`.
4. **Right anchor.** For `p` from `e` up to `re − 1`: likewise giving `pr`, `qr`. If none: return `UNINFORMATIVE`.
5. If `seq[ql]` differs from the reference base at `pl`, or `seq[qr]` differs from the reference base at `pr`: return `OTHER`.
6. Let `mid = seq[ql + 1 : qr]`. Let `lb` be the reference base at `s − 1` and `rb` the reference base at `e`.
7. Strip all leading `lb` from `mid`. Strip all trailing `rb` from what remains. Call the result `core`.
8. If `core` contains any base other than `b`: return `OTHER`.
9. Let `n = len(core)`. If `n == L`: return `REF`. Otherwise return `ALT(n − L)`.

Note for step 7: `lb ≠ b` and `rb ≠ b` by construction, so the stripping is unambiguous. Extra or missing flank bases are ignored; they belong to the neighbouring run.

### 6.8 Counts, threshold, groups

For each seeded run, classify every joined read whose original reference span overlaps `[s, e)`.

- `RAW_OUTSIDE_AR` = reads returning `OUTSIDE_AR`.
- `RAW_NO_AR` = `NO_AR` reads from §6.2 whose reference span overlaps `[s, e)`.
- `DP` = all other classified reads.
- `RO` = `REF` reads. `AO[c]` = reads returning `ALT(c)`. `OTHER_SEQ` = `OTHER` reads. `RAW_UNINFORMATIVE` = `UNINFORMATIVE` reads.
- For allele `c`: `RAW_OTHER = OTHER_SEQ + Σ AO[c'] for c' ≠ c`.
- Strand counts `SRF, SRR` split `RO` and `SAF, SAR` split `AO[c]` by `is_reverse`.
- `RAWAF = AO[c] / (AO[c] + RO)`.

**Qualifying.** Allele `c` qualifies when `AO[c] > 0` and `100 × AO[c] ≥ 3 × (AO[c] + RO)`, in integer arithmetic.

**Group.** For a run with at least one qualifying allele:

1. `Q` = qualifying changes. If `len(Q) > 4`: no group is submitted; every row of the run gets `FLOW_STATUS=NOT_EVALUABLE`, `FLOW_REASON=GROUP_CAPACITY`.
2. `N_simple = RO + Σ AO[c]` over all observed changes. `floor = max(2, ceil(0.005 × N_simple))`.
3. Fill candidates: non-qualifying changes with `AO[c] ≥ floor`, sorted by `AO` descending, then `c` ascending. Add them until the group has 4 ALTs.
4. Submitted set `S = Q ∪ fill`. `RAW_UNSUPPLIED = Σ AO[c] for c ∉ S`. `RAWAF_JOINT[c] = AO[c] / (RO + Σ AO[c'] for c' ∈ S)`.
5. If `DP > 200000`: `FLOW_STATUS=NOT_EVALUABLE`, `FLOW_REASON=DEPTH_LIMIT`, no group submitted.

Group id is `<contig>:<s>-<e>`. Write the raw class of every read at every submitted group to `work/read_classes/<contig>.tsv` as `group_id, read_name, class` where class is `REF`, the signed change, `OTHER` or `UNINFORMATIVE`.

### 6.9 Candidate file

For a group at run `[s, e)` with base `b`, let `pad` be the reference base at `s − 1`. For each submitted change `c`:

- deletion (`c < 0`, `k = −c`): `REF = pad + b×k`, `ALT = pad`
- insertion (`c > 0`, `k = c`): `REF = pad`, `ALT = pad + b×k`
- `POS = s` (the 1-based position of `pad`)

File format, tab-separated, sorted by contig order of the FASTA then POS:

```
##fileformat=VCFv4.1
##allowBlockSubstitutions=true
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	12345	.	CA	C	.	.	OID=<group_id>.<c>;OPOS=12345;OREF=CA;OALT=C;OMAPALT=C
```

One line per allele. `OPOS`, `OREF`, `OALT` repeat POS, REF, ALT. `OMAPALT` repeats ALT.

### 6.10 Evaluation BAM

For each joined read not dropped in §6.3, write the original record with only `reference_start` and `cigartuples` replaced by the clipped values. Remove the `MD` and `NM` tags. Change nothing else. Use the original BAM's header unchanged. Write unsorted, then `pysam.sort` to `work/eval.bam` and `pysam.index`.

### 6.11 Running the evaluator

`work/trace/` must exist and be empty: the evaluator appends.

```bash
HPOLY_TRACE_DIR=<outdir>/work/trace <eval> \
  --reference <fasta> --input-bam <outdir>/work/eval.bam \
  --input-vcf <outdir>/work/candidates.vcf \
  --process-input-positions-only on --use-input-allele-only on \
  --output-dir <outdir>/work/tvc_out --output-vcf sv.vcf \
  --num-threads 1 --downsample-to-coverage 200000 \
  --min-mapping-qv 0 --read-snp-limit 1000 \
  --read-max-mismatch-fraction 1.0 --min-cov-fraction 0
```

Append `--parameters-file <json>` only when `--tvc-params` was given. Do not pass any other option. A non-zero exit status is `FLOW_STATUS=ERROR`, `FLOW_REASON=ENGINE_EXIT` for every submitted group, and `run` exits 3 after writing outputs.

**`groups.tsv` columns:** `chrom, pos, ref, alts, stack, outliers, ref_count, alt_counts`. `alts` and `alt_counts` are comma-separated and parallel.

**`reads.tsv` columns:** `chrom, pos, ref, alts, read_name, is_reverse, success, class, read_count, responsibilities`. `class` is −1 for outlier, 0 for REF, `k ≥ 1` for the k-th entry of `alts`.

### 6.12 Mapping evaluator output to allele keys

For each `groups.tsv` row with `pos`, `ref`, `alts`:

1. Find the group whose run has `s == pos`. If none: record `UNMAPPED_ROW` in the manifest and skip the row.
2. For each ALT string `a`: `change = len(a) − len(ref)`. Check that `ref[0] == a[0] == pad` and that `ref[1:]` and `a[1:]` consist only of `b`. If the check fails or `change == 0`: the whole group gets `FLOW_STATUS=ERROR`, `FLOW_REASON=ALLELE_MAP`.
3. The set of changes must equal the submitted set `S`. Otherwise `ERROR / ALLELE_MAP`.
4. If a group has no row: `NOT_EVALUABLE / NO_ENGINE_ROW`. If it has more than one row: `ERROR / GROUP_SPLIT`.
5. If `stack − outliers == 0`: `NOT_EVALUABLE / NO_FLOW_DEPTH`.
6. Otherwise `FLOW_STATUS=OK`.

### 6.13 Flow fields

For an `OK` group: `FDP = ref_count + Σ alt_counts`, `FRO = ref_count`, `FAO[c]` = the matching alt count, `FXX = outliers / stack`. Strand counts come from `reads.tsv`: `FSRF, FSRR` for class 0, `FSAF, FSAR` for the allele's class, split by `is_reverse`. `AF = FAO / FDP`. `AF_REFALT = FAO / (FAO + FRO)`, missing when the denominator is 0.

### 6.14 Same-read join

For an `OK` group, join `work/read_classes` to `reads.tsv` rows of that group by read name. A read's flow class is `REF`, a change, `OUTLIER`, or `NOT_IN_STACK` when it has no row.

The **paired set** is the reads whose raw class is `REF` or a change in `S` **and** whose flow class is `REF` or a change in `S`.

For allele `c`:

- `PAIRED_DP` = size of the paired set (same for every allele of the group).
- `PAIRED_RAW_AO` = paired reads with raw class `c`. `PAIRED_FLOW_AO` = paired reads with flow class `c`.
- `PAIRED_RAW_AF`, `PAIRED_FLOW_AF` = those divided by `PAIRED_DP`.
- `PAIRED_ALT_LOST` = paired reads with raw `c` and flow not `c`. `PAIRED_ALT_GAINED` = paired reads with raw not `c` and flow `c`.
- `FLOW_SHIFT_PP = 100 × (gained − lost) / PAIRED_DP`.
- `PAIRED_ALT_OUTLIER` = reads with raw class `c` and flow class `OUTLIER`.
- `PAIRED_ALT_NOT_IN_STACK` = reads with raw class `c` and flow class `NOT_IN_STACK`.

Interval: with `n = PAIRED_DP`, compute two-sided 97.5% Clopper–Pearson intervals `[gl, gh]` for `gained/n` and `[ll, lh]` for `lost/n` using `scipy.stats.beta.ppf` with tail probability 0.0125 each side (lower bound 0 when the count is 0, upper bound 1 when the count is `n`). `FLOW_SHIFT_CI_LO_PP = 100 × (gl − lh)`, `FLOW_SHIFT_CI_HI_PP = 100 × (gh − ll)`.

When `PAIRED_DP == 0`, all paired fractions, the shift and the interval are missing.

### 6.15 Output VCF and manifest

One biallelic record per qualifying allele. Sort by FASTA contig order, then POS, REF, ALT. `ID`, `QUAL`, `FILTER` are `.`. REF and ALT are built as in §6.9. One sample column named by the read-group `SM`. Write with `pysam.VariantFile` in `wz` mode, then `pysam.tabix_index(preset="vcf")`.

Header meta lines: `##source=hpoly-vaf`, `##tvc_source_commit=5fe78dc2f22873845dbbecc7fc2ac8fc33ca42fb`, `##tvc_version=5.18-6`.

INFO fields:

| ID | Number | Type | Content |
|---|---|---|---|
| `TYPE` | 1 | String | `ins` or `del` |
| `LEN` | 1 | Integer | absolute change |
| `HPBASE` | 1 | String | `b` |
| `HPLEN` | 1 | Integer | `L` |
| `HPSTART` | 1 | Integer | `s + 1` |
| `HPEND` | 1 | Integer | `e` |
| `HPCHANGE` | 1 | Integer | signed change |
| `FLOW_GROUP` | 1 | String | group id |
| `FLOW_NALT` | 1 | Integer | size of `S`, 0 if not submitted |
| `FLOW_STATUS` | 1 | String | `OK`, `NOT_EVALUABLE`, `ERROR`, or `NOT_RUN` |
| `FLOW_REASON` | 1 | String | reason code; absent when `OK` |

FORMAT fields, all `Number=1` except `AD` (`Number=R`):

| Type | IDs |
|---|---|
| Integer | `DP, RO, AO, SRF, SRR, SAF, SAR, AD, RAW_OTHER, RAW_UNINFORMATIVE, RAW_OUTSIDE_AR, RAW_NO_AR, RAW_UNSUPPLIED, FDP, FRO, FAO, FSRF, FSRR, FSAF, FSAR, PAIRED_DP, PAIRED_RAW_AO, PAIRED_FLOW_AO, PAIRED_ALT_LOST, PAIRED_ALT_GAINED, PAIRED_ALT_OUTLIER, PAIRED_ALT_NOT_IN_STACK` |
| Float | `RAWAF, RAWAF_JOINT, AF, AF_REFALT, FXX, PAIRED_RAW_AF, PAIRED_FLOW_AF, FLOW_SHIFT_PP, FLOW_SHIFT_CI_LO_PP, FLOW_SHIFT_CI_HI_PP` |

`AD = [RO, AO]`. Fractions are written as fractions, not percentages, to 6 significant figures. Fields ending `_PP` are percentage points. A field that does not apply is omitted from the record's sample (missing), never written as 0.

`run_manifest.json` keys: `inputs` (path, size, SHA-256 of each input and of `eval.bam`), `eval_binary` (path, SHA-256), `eval_command`, `preflight` (copy of `preflight.json`), `read_counts` (`joined`, `no_ar`, `empty_after_clip`), `runs` (`in_scope`, `seeded`, `with_qualifying`), `groups_by_status`, `over_capacity_by_hplen`, `unmapped_rows`, `undefined_cases`.

Exit status of `run`: 0 on success, 2 on a failed preflight, 3 when the engine failed or any group has `FLOW_STATUS=ERROR`.

---

## 7. Unit tests (gate R1)

Reference `ACGTCAAAAGTCA`. The run is `AAAA` at `[5, 9)`, left flank `C` at 4, right flank `G` at 9. Each test builds a clipped read object directly and calls the classifier.

| # | Read (pos, CIGAR, SEQ) | Expected |
|---|---|---|
| 1 | 0, `13M`, `ACGTCAAAAGTCA` | `REF` |
| 2 | 0, `5M1D7M`, `ACGTCAAAGTCA` | `ALT(−1)` |
| 3 | 0, `8M1D4M`, `ACGTCAAAGTCA` | `ALT(−1)` |
| 4 | 0, `5M1I8M`, `ACGTCAAAAAGTCA` | `ALT(+1)` |
| 5 | 0, `9M1I4M`, `ACGTCAAAAAGTCA` | `ALT(+1)` |
| 6 | 0, `5M4D4M`, `ACGTCGTCA` | `ALT(−4)` |
| 7 | 0, `13M`, `ACGTCAAGAGTCA` | `OTHER` |
| 8 | 0, `13M`, `ACGTTAAAAGTCA` | `OTHER` (left anchor mismatch) |
| 9 | 6, `7M`, `AAAGTCA` | `OUTSIDE_AR` |
| 10 | 0, `4M3N6M`, `ACGTAAGTCA` | `UNINFORMATIVE` |
| 11 | 0, `8M2D3M`, `ACGTCAAATCA` | `UNINFORMATIVE` (right flank has no aligned base) |

Reference `GGAAAACCCTT`. Runs `AAAA` at `[2, 6)` and `CCC` at `[6, 9)`.

| # | Read | Run | Expected |
|---|---|---|---|
| 12 | 0, `6M1D4M`, `GGAAAACCTT` | `AAAA` | `REF` |
| 13 | same | `CCC` | `ALT(−1)` |
| 14 | 0, `8M1D2M`, `GGAAAACCTT` | `CCC` | `ALT(−1)` |

Further required tests:

- `clip_to_span`: `pos=100`, `20M`, span `[105, 115)` → `pos=105`, `5S10M5S`.
- `clip_to_span`: `pos=100`, `5M2I13M`, span `[105, 115)` → `pos=105`, `7S10M3S`.
- `clip_to_span`: `pos=100`, `5M3D12M`, span `[105, 115)` → `pos=108`, `5S7M5S`.
- `clip_to_span`: `pos=100`, `3S17M`, span `[100, 117)` → unchanged.
- Threshold: `AO=3, RO=97` qualifies; `AO=3, RO=98` does not; `AO=0` never qualifies.
- Group fill: `RO=1000`, `AO={−1: 200, +1: 6, −2: 5, +2: 2}` → `Q={−1}`, `floor=7`, no fill, `RAW_UNSUPPLIED=13`.
- Candidate lines for run `AAAA` with pad `C`, changes −1 and +2: `CA>C` and `C>CAA` at `POS = s`.
- Allele mapping: ref `CAA`, alts `CA,C,CAAA`, pad `C`, base `A` → changes −1, −2, +1.
- Same-read join on the 3,800-read example: raw A6 570, flow A6 520, lost 80, gained 30 → `FLOW_SHIFT_PP = −1.31579`.

---

## 8. Pruning: what was cut, what is kept, and rules for going further

### 8.1 Principle

The evaluator is TVC's code with features removed, never rewritten. A feature may be removed only when the trace and VCF body stay byte-identical on both the synthetic gates and the real reference (P5). Code that is linked but never executed under the fixed settings is left in place: removing it cannot change a result and only adds risk.

### 8.2 Feature map

Coverage figures are the share of lines executed on the synthetic gate run with the §6.11 settings. They show what the fixed settings switch off. They are from synthetic data; real data may execute more.

**Cut. Done in `kit/pruned/`, tested byte-identical (gates G8, G9).**

| Feature | What was removed | Patch |
|---|---|---|
| `tvc consensus` sub-command | dispatch in `main`; `Consensus/Consensus.cpp`, `FlowSpaceConsensus.cpp`, `ConsensusParameters.cpp` | `01` |
| VCF unification and hotspot preparation inside `tvc` | `tvcutils/unify_vcf.cpp`, `tvcutils/prepare_hotspots.cpp` and their include | `01` |
| `tvcutils` executable | build target | `pruned/CMakeLists.txt` |
| Long-indel assembly | construction and calls in `VariantCaller.cpp` and `BAMWalkerEngine.cpp`; `IndelAssembly/IndelAssembly.cpp` | `02` |
| JSON diagnostics | two call blocks in `HandleVariant.cpp`; `EnsembleEval/DiagnosticJSON.cpp` | `03` |

**Keep. This is the method.** Do not edit these beyond the trace patch.

| Role | Files | Coverage |
|---|---|---|
| Site-level model, responsibilities, hard classification | `EnsembleEval/StackEngine.cpp`, `ShortStack.cpp`, `CrossHypotheses.cpp`, `PosteriorInference.cpp`, `BiasGenerator.cpp`, `SigmaGenerator.cpp` | 33 to 91% |
| Hypothesis construction | `HypothesisEvaluator.cpp`, `Splice/SpliceVariantHypotheses.cpp`, `Splice/ClassifyVariant.cpp`, `Splice/LocalContext.cpp` | 35 to 50% |
| Flow-signal prediction | `BaseCaller/TreephaserSSE.cpp` | 40% |
| Read unpacking, prefix and start flow | `Reads/ExtendedReadInfo.cpp`, `Bookkeeping/InputStructures.cpp` | 42% |
| Per-variant driver and trace | `HandleVariant.cpp` | 48% |
| Read stack, read filters, positions | `BAMWalkerEngine.cpp`, `TargetsManager.cpp`, `SampleManager.cpp`, `VariantCaller.cpp` | 21 to 76% |
| Candidate input and read alleles | `HotspotReader.cpp`, `external/freebayes/src/AlleleParser.cpp`, `external/vcflib/Variant.cpp` | 14 to 48% |
| Counts and VCF fields used by `verify` | `Filter/VariantAssist.cpp`, `Filter/DecisionTreeData.cpp`, `Bookkeeping/VcfFormat.cpp` | 26 to 61% |
| Parameters | `Bookkeeping/ExtendParameters.cpp`, `Util/OptArgs.cpp`, `AnalysisOrg/IO/OptBase.cpp`, jsoncpp | 14 to 84% |

**Leave in place, switched off by the fixed settings.** Each except the skew term was tested by removing its source files from the build. Every one left unresolved references from kept code, so removal means editing the method's files.

| Feature | Files | Coverage | Why it stays |
|---|---|---|---|
| Molecular-tag families | `MolecularTag.cpp`, `BaseCaller/MolecularTagTrimmer.cpp` | 2%, 13% | Family structures are arguments of the evaluator's entry function |
| Flow recalibration | `Calibration/LinearCalibrationModel.cpp`, `Calibration/FlowAlignment.cpp` | 4%, 0% | Initialised from kept code; off by default (`suppress-recalibration`) |
| Reference-solver phasing | `BaseCaller/DPTreephaser.cpp`, `BaseCaller/PIDloop.cpp` | 2%, 0% | Referenced from kept code |
| Read realignment | `realignment/Realigner.cpp` | 4% | Constructed in kept code |
| Sequence-specific error motifs | `Splice/ErrorMotifs.cpp` | 15% | Queried during allele classification; empty table without a motif file |
| Skew term | `EnsembleEval/SkewGenerator.cpp` | 0% | Member of the model |
| Smith-Waterman, tabix | `external/vcflib/smithwaterman/*`, `external/vcflib/tabixpp/*` | 0% | Referenced from kept code |
| Metrics, ordered BAM writer | `MetricsManager.cpp`, `OrderedBAMWriter.cpp` | 90%, 97% | Small, executed, harmless |

### 8.3 Settings that keep those features off

These are the §6.11 options. Do not change them: `--process-input-positions-only on`, `--use-input-allele-only on` (no de novo candidates), `--downsample-to-coverage 200000` (no sampling), permissive read filters, no `--target-file` (no primer trimming), no `--error-motifs`, no `--model-file`. Recalibration is suppressed by default. Molecular tags are off when the BAM header declares none.

### 8.4 Rules for any further slice

Further pruning is optional and is not needed for any deliverable. If attempted:

1. Start from the TVC tree with the trace patch and the three kit patches applied, committed in git.
2. One feature per commit. Delete code and the call sites of deleted code. Never change a statement that remains.
3. Never touch a file in the **Keep** table except to delete a call into the feature being removed.
4. Before deleting anything in the **Leave in place** table, build with `--coverage`, run the real evaluator inputs, and confirm with `gcov` that the lines to be deleted have zero executions on real data.
5. After each commit: rebuild, then require gates G8 and G9 and gate X2. On failure run `git reset --hard HEAD~1` and record the slice as `SKIPPED`. Do not repair.
6. Stop after two consecutive skipped slices.

Candidates in order of safety: the `--postprocessed-bam` writing path in `BAMWalkerEngine.cpp`; `depth.txt` output; `tvc_metrics.json` output; the black-list VCF. Molecular tags and recalibration should not be attempted.

---

## 9. What this spec does not establish

- Results are those of TVC 5.18.6 with default parameters unless `--tvc-params` is given. The Genexus export came from 5.20-17 with a vendor parameter set and 2,000-read sampling. Agreement with it is not expected and is not tested here.
- The facts in §4.2 come from reading source. P2 and P5 test them on the real files; until then they are unconfirmed for this data.
- The kit build and all ten gates were run on Ubuntu 24.04 with OpenBLAS supplied through a compatibility shim. The Docker and conda routes in §3.4, and linking against a system or conda `libopenblas`, have not been executed.
- Fact C1 rests on one example BAM pair. Preflight repeats the check on every run.
- Pruning equality (T12) is shown on synthetic data only until gate X2 passes on real data.
- TVC source is under GPL version 2. `hpoly-eval` carries that licence if distributed.
