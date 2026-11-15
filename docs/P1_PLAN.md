# P1 execution plan: confirm the engine on synthetic data

Status: planned only; not executed locally.
Prerequisite: P0 passed on 15 November 2026; see the appended local entry in
`/home/alex-wong/vscode_projects/hpoly_vaf/REPORT.md`.

## Scope and constraints

P1 already has an implementation in the delivered kit. No new Python modules,
test framework, dependencies, or kit edits are needed. Execute the unchanged
synthetic generator and gate script in the exact Docker image used for P0.
Do not start P2 or substitute the pruned evaluator for `hpoly-eval`.

- Repository: `/home/alex-wong/vscode_projects/hpoly_vaf`.
- Work directory: `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder`.
- Image: `sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43`.
- Container paths: repository `/repo` (read-only), work directory `/work` (read-write).
- Logs and synthetic sequencing outputs remain outside the public repository.

## 1. Reconfirm the P1 input gate

Inside the same image, repeat all three P0 output checks:

1. `/work/bin/tvc 2>&1 | head -1` must exactly equal
   `tvc 5.18-6 (5fe78dc2) - Torrent Variant Caller`.
2. `/work/bin` must contain all five executable binaries: `hpoly-eval`,
   `hpoly-eval-pruned`, `tvc`, `tvc-trace`, `tvcutils`.
3. `cmp /work/bin/hpoly-eval /work/bin/tvc-trace` must exit 0.

Check the binary SHA-256 values against the P0 report and confirm that the image
ID still matches. Do not rebuild or change compiler/environment between phases.
On a HARD gate failure, record command, expected value and observed value in
the report and stop; do not run the synthetic script.

## 2. Execute the delivered synthetic gates

**Destructive scope:** the kit script deletes and recreates `/work/synth`, which
maps to `/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/synth`.
Before reruns, preserve any earlier synthetic outputs elsewhere in the external
work directory. The first run has no existing synthetic output to preserve.

Run from a Bash shell on the host, capturing the real Docker/script exit status:

```bash
set -o pipefail
docker run --rm --name hpoly-vaf-p1 \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v /home/alex-wong/vscode_projects/hpoly_vaf:/repo:ro \
  -v /home/alex-wong/vscode_projects/hpoly_vaf_workfolder:/work \
  -w /repo \
  sha256:e55e194d82e5f7442f0396def429a51c7662209c6c5f8c88ce1c9e1aa8c32b43 \
  /repo/kit/run_synth_gates.sh /work \
  2>&1 | tee /home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-gates.log
rc=${PIPESTATUS[0]}
printf '%s\n' "$rc" > /home/alex-wong/vscode_projects/hpoly_vaf_workfolder/logs/p1-gates.exit
```

The generator creates deterministic flow-bearing synthetic data with 6,000
reads: 4,800 REF and 1,200 deletion ALT reads, forward and reverse strands,
fixed random seeds, and `ZP=0,0,0`. It also creates the reference and one-ALT
and three-ALT hotspot inputs. The script then runs unmodified, traced, pruned,
and soft-clipped evaluations.

## 3. Require every P1 output gate

Success requires script exit status 0 and all ten expected `PASS` lines.

| Gate | Required result |
|---|---|
| G1 | One ALT: `FAO=1200 FDP=6000 FRO=4800` |
| G2 | Three ALTs: `FAO=1200,0,0 FDP=6000 FRO=4800` |
| G3 | Repeated unmodified run has identical VCF body |
| G4 | Default sampling gives `FDP=2000` |
| G5 | Traced VCF body equals unmodified VCF body |
| G6 | Group trace: `chrT 300 CAA CA,C,CAAA 6000 0 4800 1200,0,0` |
| G7 | Read trace classes: 4,800 class 0 and 1,200 class 1 |
| G8 | Pruned VCF body equals unmodified VCF body |
| G9 | Pruned read and group traces equal traced-build traces |
| G10 | Soft-clipped read trace equals unclipped read trace |

Check the gate identifiers, not just the total number of PASS lines. Retain
the full console log and script exit status. The delivered script performs
all ten checks in one invocation; do not rewrite its control flow or gates.

## 4. Preserve evidence and record the outcome

Synthetic inputs, per-run engine logs, VCFs, and traces will be under
`/home/alex-wong/vscode_projects/hpoly_vaf_workfolder/synth/`:

- `a1`, `a3`, `a3b`, `a3d`: unmodified engine runs.
- `b3`, `t3`: traced engine VCF and trace.
- `p3`, `t3p`: pruned engine VCF and trace.
- `c3`, `t3c`: soft-clipped engine VCF and trace.
- `ref.fa`, `reads.bam`, `clip.bam`, their indexes, and hotspot VCFs.
- Corresponding engine `.log` files.

Append a new local P1 entry to
`/home/alex-wong/vscode_projects/hpoly_vaf/REPORT.md` containing the exact
command, image ID, script exit status, all ten gate outcomes, and log location.
Do not overwrite earlier report entries or confuse historical kit-author
results with this local run. Commit and push only the text evidence; do not
publish BAMs, synthetic run directories, or compiled binaries.

If the script fails, record the failed command, expected and observed results,
and relevant engine-log excerpts, then stop. Do not edit kit files, relax
gates, or automatically retry with altered parameters.

## 5. Completion boundary

Declare P1 complete only after every required check passes and the result is
recorded. P1 does not implement the Python CLI, perform real-data preflight,
or establish real-data pruning equality. `hpoly-eval` must remain the copy of
`tvc-trace`; adopting the pruned binary requires P7 gates on real data.

The next phase, P2, needs paths to the original BAM, trimmed BAM and reference
FASTA, plus a `VERDICT: SOFT` result from the delivered trimmed-BAM checker.
No P2 work is part of this P1 plan.