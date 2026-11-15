# hpoly-vaf

Homopolymer indel VAF with flow evaluation by code extracted from the Torrent Variant Caller (TVC).

**Start with `docs/SPEC.md`.** Section 3.3 explains this repository, section 3.4 the environment, section 5 the phases and their gates.

## Status

| Part | State |
|---|---|
| `docs/SPEC.md` | Implementation spec v3 with input and output gates for phases P0 to P7 |
| `kit/` | **Working.** Builds unmodified `tvc` 5.18.6, `tvcutils`, `tvc-trace` and `hpoly-eval-pruned`. Ten gates pass on synthetic flow data |
| `kit/check_trimmed_bam.py` | **Working.** Returned `SOFT` on a real example BAM pair |
| `src/hpoly_vaf/` | **Not implemented.** Phases P2 to P6 |
| `tests/` | **Not implemented.** Cases are listed in spec section 7 |

No evaluator binary has been run on real sequencing data.

## Quick start

```bash
git init && git add -A && git commit -m "kit as delivered"

# Docker (preferred)
docker build -t hpoly-vaf-env -f kit/Dockerfile kit
docker run --rm -it -v "$PWD":/repo -v /path/to/data:/data:ro -v /path/to/work:/work -w /repo hpoly-vaf-env

# or conda
conda env create -f kit/environment.yml && conda activate hpoly-vaf

kit/build_tvc.sh /work            # P0: clones pinned sources, builds five binaries into /work/bin
kit/run_synth_gates.sh /work      # P1: expects ten PASS lines
python3 kit/check_trimmed_bam.py --bam merged.bam --trimmed-bam merged.bam.ptrim.bam
```

`build_tvc.sh` clones the full public TVC repository, which is large. Allow 10 GB in the work directory.

## Licence

TVC source is under GPL version 2. The patches in `kit/` and any binary built from them are derived works under that licence.
