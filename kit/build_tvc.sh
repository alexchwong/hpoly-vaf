#!/usr/bin/env bash
# Builds, into WORKDIR/bin:
#   tvc               unmodified TVC 5.18.6
#   tvcutils          unmodified
#   tvc-trace         tvc + hpoly_trace.patch
#   hpoly-eval-pruned tvc-trace with the tested pruning slices applied (kit/pruned/)
#   hpoly-eval        copy of tvc-trace (replace with the pruned binary only after gate X2 of the spec)
# Usage: kit/build_tvc.sh WORKDIR
# Optional env: LAPACK_LINK (default "openblas"), EXTRA_SRC (default empty), JOBS (default 4)
set -euo pipefail
KIT="$(cd "$(dirname "$0")" && pwd)"; W="$(mkdir -p "$1" && cd "$1" && pwd)"; JOBS="${JOBS:-4}"
TS_COMMIT=5fe78dc2f22873845dbbecc7fc2ac8fc33ca42fb
BAMTOOLS_TAG=v2.5.3
ARMA_COMMIT=33036ee0be22ee5f5dd6b0059c97f75be0e1d276
BH_COMMIT=f12f31c9c82b87837e78f0284af88517427a2b9f
fetch(){ # dir url commit
  if [ ! -d "$1/.git" ]; then git clone -q "$2" "$1"; fi; git -C "$1" checkout -q "$3"; }
clean_ts(){ git -C "$W/TS" checkout -q -- Analysis/VariantCaller; }
mkdir -p "$W/deps" "$W/bin"
fetch "$W/TS" https://github.com/iontorrent/TS.git $TS_COMMIT
fetch "$W/deps/bamtools" https://github.com/pezmaster31/bamtools.git $BAMTOOLS_TAG
fetch "$W/deps/RcppArmadillo" https://github.com/RcppCore/RcppArmadillo.git $ARMA_COMMIT
fetch "$W/deps/bh" https://github.com/eddelbuettel/bh.git $BH_COMMIT
clean_ts
if [ ! -f "$W/deps/bt/lib/libbamtools.a" ]; then
  mkdir -p "$W/deps/bamtools/build" && cd "$W/deps/bamtools/build"
  cmake -DCMAKE_INSTALL_PREFIX="$W/deps/bt" -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE=Release .. >/dev/null && make -j"$JOBS" install >/dev/null
  # TVC calls BamWriter::SetNumThreads, which upstream bamtools lacks. Add a no-op.
  sed -i '0,/bool SaveAlignment/s//void SetNumThreads(int) {}\n    bool SaveAlignment/' "$W/deps/bt/include/bamtools/api/BamWriter.h"
fi
build(){ # builddir cmake_source_dir
  mkdir -p "$1" && cd "$1"
  cmake -DTS="$W/TS" -DDEPS="$W/deps" -DLAPACK_LINK="${LAPACK_LINK:-openblas}" -DEXTRA_SRC="${EXTRA_SRC:-}" "$2" >/dev/null
  make -j"$JOBS" >/dev/null; }
echo "[1/3] unmodified tvc and tvcutils"
build "$W/build_unmodified" "$KIT"
cp "$W/build_unmodified/tvc" "$W/bin/tvc"; cp "$W/build_unmodified/tvcutils" "$W/bin/tvcutils"
echo "[2/3] tvc-trace"
git -C "$W/TS" apply "$KIT/hpoly_trace.patch"
build "$W/build_trace" "$KIT"
cp "$W/build_trace/tvc" "$W/bin/tvc-trace"; cp "$W/bin/tvc-trace" "$W/bin/hpoly-eval"
echo "[3/3] hpoly-eval-pruned"
for p in "$KIT"/pruned/0*.patch; do git -C "$W/TS" apply "$p"; done
build "$W/build_pruned" "$KIT/pruned"
cp "$W/build_pruned/tvc" "$W/bin/hpoly-eval-pruned"
clean_ts
"$W/bin/tvc" 2>&1 | head -1; ls "$W/bin"
