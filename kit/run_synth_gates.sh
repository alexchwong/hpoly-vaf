#!/usr/bin/env bash
# Usage: ./run_synth_gates.sh WORKDIR    (after build_tvc.sh). Prints PASS/FAIL per gate; exit 1 on any FAIL.
set -uo pipefail
KIT="$(cd "$(dirname "$0")" && pwd)"; W="$(cd "$1" && pwd)"; S="$W/synth"; rm -rf "$S"; mkdir -p "$S"; cd "$S"
FAIL=0; ok(){ if [ "$2" = "$3" ]; then echo "PASS $1"; else echo "FAIL $1: got [$2] want [$3]"; FAIL=1; fi; }
OPTS="--reference ref.fa --process-input-positions-only on --use-input-allele-only on --output-vcf sv.vcf --num-threads 1"
NOSAMP="--downsample-to-coverage 200000"
counts(){ grep -v '^#' "$1/sv.vcf" | cut -f8 | tr ';' '\n' | grep -E '^(FDP|FRO|FAO)=' | sort | tr '\n' ' '; }
body(){ grep -v '^##' "$1/sv.vcf" | md5sum | cut -d' ' -f1; }
python3 "$KIT/make_synth.py" 6000 >/dev/null
"$W/bin/tvc" $OPTS $NOSAMP --input-bam reads.bam --input-vcf hs1.vcf --output-dir a1 >a1.log 2>&1
ok G1_single_allele "$(counts a1)" "FAO=1200 FDP=6000 FRO=4800 "
"$W/bin/tvc" $OPTS $NOSAMP --input-bam reads.bam --input-vcf hs3.vcf --output-dir a3 >a3.log 2>&1
ok G2_three_alleles "$(counts a3)" "FAO=1200,0,0 FDP=6000 FRO=4800 "
"$W/bin/tvc" $OPTS $NOSAMP --input-bam reads.bam --input-vcf hs3.vcf --output-dir a3b >a3b.log 2>&1
ok G3_rerun_identical "$(body a3b)" "$(body a3)"
"$W/bin/tvc" $OPTS --input-bam reads.bam --input-vcf hs3.vcf --output-dir a3d >a3d.log 2>&1
ok G4_default_downsample "$(counts a3d | grep -o 'FDP=[0-9]*')" "FDP=2000"
mkdir -p t3; HPOLY_TRACE_DIR="$S/t3" "$W/bin/tvc-trace" $OPTS $NOSAMP --input-bam reads.bam --input-vcf hs3.vcf --output-dir b3 >b3.log 2>&1
ok G5_trace_vcf_equals_unmodified "$(body b3)" "$(body a3)"
ok G6_groups_tsv "$(cat t3/groups.tsv | tr '\t' ' ')" "chrT 300 CAA CA,C,CAAA 6000 0 4800 1200,0,0"
ok G7_reads_tsv_classes "$(cut -f8 t3/reads.tsv | sort | uniq -c | tr -s ' ' | tr '\n' ';')" " 4800 0; 1200 1;"
mkdir -p t3p; HPOLY_TRACE_DIR="$S/t3p" "$W/bin/hpoly-eval-pruned" $OPTS $NOSAMP --input-bam reads.bam --input-vcf hs3.vcf --output-dir p3 >p3.log 2>&1
ok G8_pruned_vcf_equals_unmodified "$(body p3)" "$(body a3)"
ok G9_pruned_trace_equals_trace "$(cat t3p/reads.tsv t3p/groups.tsv | md5sum | cut -d' ' -f1)" "$(cat t3/reads.tsv t3/groups.tsv | md5sum | cut -d' ' -f1)"
python3 - <<'PY' >/dev/null
import pysam
i=pysam.AlignmentFile("reads.bam"); o=pysam.AlignmentFile("clip.bam","wb",template=i); C=15
for a in i:
    cig=list(a.cigartuples)
    new=[(4,C),(0,cig[0][1]-2*C),(4,C)] if len(cig)==1 else [(4,C),(0,cig[0][1]-C)]+cig[1:-1]+[(0,cig[-1][1]-C),(4,C)]
    a.reference_start+=C; a.cigartuples=new; o.write(a)
o.close(); pysam.index("clip.bam")
PY
mkdir -p t3c; HPOLY_TRACE_DIR="$S/t3c" "$W/bin/tvc-trace" $OPTS $NOSAMP --input-bam clip.bam --input-vcf hs3.vcf --output-dir c3 >c3.log 2>&1
ok G10_soft_clipped_reads_same_trace "$(md5sum < t3c/reads.tsv)" "$(md5sum < t3/reads.tsv)"
exit $FAIL
