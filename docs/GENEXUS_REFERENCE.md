# Explicit hg19 reference preparation for Genexus BAMs

Some Genexus/Ion Torrent hg19 BAM headers require chrM length 16,569, whereas
the supplied UCSC hg19 archive contains the older 16,571-base chrM. The supplied
reference-preparation script is integrated as
`/home/alex-wong/vscode_projects/hpoly_vaf/src/hpoly_vaf/genexus_reference.py`.
It replaces the mitochondrial sequence with NC_012920.1, rather than trimming
two bases or editing BAM headers. F1 remains unchanged and must subsequently
pass against the prepared FASTA.

Run the explicit preparation stage before preflight, using the fixed P0/P1
container. With the repository at `/repo`, inputs read-only at `/data`, and
external output at `/work`, the command is:

```bash
PYTHONPATH=/repo/src python3 -m hpoly_vaf prepare-genexus-reference \
  --bam /data/merged.bam --trimmed-bam /data/merged.bam.ptrim.bam \
  --hg19 /data/reference/hg19.fa.gz \
  --hg19-md5-file /data/reference/hg19.md5sum.txt \
  --mt /data/reference/NC_012920.1.fa \
  --out /work/reference/hg19_genexus_rcrs.fa
```

The supplied local files can instead be downloaded using the script's original
UCSC/NCBI endpoints by omitting `--hg19` and/or `--mt`. Downloaded files use
`reference_cache` beside the output unless `--cache` is supplied. Downloaded
hg19 archives must match UCSC's checksum manifest. Local archives can be
verified with `--hg19-md5-file`; omission is explicitly recorded as unverified,
not reported as a checksum pass. rCRS sequence length and uppercase sequence MD5
`c68f52674c9fb33aef52dcf399755519` are always checked.

The builder keeps only original BAM contigs, in BAM order; checks every length
and all available original and trimmed BAM M5 tags; writes a FASTA, `.fai`, and
`.manifest.json`; and refuses to overwrite existing outputs. It uppercases the
reference by default. `--keep-case` preserves nuclear soft masking. Input BAMs
and source FASTAs are not modified. Neither the kit nor TVC source is changed.

The manifest records source and output SHA-256 hashes, accession and sequence
MD5, contig lengths/digests, and the number of BAM M5 checks. If BAM M5 tags are
absent, this preparation establishes length compatibility and an explicit
reference recipe, **not sequence identity to an unknown alignment reference**.
Use the original Genexus alignment reference whenever available; a nuclear
length or M5 mismatch is a reason to stop, not an invitation to alter the gate.

Then invoke the existing `preflight` command with the prepared FASTA as
`--reference` and a new empty/absent `--outdir`. A preparation success does not
replace any P2 input or output gate. The later `run` command must use the same
prepared FASTA and shared preflight implementation.

Authoritative background is documented by UCSC in its hg19 description's chrM
note and by Thermo Fisher's Ion Reporter help under chromosome M. The original
supplied explanation has been consolidated here, allowing its temporary root
copy to be deleted after integration validation.