"""Command-line entry points; only implemented phases are exposed."""

import argparse

from .preflight import run_preflight
from .genexus_reference import add_arguments, prepare_reference


def main(argv=None):
    parser = argparse.ArgumentParser(prog="hpoly-vaf")
    commands = parser.add_subparsers(dest="command", required=True)
    preflight = commands.add_parser("preflight", help="validate original and trimmed BAMs")
    for option in ("bam", "trimmed-bam", "reference", "outdir"):
        preflight.add_argument("--" + option, required=True)
    reference = commands.add_parser("prepare-genexus-reference", help="prepare separate hg19/rCRS FASTA for Genexus BAMs")
    add_arguments(reference)
    args = parser.parse_args(argv)
    if args.command == "prepare-genexus-reference":
        return prepare_reference(args)
    return run_preflight(args.bam, args.trimmed_bam, args.reference, args.outdir)


if __name__ == "__main__":
    raise SystemExit(main())