"""Offline reference preparation tests, including failure without output."""

import hashlib
import json

import pysam
import pytest

from hpoly_vaf.cli import main
from hpoly_vaf import genexus_reference as gr


@pytest.fixture
def reference_inputs(tmp_path, monkeypatch):
    mt = "ACGTACG"
    monkeypatch.setattr(gr, "RCRS_LEN", len(mt))
    monkeypatch.setattr(gr, "RCRS_MD5", hashlib.md5(mt.encode()).hexdigest())
    bam = tmp_path / "input.bam"
    header = dict(HD=dict(SO="coordinate"), SQ=[
        dict(SN="chrM", LN=7, M5=hashlib.md5(mt.encode()).hexdigest()),
        dict(SN="chr1", LN=4, M5=hashlib.md5(b"ACGT").hexdigest())])
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    hg = tmp_path / "hg.fa"
    hg.write_text(">chr1\nacgt\n>chrM\nAAAAA\n>unwanted\nNNN\n")
    mtfile = tmp_path / "mt.fa"
    mtfile.write_text(">NC_012920.1\n" + mt + "\n")
    out = tmp_path / "prepared.fa"
    args = ["prepare-genexus-reference", "--bam", str(bam), "--trimmed-bam", str(bam),
            "--out", str(out), "--hg19", str(hg), "--mt", str(mtfile)]
    return args, out, mtfile, hg


def test_reference_cli_and_index(reference_inputs, tmp_path):
    args, out, _, _ = reference_inputs
    assert main(args) == 0
    with pysam.FastaFile(str(out)) as fasta:
        assert fasta.references == ["chrM", "chr1"]
        assert fasta.fetch("chrM") == "ACGTACG"
        assert fasta.fetch("chr1") == "ACGT"
    original_index = (tmp_path / "prepared.fa.fai").read_text()
    pysam.faidx(str(out))
    assert (tmp_path / "prepared.fa.fai").read_text() == original_index
    manifest = json.loads((tmp_path / "prepared.fa.manifest.json").read_text())
    assert manifest["original_m5_checked"] == 2
    assert manifest["trimmed_m5_checked"] == 2
    assert manifest["output"]["sha256"] == hashlib.sha256(out.read_bytes()).hexdigest()
    with pytest.raises(SystemExit):
        main(args)


@pytest.mark.parametrize("case", ["bad_mt", "multiple_mt", "missing_contig", "bad_m5", "bad_checksum"])
def test_reference_rejects_bad_inputs(reference_inputs, tmp_path, case):
    args, out, mtfile, hg = reference_inputs
    if case == "bad_mt":
        mtfile.write_text(">wrong\nAAAAAAA\n")
    elif case == "multiple_mt":
        mtfile.write_text(">a\nACG\n>b\nTACG\n")
    elif case == "missing_contig":
        hg.write_text(">chrM\nAAAAA\n")
    elif case == "bad_m5":
        hg.write_text(">chr1\nAAAA\n>chrM\nAAAAA\n")
    else:
        checksum = tmp_path / "md5sum.txt"
        checksum.write_text("0" * 32 + "  hg19.fa.gz\n")
        args += ["--hg19-md5-file", str(checksum)]
    with pytest.raises(SystemExit):
        main(args)
    assert not out.exists()


def test_checksum_entry_required(tmp_path):
    p = tmp_path / "md5.txt"
    p.write_text("abc  another.fa.gz\n")
    with pytest.raises(SystemExit):
        gr.ucsc_md5(p)