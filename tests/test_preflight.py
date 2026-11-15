"""Small synthetic fixtures for P2; no real input data are used here."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pysam
import pytest

from hpoly_vaf.cigar import clip_to_span
from hpoly_vaf.islands import Island, iter_islands, maximum_depth
from hpoly_vaf.join import join_island
from hpoly_vaf.preflight import CheckLog, HardFailure, clip_agrees, ensure_index, run_preflight, sort_observation


def cigar(text):
    import re
    return [("MIDNSHP=X".index(op), int(n)) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", text)]


@pytest.mark.parametrize("pos,text,span,want_pos,want", [
    (100, "20M", (105, 115), 105, "5S10M5S"),
    (100, "5M2I13M", (105, 115), 105, "7S10M3S"),
    (100, "5M3D12M", (105, 115), 108, "5S7M5S"),
    (100, "3S17M", (100, 117), 100, "3S17M"),
    (100, "5M3D12M", (100, 107), 100, "5M12S"),
    (100, "5M2I5M", (100, 105), 100, "5M7S"),
    (100, "5M2I5M", (105, 110), 105, "7S5M"),
    (100, "2H10M2H", (102, 108), 102, "2S6M2S"),
    (100, "5M3N5M", (102, 111), 102, "2S3M3N3M2S"),
    (100, "5M3D5M", (105, 108), 108, "5S5S"),
    (100, "5M2I5M", (100, 110), 100, "5M2I5M"),
])
def test_clipping(pos, text, span, want_pos, want):
    original = cigar(text)
    assert clip_to_span(pos, original, *span) == (want_pos, cigar(want))
    assert original == cigar(text)


def test_unsupported_clipping():
    with pytest.raises(ValueError):
        clip_to_span(0, [(6, 1)], 0, 1)


def header():
    return dict(HD=dict(VN="1.6", SO="coordinate"),
                SQ=[dict(SN="chrT", LN=1000), dict(SN="chrU", LN=1000)],
                RG=[dict(ID="RUN.bc", SM="synthetic", FO="TACG", KS="TCAG")])


def read(h, name="RUN:1", pos=100, text="20M", flag=0, seq=None, tags=True):
    r = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(h))
    r.query_name, r.reference_id, r.reference_start = name, 0, pos
    r.flag, r.cigartuples, r.mapping_quality = flag, cigar(text), 60
    length = sum(n for op, n in r.cigartuples if op in {0, 1, 4, 7, 8})
    r.query_sequence = seq or "A" * length
    if tags:
        from array import array
        r.set_tag("RG", "RUN.bc")
        r.set_tag("ZM", array("h", [256] * 40))
        r.set_tag("ZP", array("f", [0, 0, 0]))
        r.set_tag("ZF", 0)
    return r


def write_bam(path, h, reads):
    with pysam.AlignmentFile(str(path), "wb", header=h) as out:
        for r in reads:
            out.write(r)


@pytest.fixture
def pair(tmp_path):
    h = header()
    bam, trimmed = tmp_path / "original.bam", tmp_path / "trimmed.bam"
    write_bam(bam, h, [read(h), read(h, "RUN:absent", 300)])
    write_bam(trimmed, h, [read(h, pos=105, text="5S10M5S")])
    ref = tmp_path / "ref.fa"
    ref.write_text(">chrT\n" + "A" * 1000 + "\n>chrU\n" + "A" * 1000 + "\n")
    pysam.faidx(str(ref))
    return bam, trimmed, ref, h


def entries(out):
    return json.loads((out / "work/preflight.json").read_text())


def test_full_preflight(pair, tmp_path):
    bam, trimmed, ref, _ = pair
    out = tmp_path / "out"
    before = {p: p.read_bytes() for p in (bam, trimmed, ref)}
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 0
    checks = entries(out)
    assert len(checks) == 11
    assert all(c["status"] == "PASS" for c in checks)
    assert all({"name", "status", "observed", "expected"} <= c.keys() for c in checks)
    assert checks[9]["observed"] == dict(original_without_trimmed=1, minimum_trimmed_mapq=60)
    assert checks[10]["observed"]["maximum_island_depth"] == 1
    assert (out / "work/islands.bed").read_text() == "chrT\t105\t115\n"
    assert all(p.read_bytes() == data for p, data in before.items())
    assert not Path(str(bam) + ".bai").exists()
    assert (out / "work/original.csi").exists()
    assert (out / "work/trimmed.csi").exists()
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2


@pytest.mark.parametrize("case,gate", [
    ("contig", "F1"), ("sort", "F2"), ("rg", "F3"), ("samples", "F3"),
    ("prefix", "F4"), ("ZM", "F5"), ("ZP", "F5"), ("ZF", "F5"),
    ("hardclip", "F6"), ("missing", "F7"), ("strand", "F7"),
    ("seq", "F7"), ("containment", "F7"), ("duplicate_original", "F7"),
    ("duplicate_trimmed", "F7"), ("clip", "F9"),
])
def test_hard_failures(pair, tmp_path, case, gate):
    bam, trimmed, ref, h = pair
    o, t = read(h), read(h, pos=105, text="5S10M5S")
    if case == "contig":
        h["SQ"][0]["LN"] = 999
    if case == "sort":
        h["HD"]["SO"] = "unsorted"
    if case == "rg":
        del h["RG"][0]["KS"]
    if case == "samples":
        h["RG"].append(dict(ID="OTHER.bc", SM="other", FO="TACG", KS="TCAG"))
    if case == "prefix":
        o.query_name = "WRONG:1"
    if case in {"ZM", "ZP", "ZF"}:
        o.set_tag(case, None)
    if case == "hardclip":
        t.cigartuples = cigar("1H5S10M5S")
    if case == "missing":
        t.query_name = "RUN:missing"
    if case == "strand":
        t.flag = 16
    if case == "seq":
        t.query_sequence = "C" * 20
    if case == "containment":
        t.reference_start = 95
    if case == "clip":
        t.cigartuples = cigar("4S10M6S")
    write_bam(bam, h, [o, o] if case == "duplicate_original" else [o])
    write_bam(trimmed, h, [t, t] if case == "duplicate_trimmed" else [t])
    out = tmp_path / case
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2
    assert entries(out)[-1]["name"].startswith(gate)
    assert entries(out)[-1]["status"] == "FAIL"


def test_sort_regressions(tmp_path):
    h = header()
    path = tmp_path / "bad.bam"
    write_bam(path, h, [read(h, pos=200), read(h, pos=100)])
    assert sort_observation(str(path))["regressions"] == 1
    a, b = read(h), read(h)
    a.reference_id = 1
    write_bam(path, h, [a, b])
    assert sort_observation(str(path))["regressions"] == 1


def test_islands_and_flags():
    h = header()
    reads = [read(h, pos=100), read(h, pos=110), read(h, pos=130),
             read(h, pos=140, flag=0x100), read(h, pos=150, flag=0x200),
             read(h, pos=160, flag=0x800), read(h, pos=170, flag=4), read(h, pos=200)]
    other = read(h, pos=20)
    other.reference_id = 1
    assert list(iter_islands(reads + [other])) == [Island("chrT", 100, 150), Island("chrT", 200, 220), Island("chrU", 20, 40)]


def test_depth_and_agreement_boundaries():
    assert maximum_depth([(0, 10), (10, 20)]) == 1
    assert maximum_depth([(0, 10)] * 200000) == 200000
    assert maximum_depth([(0, 10)] * 200001) == 200001
    assert not clip_agrees(998, 1000)
    assert clip_agrees(999, 1000)
    assert clip_agrees(1000, 1000)
    assert not clip_agrees(0, 0)


def test_report_failure_does_not_block(tmp_path):
    log = CheckLog(tmp_path)
    log.add("report", False, "observed", "expected", "REPORT")
    with pytest.raises(HardFailure):
        log.add("hard", False, "observed", "expected")


def test_index_fresh_and_stale(pair, tmp_path):
    bam, _, _, _ = pair
    pysam.index(str(bam))
    index = Path(str(bam) + ".bai")
    assert ensure_index(bam, tmp_path, "a") == str(index)
    os.utime(index, (1, 1))
    assert ensure_index(bam, tmp_path, "a") == str(tmp_path / "a.csi")
    assert index.stat().st_mtime == 1


def test_diagnostic_cap(pair, tmp_path):
    bam, trimmed, _, h = pair
    write_bam(bam, h, [read(h, f"RUN:{n}") for n in range(30)])
    write_bam(trimmed, h, [read(h, f"RUN:{n}", 105, "4S10M6S") for n in range(30)])
    bi = ensure_index(bam, tmp_path, "o")
    ti = ensure_index(trimmed, tmp_path, "t")
    result = join_island((str(bam), str(trimmed), bi, ti, Island("chrT", 105, 115)))
    assert result["joined"] == 30 and result["agreements"] == 0
    assert len(result["mismatches"]) == 20


def test_missing_reference(pair, tmp_path):
    bam, trimmed, _, _ = pair
    out = tmp_path / "missing"
    assert run_preflight(str(bam), str(trimmed), str(tmp_path / "absent.fa"), out) == 2
    assert entries(out)[-1]["status"] == "FAIL"


@pytest.mark.parametrize("field", ["FO", "KS", "SM"])
def test_each_required_rg_field(pair, tmp_path, field):
    bam, trimmed, ref, h = pair
    r = read(h)
    del h["RG"][0][field]
    write_bam(bam, h, [r])
    out = tmp_path / field
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2
    assert entries(out)[-1]["name"] == "F3 read groups"


def test_missing_rg_tag(pair, tmp_path):
    bam, trimmed, ref, h = pair
    r = read(h)
    r.set_tag("RG", None)
    write_bam(bam, h, [r])
    out = tmp_path / "missing-rg"
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2
    assert entries(out)[-1]["name"] == "F4 run id"


def test_empty_inputs_fail_closed(pair, tmp_path):
    bam, trimmed, ref, h = pair
    write_bam(bam, h, [])
    out = tmp_path / "empty"
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2
    assert entries(out)[-1]["name"] == "F4 run id"


def test_empty_after_clip(pair, tmp_path):
    bam, trimmed, _, h = pair
    write_bam(bam, h, [read(h, text="5M3D5M")])
    write_bam(trimmed, h, [read(h, pos=105, text="5S3D5S")])
    bi = ensure_index(bam, tmp_path, "o")
    ti = ensure_index(trimmed, tmp_path, "t")
    result = join_island((str(bam), str(trimmed), bi, ti, Island("chrT", 105, 108)))
    assert result["joined"] == 1 and result["empty_after_clip"] == 1


def test_no_ar_whole_input_not_double_counted(pair, tmp_path):
    bam, trimmed, ref, h = pair
    write_bam(bam, h, [read(h, "RUN:long", 80, "300M"), read(h, "RUN:a", 100), read(h, "RUN:b", 200)])
    write_bam(trimmed, h, [read(h, "RUN:a", 105, "5S10M5S"), read(h, "RUN:b", 205, "5S10M5S")])
    out = tmp_path / "two-islands"
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 0
    assert entries(out)[9]["observed"]["original_without_trimmed"] == 1


def test_f8_duplicate_without_trimmed_counterpart(pair, tmp_path):
    bam, trimmed, ref, h = pair
    write_bam(bam, h, [read(h), read(h, "RUN:extra"), read(h, "RUN:extra")])
    out = tmp_path / "f8"
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2
    assert entries(out)[-1]["name"] == "F8 unique names"


def test_f6_missing_sequence(pair, tmp_path):
    bam, trimmed, ref, h = pair
    r = read(h, pos=105, text="5S10M5S")
    r.query_sequence = None
    write_bam(trimmed, h, [r])
    out = tmp_path / "missing-seq"
    assert run_preflight(str(bam), str(trimmed), str(ref), out) == 2
    assert entries(out)[-1]["name"] == "F6 trimmed SEQ"


@pytest.mark.parametrize("module", ["hpoly_vaf", "hpoly_vaf.cli"])
def test_entry_points(module):
    result = subprocess.run([sys.executable, "-m", module, "--help"], capture_output=True, text=True)
    assert result.returncode == 0 and "preflight" in result.stdout
    result = subprocess.run([sys.executable, "-m", module, "preflight"], capture_output=True)
    assert result.returncode == 2