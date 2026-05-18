# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_dpm_writer.py — Route R2 (B 案 dpm patch + Virtual.mom + Message.txt)."""
from __future__ import annotations

import re
from pathlib import Path

import pytest


# --- Virtual.mom propagator ------------------------------------------------

def test_propagate_virtual_mom_creates_per_segment_dirs(
    tmp_path, virtual_mom_template_path
):
    from abmptools.cg.dpd import propagate_virtual_mom
    segments = ["segA", "segB", "WAT"]
    paths = propagate_virtual_mom(
        template_path=virtual_mom_template_path,
        segments=segments,
        output_dir=tmp_path,
    )
    assert len(paths) == 3
    for seg, p in zip(segments, paths):
        assert p == tmp_path / seg / "Virtual.mom"
        assert p.exists()
        # 全 segment の内容は同一
    contents = [p.read_text() for p in paths]
    assert all(c == contents[0] for c in contents)


def test_propagate_virtual_mom_missing_template(tmp_path):
    from abmptools.cg.dpd import propagate_virtual_mom
    with pytest.raises(FileNotFoundError):
        propagate_virtual_mom(
            template_path=tmp_path / "nonexistent.mom",
            segments=["A"], output_dir=tmp_path,
        )


# --- #Message.txt -----------------------------------------------------------

def test_write_message_txt(tmp_path):
    from abmptools.cg.dpd import write_message_txt
    p = write_message_txt(tmp_path)
    assert p == tmp_path / "#Message.txt"
    assert p.read_text() == "RUN\n"


def test_write_message_txt_custom(tmp_path):
    from abmptools.cg.dpd import write_message_txt
    p = write_message_txt(tmp_path, content="CUSTOM\n")
    assert p.read_text() == "CUSTOM\n"


# --- patch_dpm (B 案 main test) --------------------------------------------

def test_patch_dpm_cholesterol_e2e(
    tmp_path, cholesterol_cg, sample_aij_a_mode,
    dpm_template_path, virtual_mom_template_path,
):
    """cholesterol → cg_segmenter → CGDpdBuilder.build_dpm の full e2e (B 案)."""
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    out_dir = tmp_path / "chol_project"
    dpm = builder.build_dpm(
        template=dpm_template_path,
        output_dir=out_dir,
        virtual_mom_template=virtual_mom_template_path,
        dpm_filename="chol.dpm",
    )
    assert dpm.exists()
    text = dpm.read_text()

    # B 案: def section (J-OCTA spec の class 定義) は template のまま温存
    assert "\\begin{def}" in text
    assert "\\end{def}" in text
    assert "DPDPair:{" in text  # class 定義 (型宣言) 残存
    assert "DPDInput:{" in text  # 5 block (data section の field) は patch 済

    # 5 segment が SegmentModel に登録 (P0..P4)
    for i in range(5):
        assert f'"P{i}"' in text, f"P{i} missing in SegmentModel"

    # DpdBond 4 (cholesterol bond12 = 4 pair)
    bonds = re.findall(r'\{"P\d","P\d","Harmonic"', text)
    assert len(bonds) == 4, f"expected 4 DpdBond, got {len(bonds)}: {bonds}"

    # DpdAngle 3 (cholesterol angle13 = 3 entry)
    angles = re.findall(r'\{"P\d","P\d","P\d"', text)
    assert len(angles) == 3, f"expected 3 DpdAngle, got {len(angles)}"

    # 'chol' polymer 1 entry
    assert text.count('"chol"') >= 1

    # monomer-lib/<seg>/Virtual.mom が全 5 segment 配置されている
    for i in range(5):
        vmom = out_dir / "chol" / "monomer-lib" / f"P{i}" / "Virtual.mom"
        assert vmom.exists(), f"Virtual.mom missing for P{i}"

    # #Message.txt が生成されている
    assert (out_dir / "chol" / "#Message.txt").exists()
    assert (out_dir / "chol" / "#Message.txt").read_text() == "RUN\n"


def test_patch_dpm_missing_template(tmp_path, cholesterol_cg, sample_aij_a_mode):
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    with pytest.raises(FileNotFoundError, match="dpm template not found"):
        builder.build_dpm(
            template=tmp_path / "nonexistent.dpm",
            output_dir=tmp_path / "out",
        )


def test_patch_dpm_invalid_field(
    tmp_path, cholesterol_cg, sample_aij_a_mode, dpm_template_path,
):
    from abmptools.cg.dpd import patch_dpm, CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    with pytest.raises(ValueError, match="unknown field"):
        patch_dpm(
            template_path=dpm_template_path,
            output_path=tmp_path / "out.dpm",
            spec=builder.spec,
            replace_fields=["UnknownField"],
        )
