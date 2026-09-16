"""``abmptools.trajectory.gen_for_udf`` -- stage 検出と 2 点セット出力。

gmx は呼ばない (monkeypatch)。 ここで守りたいのは
「amorphous の 05_npt_final 決め打ちではない」という 1 点なので、
検証も stage の決まり方と、 決まらない時に黙って進まないこと。
"""

import pytest

from abmptools.trajectory import postprocess as pp


def _stage_files(d, stage, exts=("tpr", "edr", "xtc")):
    for ext in exts:
        (d / f"{stage}.{ext}").write_text("")


# ---------------------------------------------------------------------------
# find_stage
# ---------------------------------------------------------------------------

def test_find_stage_single(tmp_path):
    """stage が 1 つなら名前が何であれそれを使う (amorphous 以外も通る)。"""
    _stage_files(tmp_path, "prod")
    assert pp.find_stage(tmp_path) == "prod"


def test_find_stage_arbitrary_name(tmp_path):
    """Tg 計算のような任意名 1 つでも検出できる。"""
    _stage_files(tmp_path, "tg_cool_300K")
    assert pp.find_stage(tmp_path) == "tg_cool_300K"


def test_find_stage_prefers_amorphous_production(tmp_path):
    """amorphous 5-stage が並ぶときは production stage を選ぶ。"""
    for stage in ("01_em", "02_nvt_highT", "03_npt_highT",
                  "04_anneal", "05_npt_final"):
        _stage_files(tmp_path, stage, exts=("tpr", "edr"))
    assert pp.find_stage(tmp_path) == "05_npt_final"


def test_find_stage_ambiguous_raises_and_lists(tmp_path):
    """決められない時は黙って 1 つ選ばず、 候補を挙げて止まる。"""
    _stage_files(tmp_path, "runA", exts=("tpr", "edr"))
    _stage_files(tmp_path, "runB", exts=("tpr", "edr"))
    with pytest.raises(ValueError) as exc:
        pp.find_stage(tmp_path)
    assert "runA" in str(exc.value) and "runB" in str(exc.value)


def test_find_stage_ignores_lone_tpr(tmp_path):
    """.tpr だけの stage は後処理対象にならない。"""
    (tmp_path / "grompp_only.tpr").write_text("")
    with pytest.raises(FileNotFoundError):
        pp.find_stage(tmp_path)


# ---------------------------------------------------------------------------
# find_ndx
# ---------------------------------------------------------------------------

def test_find_ndx_from_amorphous_build_dir(tmp_path):
    md = tmp_path / "md"
    md.mkdir()
    build = tmp_path / "build"
    build.mkdir()
    (build / "system.ndx").write_text("")
    assert pp.find_ndx(md) == (build / "system.ndx").resolve()


def test_find_ndx_absent_is_none(tmp_path):
    """index の無い run (Tg 等) でも None が返るだけで例外にしない。"""
    assert pp.find_ndx(tmp_path) is None


# ---------------------------------------------------------------------------
# gen_for_udf
# ---------------------------------------------------------------------------

@pytest.fixture
def fake_gmx(monkeypatch):
    """gmx_energy / nojump を呼び出し記録に差し替える。"""
    calls = {}

    def _energy(*, edr, output, terms, gmx="gmx"):
        calls["energy"] = dict(edr=str(edr), output=str(output),
                               terms=list(terms))
        return output

    def _nojump(*, trajectory, tpr, output, group="System", ndx=None,
                gmx="gmx"):
        calls["nojump"] = dict(trajectory=str(trajectory), tpr=str(tpr),
                               output=str(output), group=group,
                               ndx=None if ndx is None else str(ndx))
        return output

    monkeypatch.setattr(pp, "gmx_energy", _energy)
    monkeypatch.setattr(pp, "nojump", _nojump)
    return calls


def test_gen_for_udf_autodetects_stage(tmp_path, fake_gmx):
    _stage_files(tmp_path, "prod")
    res = pp.gen_for_udf(directory=tmp_path)
    assert res["stage"] == "prod"
    assert res["energy"].name == "prod_energy.xvg"
    assert res["trajectory"].name == "prod_nojump.gro"


def test_gen_for_udf_explicit_stage_wins(tmp_path, fake_gmx):
    _stage_files(tmp_path, "05_npt_final")
    _stage_files(tmp_path, "04_anneal")
    res = pp.gen_for_udf(stage="04_anneal", directory=tmp_path)
    assert res["stage"] == "04_anneal"


def test_gen_for_udf_prefers_trr_over_xtc(tmp_path, fake_gmx):
    """.trr は速度も持つので、 両方あれば .trr。"""
    _stage_files(tmp_path, "prod", exts=("tpr", "edr", "xtc", "trr"))
    pp.gen_for_udf(directory=tmp_path)
    assert fake_gmx["nojump"]["trajectory"].endswith("prod.trr")


def test_gen_for_udf_energy_only_when_no_trajectory(tmp_path, fake_gmx):
    """.edr だけの stage は energy を出し、 trajectory は None。"""
    _stage_files(tmp_path, "npt", exts=("tpr", "edr"))
    res = pp.gen_for_udf(directory=tmp_path)
    assert res["energy"] is not None
    assert res["trajectory"] is None
    assert "nojump" not in fake_gmx


def test_gen_for_udf_nothing_to_do_raises(tmp_path, fake_gmx):
    """何も作れない時に RC 0 で「完了」と言わない。"""
    (tmp_path / "prod.tpr").write_text("")
    with pytest.raises(FileNotFoundError):
        pp.gen_for_udf(stage="prod", directory=tmp_path)


def test_gen_for_udf_auto_ndx(tmp_path, fake_gmx):
    md = tmp_path / "md"
    md.mkdir()
    build = tmp_path / "build"
    build.mkdir()
    (build / "system.ndx").write_text("")
    _stage_files(md, "05_npt_final")
    res = pp.gen_for_udf(directory=md)
    assert res["ndx"] == (build / "system.ndx").resolve()
    assert fake_gmx["nojump"]["ndx"].endswith("system.ndx")


def test_gen_for_udf_no_ndx_disables_lookup(tmp_path, fake_gmx):
    """index を持たない run 用に、 自動検出を止められる。"""
    md = tmp_path / "md"
    md.mkdir()
    build = tmp_path / "build"
    build.mkdir()
    (build / "system.ndx").write_text("")
    _stage_files(md, "05_npt_final")
    res = pp.gen_for_udf(directory=md, auto_ndx=False)
    assert res["ndx"] is None
    assert fake_gmx["nojump"]["ndx"] is None


def test_gen_for_udf_terms_range(tmp_path, fake_gmx):
    _stage_files(tmp_path, "prod", exts=("tpr", "edr"))
    pp.gen_for_udf(directory=tmp_path, n_energy_terms=80)
    assert fake_gmx["energy"]["terms"] == list(range(1, 81))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_cli_gen_for_udf_reports_stage(tmp_path, monkeypatch, capsys):
    from abmptools.trajectory import cli

    _stage_files(tmp_path, "prod")

    def _fake(**kwargs):
        return {"stage": "prod", "energy": tmp_path / "prod_energy.xvg",
                "trajectory": None, "ndx": None}

    (tmp_path / "prod_energy.xvg").write_text("")
    monkeypatch.setattr(cli, "gen_for_udf", _fake)
    rc = cli.main(["gen_for_udf", "--dir", str(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "stage: prod" in out
    assert "skipped" in out          # trajectory が無いことを黙らない


def test_cli_gen_for_udf_ambiguous_is_message_not_traceback(tmp_path, capsys):
    from abmptools.trajectory import cli

    _stage_files(tmp_path, "runA", exts=("tpr", "edr"))
    _stage_files(tmp_path, "runB", exts=("tpr", "edr"))
    rc = cli.main(["gen_for_udf", "--dir", str(tmp_path)])
    err = capsys.readouterr().err
    assert rc == 1
    assert "runA" in err and "--stage" in err
