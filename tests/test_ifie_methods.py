"""IFIE 表の手法別の列構成 (MP2 / MP3 / MP4(CCPT)) の読み取り。

ABINIT-MP は `Method` によって IFIE 表の列の顔ぶれ自体を変える。MP4 だけ 1 列
多い (GRIMME-MP4) 点と、共有結合対のマスクが手法によらず全数値列に効く点を、
実ログを持ち込まずに小さな合成ログで確かめる。
"""

import logging
import os

import pytest

pytest.importorskip("pandas")

from abmptools.anlfmo import anlfmo, _icolumn_for  # noqa: E402

HARTREE = 627.5095

# 3 行だけの IFIE 表。1 行目は HF-IFIE が -2 Hartree を下回る共有結合対で、
# 読み取り時にゼロへ落とされる。最後の `## Mulliken` は表の終端で、PIEDA 節を
# 持たないログ (MP4 で実在) ではここまで表が続く。
_ROWS = {
    "MP2": (
        "MP2-IFIE",
        "IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   PR-TYPE1"
        "   GRIMME     JUNG       HILL",
        [
            "    2    1    0.000000   F      -15.052480  -0.058601"
            "  -0.056584  -0.049969  -0.045653  -0.041329",
            "    3    1    3.694085   F        0.002822  -0.001209"
            "  -0.001166  -0.000981  -0.000866  -0.000956",
            "    3    2    5.100000   T        0.000100   0.000000"
            "   0.000000   0.000000   0.000000   0.000000",
        ],
    ),
    "MP3": (
        "MP3-IFIE",
        "IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   USER-MP2"
        "   MP3-IFIE   USER-MP3   PADE[2/1]",
        [
            "    2    1    0.000000   F      -15.052480  -0.058601"
            "  -0.058601  -0.050030  -0.054315  -0.049080",
            "    3    1    3.694085   F        0.002822  -0.001209"
            "  -0.001209  -0.000924  -0.001067  -0.000883",
            "    3    2    5.100000   T        0.000100   0.000000"
            "   0.000000   0.000000   0.000000   0.000000",
        ],
    ),
    "CCPT": (
        "MP4-IFIE",
        "IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE  GRIMME-MP2"
        "  MP3-IFIE  GRIMME-MP3  MP4-IFIE  GRIMME-MP4",
        [
            "    2    1    0.000000   F      -15.052480  -0.058601"
            "  -0.058601  -0.050030  -0.054315  -0.050033  -0.054317",
            "    3    1    3.694085   F        0.002822  -0.001209"
            "  -0.001209  -0.000924  -0.001067  -0.001002  -0.001106",
            "    3    2    5.100000   T        0.000100   0.000000"
            "   0.000000   0.000000   0.000000   0.000000   0.000000",
        ],
    ),
}


def _write_log(tmp_path, method):
    heading, header, rows = _ROWS[method]
    text = [
        "     ## READ NAMELIST",
        "         Method                  = {}".format(method),
        "",
        "     ## {}".format(heading),
        "",
        "           " + header,
        "                      / A      APPROX.   / Hartree",
        "        " + "-" * 90,
    ]
    text += ["        " + r for r in rows]
    text += ["", "   ## Mulliken atomic population", "", "     1   N   7.85  -0.85"]
    path = os.path.join(str(tmp_path), "{}.log".format(method))
    with open(path, "w") as f:
        f.write("\n".join(text) + "\n")
    return path


def _read(path):
    obj = anlfmo()
    obj.f90soflag = False
    obj.anlmode = "frag"
    obj.tgt1frag = [1]
    obj.logMethod = obj.getlogmethod(path)
    obj.icolumn = _icolumn_for(obj.logMethod)
    ifie = obj.read_ifiepieda(path)[0]
    return obj, obj.getifiedf(ifie)


@pytest.mark.parametrize("method", ["MP2", "MP3", "CCPT"])
def test_method_is_detected(tmp_path, method):
    obj, _ = _read(_write_log(tmp_path, method))
    assert obj.logMethod == method


@pytest.mark.parametrize(
    "method,columns",
    [
        ("MP2", ["HF-IFIE", "MP2-IFIE", "PR-TYPE1", "GRIMME", "JUNG", "HILL"]),
        ("MP3", ["HF-IFIE", "MP2-IFIE", "USER-MP2", "MP3-IFIE", "USER-MP3",
                 "PADE[2/1]"]),
        ("CCPT", ["HF-IFIE", "MP2-IFIE", "GRIMME-MP2", "MP3-IFIE",
                  "GRIMME-MP3", "MP4-IFIE", "GRIMME-MP4"]),
    ],
)
def test_columns_match_the_method(tmp_path, method, columns):
    _, df = _read(_write_log(tmp_path, method))
    assert list(df.columns) == ["I", "J", "DIST", "DIMER-ES"] + columns


@pytest.mark.parametrize("method", ["MP2", "MP3", "CCPT"])
def test_all_rows_are_read(tmp_path, method):
    _, df = _read(_write_log(tmp_path, method))
    # 終端の `## Mulliken` 見出しが表の行として混入しないこと。
    assert len(df) == 3


@pytest.mark.parametrize("method", ["MP2", "MP3", "CCPT"])
def test_values_are_converted_to_kcal_per_mol(tmp_path, method):
    _, df = _read(_write_log(tmp_path, method))
    row = df[(df["I"] == 3) & (df["J"] == 1)].iloc[0]
    raw = [float(x) for x in _ROWS[method][2][1].split()[4:]]
    for col, value in zip(list(df.columns)[4:], raw):
        assert row[col] == pytest.approx(value * HARTREE, abs=1e-6)


@pytest.mark.parametrize("method", ["MP2", "MP3", "CCPT"])
def test_bonded_pair_is_zeroed_across_every_value_column(tmp_path, method):
    """HF-IFIE < -2 Hartree の対は分散補正や高次項も含めて全列落とす。"""
    _, df = _read(_write_log(tmp_path, method))
    row = df[(df["I"] == 2) & (df["J"] == 1)].iloc[0]
    for col in list(df.columns)[4:]:
        assert row[col] == 0.0, col


def test_mp4_keeps_the_column_the_fortran_reader_drops(tmp_path):
    """GRIMME-MP4 は MP4 表の 7 列目。Fortran リーダは 6 値しか読まない。"""
    _, df = _read(_write_log(tmp_path, "CCPT"))
    assert "GRIMME-MP4" in df.columns
    row = df[(df["I"] == 3) & (df["J"] == 1)].iloc[0]
    assert row["GRIMME-MP4"] == pytest.approx(-0.001106 * HARTREE, abs=1e-6)


def test_unreadable_ifie_table_warns(tmp_path, caplog):
    """見出しを取り違えて 0 行になったら黙って通さない。"""
    path = os.path.join(str(tmp_path), "odd.log")
    with open(path, "w") as f:
        f.write(
            "     ## READ NAMELIST\n"
            "         Method                  = MP3\n"
            "\n"
            "     ## MP9-IFIE\n"
            "\n"
            "           IJ-PAIR    DIST\n"
            "                      / A\n"
            "        ----------\n"
            "            2    1    0.000000   F      -1.0  -1.0  -1.0  -1.0  -1.0  -1.0\n"
            "\n"
            "     ## PIEDA\n"
            "\n"
            "           IJ-PAIR    ES\n"
            "                      / kcal\n"
            "        ----------\n"
            "            2    1    0.1   0.2   0.3   0.4   0.5\n"
        )
    obj = anlfmo()
    obj.f90soflag = False
    obj.anlmode = "frag"
    obj.tgt1frag = [1]
    obj.logMethod = obj.getlogmethod(path)
    with caplog.at_level(logging.WARNING):
        ifie, pieda = obj.read_ifiepieda(path)[:2]

    # PIEDA が読めているので「何も読めなかった」の警告条件はすり抜ける。
    assert not ifie
    assert pieda
    assert any("no IFIE rows" in r.getMessage() for r in caplog.records)
