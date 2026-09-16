"""``[ atomtypes ]`` を畳んでも力場が変わらないこと。

畳み方を間違えると **エラーは出ずに力場だけ変わる**。形式チェックでは
捕まらないので、ここでは

1. 畳んだ結果の形 (型数・参照の整合) を構文レベルで見る
2. GROMACS がある環境では ``grompp`` + 0 step を流し、**全エネルギー項が
   一致する**ことを見る (無ければ skip)

の 2 段で確かめる。
"""
import shutil
import subprocess
import textwrap

import pytest

from abmptools.core.top_atomtypes import fold_atomtypes

# 2 原子が同じパラメータ、1 原子だけ違う最小の top
MINIMAL = textwrap.dedent("""\
    [ defaults ]
    1 2 yes 0.5 0.8333

    [ atomtypes ]
    ;type  at.num  mass     charge  ptype  sigma   epsilon
    M_0    6       12.011   0.0     A      0.34    0.36
    M_1    6       12.011   0.0     A      0.34    0.36
    M_2    8       15.999   0.0     A      0.30    0.88

    [ moleculetype ]
    MOL  3

    [ atoms ]
    ;  nr  type  resnr  residue  atom  cgnr  charge   mass
       1   M_0   1      MOL      C1    1     -0.1     12.011
       2   M_1   1      MOL      C2    2     -0.2     12.011
       3   M_2   1      MOL      O1    3      0.3     15.999

    [ system ]
    x

    [ molecules ]
    MOL 1
    """)


def _atomtype_names(text):
    names, sec = [], None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            continue
        if sec != "atomtypes":
            continue
        body = s.split(";", 1)[0].strip()
        if body:
            names.append(body.split()[0])
    return names


def _atom_types_used(text):
    used, sec = [], None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            continue
        if sec != "atoms":
            continue
        body = s.split(";", 1)[0].strip()
        if body and not body.startswith(";"):
            f = body.split()
            if len(f) >= 2:
                used.append(f[1])
    return used


def test_identical_parameters_collapse_to_one_type():
    out, mapping = fold_atomtypes(MINIMAL)
    assert _atomtype_names(out) == ["C1", "O1"], "同一パラメータが 1 型に畳まれる"
    assert mapping["M_0"] == mapping["M_1"] == "C1"
    assert mapping["M_2"] == "O1"


def test_atoms_section_is_remapped():
    """畳んだあと [ atoms ] が存在しない型を指していないこと。"""
    out, _ = fold_atomtypes(MINIMAL)
    defined = set(_atomtype_names(out))
    used = set(_atom_types_used(out))
    assert used <= defined, f"未定義の型を参照している: {used - defined}"
    assert used == {"C1", "O1"}


def test_per_atom_charges_survive():
    """電荷は [ atoms ] 側なので、型を畳んでも原子ごとに残る。"""
    out, _ = fold_atomtypes(MINIMAL)
    charges = [l.split()[6] for l in out.splitlines()
               if l.strip().startswith(("1 ", "2 ", "3 ")) and len(l.split()) >= 7]
    assert charges == ["-0.1", "-0.2", "0.3"]


def test_top_with_type_keyed_sections_is_left_alone():
    """bondtypes 等があると型名の引き先が変わるので、触らない。"""
    risky = MINIMAL.replace(
        "[ moleculetype ]",
        "[ bondtypes ]\nM_0  M_1  1  0.15  200000\n\n[ moleculetype ]")
    out, mapping = fold_atomtypes(risky)
    assert mapping == {}, "型名で引くセクションがある top は畳まない"
    assert out == risky


def test_no_duplicates_is_a_no_op():
    single = MINIMAL.replace("M_1    6       12.011   0.0     A      0.34    0.36",
                             "M_1    7       14.007   0.0     A      0.32    0.71")
    out, mapping = fold_atomtypes(single)
    assert mapping == {}
    assert out == single


@pytest.mark.skipif(shutil.which("gmx") is None, reason="gmx が無い")
def test_energy_is_unchanged(tmp_path):
    """**本命**: grompp + 0 step で全エネルギー項が一致すること。

    畳み方の誤りはエネルギーにしか出ない。形式チェックでは通ってしまう。
    """
    import glob
    import os

    here = os.path.dirname(__file__)
    src = os.path.join(here, "data", "fold_atomtypes")
    if not os.path.isdir(src):
        pytest.skip("参照データ (tests/data/fold_atomtypes) が無い")

    gro = os.path.join(src, "system.gro")
    top = os.path.join(src, "system.top")
    orig = tmp_path / "orig.top"
    folded = tmp_path / "folded.top"
    shutil.copy(top, orig)
    text, mapping = fold_atomtypes(open(top, encoding="utf-8").read())
    assert mapping, "このデータは畳めるはず"
    folded.write_text(text, encoding="utf-8")

    mdp = tmp_path / "zero.mdp"
    mdp.write_text(
        "integrator=md\nnsteps=0\ndt=0.002\ncutoff-scheme=Verlet\nnstlist=20\n"
        "pbc=xyz\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\n"
        "DispCorr=EnerPres\nconstraints=h-bonds\n", encoding="utf-8")

    env = dict(os.environ)
    env.pop("OMP_NUM_THREADS", None)   # -ntomp と衝突すると mdrun が落ちる

    energies = {}
    for tag, path in (("orig", orig), ("folded", folded)):
        tpr = tmp_path / f"{tag}.tpr"
        subprocess.run(["gmx", "grompp", "-f", str(mdp), "-c", gro, "-p", str(path),
                        "-o", str(tpr), "-maxwarn", "3"],
                       cwd=tmp_path, check=True, capture_output=True, env=env)
        subprocess.run(["gmx", "mdrun", "-deffnm", tag, "-ntmpi", "1", "-ntomp", "1"],
                       cwd=tmp_path, check=True, capture_output=True, env=env)
        xvg = tmp_path / f"{tag}.xvg"
        subprocess.run(["gmx", "energy", "-f", f"{tag}.edr", "-o", xvg.name],
                       cwd=tmp_path, input=b"1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n0\n",
                       check=True, capture_output=True, env=env)
        rows = [l.split() for l in xvg.read_text().splitlines()
                if l and not l[0] in "#@"]
        energies[tag] = [float(v) for v in rows[0][1:]]

    assert energies["orig"], "エネルギーが読めていない"
    diffs = [abs(a - b) for a, b in zip(energies["orig"], energies["folded"])]
    assert max(diffs) == 0.0, (
        "畳んだら力場が変わった (最大差 %.6g)。**エラーは出ないので、"
        "この検査でしか捕まらない**" % max(diffs))
