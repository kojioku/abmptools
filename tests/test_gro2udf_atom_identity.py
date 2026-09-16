# -*- coding: utf-8 -*-
"""UDF の原子まわりを、受け取る側の流儀に合わせる 2 点。

1. **atom type を畳む** —— openff-interchange は SMIRNOFF に atom type の
   概念が無いため原子 1 個につき 1 型を書く。 PVA 10-mer では 75 型出るが
   中身は 5 種類。 そのまま UDF に持ち込むと受け取る側が扱えない。
   2.16 以降は ``.top`` を書く時点で畳むが、 **それ以前に組んだ ``.top``
   が手元にある**ので、 読み込み時にも畳む。

2. **Atom_ID を全系の通し番号にする** —— 外部ツールが書いた UDF では、
   354 原子の分子で molecule 0 が 0、 1 が 354、 2 が 708 から始まる。
   分子ごとに 0 から振り直すと ID が分子の数だけ重複する。
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "fold_atomtypes"


def _atom_types_in(raw) -> set:
    return {atom[1] for mol in raw.atomlist for atom in mol}


class TestParserFoldsAtomTypes:
    def test_per_atom_types_are_folded_on_read(self):
        """75 型の ``.top`` を読むと、 パラメータの種類数まで畳まれる。"""
        from abmptools.gro2udf.top_parser import TopParser

        raw = TopParser().parse(str(DATA / "system.top"))
        types = _atom_types_in(raw)
        # 元は per-atom unique。 畳めば元素 + 連番の数種類だけが残る。
        assert len(types) < 10, types
        assert all(re.fullmatch(r"[A-Z][a-z]?\d+", t) for t in types), types

    def test_folding_is_idempotent(self, tmp_path):
        """畳み済みの ``.top`` を読んでも、 型名は変わらない (no-op)。"""
        from abmptools.core.top_atomtypes import fold_atomtypes
        from abmptools.gro2udf.top_parser import TopParser

        folded, mapping = fold_atomtypes((DATA / "system.top").read_text())
        assert mapping, "テストデータが既に畳まれている -- 前提が崩れている"
        once = tmp_path / "once.top"
        once.write_text(folded)

        assert _atom_types_in(TopParser().parse(str(once))) == \
               _atom_types_in(TopParser().parse(str(DATA / "system.top")))


class TestAtomIdIsGlobal:
    @pytest.fixture(scope="class")
    def udf_text(self, tmp_path_factory):
        out = tmp_path_factory.mktemp("udf") / "out.udf"
        rc = subprocess.run(
            [sys.executable, "-m", "abmptools.gro2udf", "--from-top",
             str(DATA / "system.top"), str(DATA / "system.gro"),
             "--out", str(out)],
            capture_output=True, text=True,
        )
        assert rc.returncode == 0, rc.stderr[-2000:]
        return out.read_text(errors="replace")

    @staticmethod
    def _atom_ids(text):
        """UDF 中の全 Atom_ID を出現順に返す。

        atom record は ``Atom_ID`` の直後に ``Atom_Name`` (元素記号) が
        並ぶので、 その 2 行の並びで拾う。 他の整数と取り違えない。
        """
        seg = text[text.find("Set_of_Molecules"):]
        return [int(m.group(1)) for m in
                re.finditer(r'^\s*(\d+),\s*\n\s*"([A-Z][a-z]?)",', seg, re.M)]

    def test_ids_do_not_repeat(self, udf_text):
        """分子ごとに 0 から振り直していれば、 ID は分子の数だけ重複する。"""
        ids = self._atom_ids(udf_text)
        assert len(ids) > 100, f"atom record を拾えていない ({len(ids)} 件)"
        assert len(set(ids)) == len(ids), \
            f"Atom_ID が重複している: {len(ids)} 原子に対し {len(set(ids))} 通り"

    def test_ids_are_a_single_running_count(self, udf_text):
        """外部ツールの書き方と同じく、 全系で 0 から連番。"""
        ids = self._atom_ids(udf_text)
        assert ids == list(range(len(ids))), ids[:5] + ["..."] + ids[-3:]

    def test_type_names_are_not_per_atom(self, udf_text):
        """UDF に 75 個の型名が並ばない。"""
        seg = udf_text[udf_text.find("Molecular_Attributes"):]
        names = re.findall(r'"([A-Z][a-z]?\d+)"', seg)
        assert len(set(names)) < 10, sorted(set(names))

class TestBondedTypesCollapseToo:
    """atom type を畳むと、 結合 / 角度 / 二面角の型も一緒に畳まれる。

    parser は結合項の型を **(型名, funct, パラメータ)** で重複判定する
    (``_put_bond_type`` ほか)。 per-atom unique な型名のままだと、 同じ
    パラメータでも名前が違うので 1 本ごとに別の型になり、 UDF に何百もの
    ``Bond_Potential`` / ``Torsion_Potential`` が並ぶ。 名前を畳めばここも
    畳まれる。

    **パラメータは判定キーに入っているので、 中身が違うものは決して
    混ざらない。** ここではその不変条件を、 畳む前後でパラメータの集合が
    変わらないことで押さえる。
    """

    @staticmethod
    def _parsed(fold: bool):
        import abmptools.gro2udf.top_parser as tp
        from abmptools.gro2udf.top_parser import TopParser

        original = tp.fold_atomtypes
        if not fold:
            tp.fold_atomtypes = lambda text: (text, {})
        try:
            return TopParser().parse(str(DATA / "system.top"))
        finally:
            tp.fold_atomtypes = original

    def test_bonded_type_counts_drop(self):
        before, after = self._parsed(False), self._parsed(True)
        for kind in ("bond_types_from_mol", "angle_types_from_mol",
                     "torsion_types_from_mol"):
            n_before = len(getattr(before, kind))
            n_after = len(getattr(after, kind))
            assert n_after < n_before, f"{kind}: {n_before} -> {n_after}"
            assert n_after < 20, f"{kind} still has {n_after} types"

    def test_parameters_are_untouched(self):
        """畳んでもパラメータの集合は 1 つも増減しない。"""
        before, after = self._parsed(False), self._parsed(True)
        for kind in ("bond_types_from_mol", "angle_types_from_mol",
                     "torsion_types_from_mol"):
            sets = [{tuple(t[-1]) for t in getattr(raw, kind)}
                    for raw in (before, after)]
            assert sets[0] == sets[1], (
                f"{kind}: パラメータが変わった\n"
                f"  畳む前のみ: {sets[0] - sets[1]}\n"
                f"  畳んだ後のみ: {sets[1] - sets[0]}"
            )

    def test_differing_parameters_stay_separate(self):
        """同じ型名の組でもパラメータが違えば別の型のまま。"""
        from abmptools.gro2udf.top_parser import get_bond_name

        raw = self._parsed(True)
        seen = {}
        for a1, a2, _funct, params in raw.bond_types_from_mol:
            name = get_bond_name(a1, a2)
            if name in seen:
                # 名前が重なるなら中身は必ず違う (同じなら畳まれているはず)
                assert seen[name] != tuple(params), name
            seen[name] = tuple(params)
