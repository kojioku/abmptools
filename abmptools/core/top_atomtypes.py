"""GROMACS ``.top`` の ``[ atomtypes ]`` から重複を畳む。

OpenFF Interchange は SMIRNOFF に atom type の概念が無いので、**原子 1 個に
つき 1 型**を書く (``MOL0_0`` … ``MOL0_74``)。パラメータが同じでも別名に
なるため、75 原子の分子なら 75 型が並ぶ。中身は数種類しかない。

これを受け取る側が壊れる。J-OCTA の ``import_gromacs`` は変な型名として
扱い、``gro2udf`` も 1 原子 1 型をそのまま運ぶ。

パラメータ列が**完全に一致する**型だけを 1 つに畳み、名前を元素記号 +
連番 (``C1``, ``O1``, ``H1`` …) に付け替える。電荷は ``[ atoms ]`` 側の
列なので、型を畳んでも原子ごとの電荷は変わらない。

**畳んで安全な条件**: 型名が ``[ atomtypes ]`` の定義と ``[ atoms ]`` の
2 列目にしか現れないこと。``bondtypes`` / ``angletypes`` /
``dihedraltypes`` / ``pairtypes`` / ``nonbond_params`` 等、型名で引く
セクションがある top は**触らない** (畳むと引き先が変わり、エラーを出さずに
力場が変わるため)。OpenFF の出力は結合パラメータを ``[ bonds ]`` 等へ
インラインで書くので、これらのセクションを持たない。

名前について: GAFF 風の名前は騙らない。SMIRNOFF から本物の GAFF 型は
復元できないので、元素記号 + 連番という「由来を主張しない」名前にする。
"""
from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

#: 型名で引くセクション。1 つでもあれば畳まない。
_TYPE_KEYED_SECTIONS = (
    "bondtypes", "angletypes", "dihedraltypes", "pairtypes",
    "constrainttypes", "nonbond_params", "implicit_genborn_params",
)

# 畳んだ型に付ける名前のもとになる元素記号。atomic_number から引く。
_SYMBOLS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 26: "Fe", 29: "Cu",
    30: "Zn", 35: "Br", 53: "I",
}

_SECTION_RE = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")


def _strip_comment(line: str) -> str:
    return line.split(";", 1)[0]


def _section_of(line: str) -> Optional[str]:
    m = _SECTION_RE.match(line)
    return m.group(1).lower() if m else None


def _symbol(atomic_number: str) -> str:
    try:
        return _SYMBOLS.get(int(float(atomic_number)), "X")
    except (TypeError, ValueError):
        return "X"


def fold_atomtypes(text: str) -> Tuple[str, Dict[str, str]]:
    """``[ atomtypes ]`` の重複を畳んだ top を返す。

    Returns
    -------
    (text, mapping)
        ``mapping`` は 旧型名 -> 新型名。畳まなかった場合は入力そのままと
        空の dict を返す。
    """
    lines = text.splitlines()

    # --- 型名で引くセクションがあれば触らない ---
    present = {s for s in (_section_of(l) for l in lines) if s}
    blocked = present.intersection(_TYPE_KEYED_SECTIONS)
    if blocked:
        logger.info(
            "atomtypes は畳まない: 型名で引くセクションがある (%s)。"
            "畳むと引き先が変わる", ", ".join(sorted(blocked)))
        return text, {}

    # --- [ atomtypes ] を読む ---
    # 列は  name [bonded_type] at.num mass charge ptype sigma eps
    # bonded_type は省略されることがあるので、name 以外の全列を同一性の鍵にする。
    section = None
    entries: List[Tuple[int, str, Tuple[str, ...]]] = []   # (行番号, 旧名, 鍵)
    for i, line in enumerate(lines):
        sec = _section_of(line)
        if sec is not None:
            section = sec
            continue
        if section != "atomtypes":
            continue
        body = _strip_comment(line).strip()
        if not body:
            continue
        fields = body.split()
        if len(fields) < 2:
            continue
        entries.append((i, fields[0], tuple(fields[1:])))

    if not entries:
        return text, {}

    # --- 鍵ごとに 1 型へ。名前は元素記号 + 連番 ---
    key_to_new: Dict[Tuple[str, ...], str] = {}
    used: Dict[str, int] = {}
    mapping: Dict[str, str] = {}
    for _, old, key in entries:
        new = key_to_new.get(key)
        if new is None:
            # at.num は bonded_type の有無で位置が動く。数値として読めるほうを取る
            atnum = key[0]
            if not atnum.lstrip("-").isdigit() and len(key) > 1:
                atnum = key[1]
            sym = _symbol(atnum)
            used[sym] = used.get(sym, 0) + 1
            new = "%s%d" % (sym, used[sym])
            key_to_new[key] = new
        mapping[old] = new

    if len(key_to_new) == len(entries):
        logger.debug("atomtypes に重複なし (%d 型)", len(entries))
        return text, {}

    # --- 書き戻し: atomtypes は代表 1 行ずつ、atoms は 2 列目を差し替え ---
    keep_line = {}
    for idx, old, key in entries:
        keep_line.setdefault(key, idx)
    drop = {idx for idx, _, key in entries if keep_line[key] != idx}

    out: List[str] = []
    section = None
    for i, line in enumerate(lines):
        sec = _section_of(line)
        if sec is not None:
            section = sec
            out.append(line)
            continue
        if section == "atomtypes":
            if i in drop:
                continue
            body = _strip_comment(line).strip()
            if body:
                f = body.split()
                f[0] = mapping[f[0]]
                out.append("  ".join(f))
                continue
        elif section == "atoms":
            body = _strip_comment(line)
            if body.strip():
                f = body.split()
                # nr type resnr residue atom cgnr charge [mass]
                if len(f) >= 2 and f[1] in mapping:
                    f[1] = mapping[f[1]]
                    comment = line.split(";", 1)
                    rebuilt = "  ".join(f)
                    out.append(rebuilt + (" ;" + comment[1] if len(comment) > 1 else ""))
                    continue
        out.append(line)

    logger.info("atomtypes を %d -> %d 型に畳んだ", len(entries), len(key_to_new))
    return "\n".join(out) + ("\n" if text.endswith("\n") else ""), mapping


def fold_atomtypes_in_file(path: str) -> Dict[str, str]:
    """``path`` の top をその場で畳む。旧名 -> 新名 の dict を返す。"""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    folded, mapping = fold_atomtypes(text)
    if mapping:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(folded)
    return mapping
