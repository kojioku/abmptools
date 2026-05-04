# -*- coding: utf-8 -*-
"""
abmptools.cg
-------------
Coarse-grained system builders.

abmptools の MO-AAMD-CGMD マルチスケール基盤の CG (粗視化) 系統。
全原子系の ``abmptools.amorphous`` / ``abmptools.membrane`` と対称な位置付けで、
Martini 3 等の粗視化力場向けビルダーを提供する。

サブパッケージ:
    peptide
        Martini 3 ペプチドアモルファス系の end-to-end builder
        (vermouth-martinize 経由、Apache-2.0 互換)

Roadmap (将来予定):
    polymer    -- polyply 経由の高分子鎖
    smallmol   -- Auto-Martini / Bartender 経由の小分子
"""
