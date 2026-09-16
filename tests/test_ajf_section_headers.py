# -*- coding: utf-8 -*-
"""Namelist headers in a generated ajf must carry no trailing whitespace.

ABINIT-MP's ajf -> inp converters differ in how they find a section:

    V1DD2024      line.find("&ANALYSIS")             substring, tolerant
    Ver.2 Rev.8   re.search("&ANALYSIS$", line, I)   anchored at end of line

``generateajf`` emitted ``"&ANALYSIS "`` with a trailing space. Rev.8's
converter therefore did not recognise the section and **dropped every key in
it**, without a word: the job ran, ``&ANALYSIS`` came out empty in the .inp,
and the log echoed ``ES_RESP = NO`` although the ajf asked for YES. Seen on
18 BRD4 jobs on Fugaku; V1DD2024 had hidden it because its converter is
tolerant.

The source is read by path rather than through ``import abmptools`` on
purpose: an installed copy shadows the working tree, and this invariant is
about the file in this repository.
"""
from __future__ import annotations

import pathlib
import re

SRC = pathlib.Path(__file__).resolve().parent.parent / "abmptools" / "abinit_io.py"
TEXT = SRC.read_text(encoding="utf-8")

#: A header written into a triple-quoted block, e.g. ``&ANALYSIS """``.
TRAILING = re.compile(r'&[A-Za-z0-9_]+[ \t]+"""')
#: A header on its own line inside a block, e.g. ``&POP  \n``.
TRAILING_LINE = re.compile(r'^&[A-Za-z0-9_]+[ \t]+$', re.MULTILINE)


def test_no_section_header_is_written_with_trailing_whitespace():
    bad = TRAILING.findall(TEXT) + TRAILING_LINE.findall(TEXT)
    assert not bad, (
        f"ajf section header(s) emitted with trailing whitespace: {bad}. "
        "Ver.2 Rev.8's mkinp anchors its section match at end of line and "
        "silently drops the whole namelist when it does not match."
    )


def test_analysis_header_specifically():
    assert '&ANALYSIS"""' in TEXT
    assert '&ANALYSIS """' not in TEXT
