# -*- coding: utf-8 -*-
"""
Tests for abmptools.amorphous.pubchem.

The tests use stubbed ``urllib.request.urlopen`` responses so no real
network access is required.
"""
from __future__ import annotations

import io
import urllib.error
from unittest.mock import patch

import pytest

from abmptools.amorphous import pubchem


class _FakeResp:
    def __init__(self, body: bytes):
        self._body = body

    def read(self) -> bytes:
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _fake_urlopen(body: bytes):
    def _fn(req, timeout=30):
        return _FakeResp(body)
    return _fn


def _fake_urlopen_http_error(code: int):
    def _fn(req, timeout=30):
        raise urllib.error.HTTPError(
            url=getattr(req, "full_url", "<url>"),
            code=code,
            msg=f"HTTP {code}",
            hdrs=None,
            fp=io.BytesIO(b""),
        )
    return _fn


def _fake_urlopen_url_error():
    def _fn(req, timeout=30):
        raise urllib.error.URLError("name resolution failed")
    return _fn


# ---------------------------------------------------------------- happy path


def test_fetch_3d_sdf_returns_body():
    body = b"fake\nSDF\ncontent\nM  END\n"
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen(body)):
        text = pubchem.fetch_3d_sdf(3825)
    assert text == body.decode()


def test_fetch_smiles_strips_whitespace():
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen(b"CCO\n")):
        assert pubchem.fetch_smiles(702) == "CCO"


def test_fetch_3d_sdf_supports_name_lookup():
    url_captured = {}

    def _fn(req, timeout=30):
        url_captured["url"] = req.full_url
        return _FakeResp(b"x")

    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fn):
        pubchem.fetch_3d_sdf("ketoprofen", by="name")
    assert "/compound/name/ketoprofen/SDF" in url_captured["url"]
    assert "record_type=3d" in url_captured["url"]


# ---------------------------------------------------------------- error paths


def test_fetch_3d_sdf_404_becomes_no3d():
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen_http_error(404)):
        with pytest.raises(pubchem.PubChemNo3DError) as exc:
            pubchem.fetch_3d_sdf(99999999)
    # The raised error should carry a helpful message
    assert "No 3D conformer" in str(exc.value)


def test_fetch_smiles_404_becomes_notfound():
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen_http_error(404)):
        with pytest.raises(pubchem.PubChemNotFoundError):
            pubchem.fetch_smiles("this_compound_does_not_exist", by="name")


def test_other_http_error_becomes_generic():
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen_http_error(500)):
        with pytest.raises(pubchem.PubChemError) as exc:
            pubchem.fetch_3d_sdf(3825)
    # HTTP 500 -> generic PubChemError, but NOT a 404-style subclass
    assert not isinstance(exc.value, pubchem.PubChemNotFoundError)


def test_url_error_becomes_pubchem_error():
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen_url_error()):
        with pytest.raises(pubchem.PubChemError):
            pubchem.fetch_3d_sdf(3825)


def test_invalid_by_is_rejected():
    with pytest.raises(ValueError):
        pubchem.fetch_3d_sdf(3825, by="not_a_valid_key")


# ---------------------------------------------------------------- download


def test_download_3d_sdf_writes_file(tmp_path):
    body = b"CID3825\n fake 3D SDF \nM  END\n$$$$\n"
    target = tmp_path / "nested" / "ketoprofen_3d.sdf"
    with patch("abmptools.amorphous.pubchem.urllib.request.urlopen",
               side_effect=_fake_urlopen(body)):
        written = pubchem.download_3d_sdf(3825, str(target))
    assert target.exists()
    assert target.read_bytes() == body
    # Returned path is absolute
    assert written == str(target.resolve())


# ---------------------------------------------------------------- lazy export


def test_lazy_reexport_from_amorphous_namespace():
    # Access through the package-level re-export path
    from abmptools import amorphous as amm
    assert amm.fetch_3d_sdf is pubchem.fetch_3d_sdf
    assert amm.PubChemNo3DError is pubchem.PubChemNo3DError


def test_unknown_attribute_raises():
    from abmptools import amorphous as amm
    with pytest.raises(AttributeError):
        _ = amm.does_not_exist
