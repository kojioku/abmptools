"""Shared fixtures for abmptools test suite."""
import os
import pytest


@pytest.fixture
def sample_dir():
    """Path to the sample data directory."""
    return os.path.join(os.path.dirname(__file__), os.pardir, "sample")


@pytest.fixture
def tmp_workdir(tmp_path, monkeypatch):
    """Change to a temporary directory for tests that write files."""
    monkeypatch.chdir(tmp_path)
    return tmp_path
