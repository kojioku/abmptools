# -*- coding: utf-8 -*-
import pytest
import os
import tempfile

from abmptools.udfcreate import udfcreate


class TestInit:
    def test_defaults(self):
        uc = udfcreate()
        assert uc.algo == ['NPT_Andersen_Kremer_Grest']
        assert uc.cellsize == [30, 30, 30]
        assert uc.totalstep == 20000
        assert uc.outstep == 200
        assert uc.tempes == [300]


class TestSetUdfParam:
    @pytest.fixture
    def required_params(self):
        return {
            'algo': ['NPT_Berendsen'],
            'temperature': [400, 500],
            'octahome': '/opt/octa',
            'cognacpath': '/usr/local/bin/cognac',
        }

    def test_all_keys_provided(self, required_params):
        param_udf = {
            **required_params,
            'nvtalgo': ['NVT_Langevin'],
            'packmode': 'cell',
            'density': 0.8,
            'cellsize': [40, 40, 40],
            'totalstep': 10000,
            'outstep': 100,
            'nvtstep': 2000,
            'nvtoutstep': 20,
            'timestep': 0.002,
            'pressure': 2.0,
            'mdlocalnp': 8,
            'mdnp': 16,
        }
        uc = udfcreate()
        uc.setudfparam(param_udf)

        assert uc.algo == ['NPT_Berendsen']
        assert uc.tempes == [400, 500]
        assert uc.octahome == '/opt/octa'
        assert uc.cognacpath == '/usr/local/bin/cognac'
        assert uc.nvtalgo == ['NVT_Langevin']
        assert uc.packmode == 'cell'
        assert uc.density == 0.8
        assert uc.cellsize == [40, 40, 40]
        assert uc.totalstep == 10000
        assert uc.outstep == 100
        assert uc.nvtstep == 2000
        assert uc.nvtoutstep == 20
        assert uc.timestep == 0.002
        assert uc.pressure == 2.0
        assert uc.mdlocalnp == 8
        assert uc.mdnp == 16

    def test_only_required_keys_applies_defaults(self, required_params):
        uc = udfcreate()
        uc.setudfparam(required_params)

        # Required keys are set correctly
        assert uc.algo == ['NPT_Berendsen']
        assert uc.tempes == [400, 500]
        assert uc.octahome == '/opt/octa'
        assert uc.cognacpath == '/usr/local/bin/cognac'

        # Optional keys fall back to defaults
        assert uc.nvtalgo == ['NVT_Nose_Hoover']
        assert uc.packmode == 'density'
        assert uc.density == 0.5
        assert uc.cellsize == [30, 30, 30]
        assert uc.totalstep == 5000
        assert uc.outstep == 50
        assert uc.nvtstep == 1000
        assert uc.nvtoutstep == 10
        assert uc.timestep == 0.001
        assert uc.pressure == 1.0
        assert uc.mdlocalnp == 24
        assert uc.mdnp == 48

    def test_missing_required_key_raises(self):
        """Missing a required key (e.g. 'algo') should raise KeyError."""
        uc = udfcreate()
        with pytest.raises(KeyError):
            uc.setudfparam({
                'temperature': [300],
                'octahome': '/opt/octa',
                'cognacpath': '/usr/local/bin/cognac',
            })


class TestGetConnectData:
    def test_reads_conect_records(self, tmp_path):
        pdb_content = (
            "ATOM      1  C   GLY A   1       1.000   2.000   3.000  1.00  0.00           C\n"
            "ATOM      2  N   GLY A   1       4.000   5.000   6.000  1.00  0.00           N\n"
            "ATOM      3  O   GLY A   1       7.000   8.000   9.000  1.00  0.00           O\n"
            "ATOM      4  H   GLY A   1      10.000  11.000  12.000  1.00  0.00           H\n"
            "CONECT    1    2    3    4\n"
            "CONECT    2    1\n"
            "CONECT    3    1\n"
            "END\n"
        )
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(pdb_content)

        uc = udfcreate()
        data = uc.getconnectdata(str(pdb_file))

        # CONECT lines with >= 4 items are kept; "CONECT    2    1" and "CONECT    3    1"
        # have only 3 items each, so they are excluded.
        assert len(data) == 1
        # Values are converted to 0-indexed integers; "CONECT" remains as-is at index 0
        assert data[0][0] == "CONECT"
        assert data[0][1] == 0   # 1 - 1
        assert data[0][2] == 1   # 2 - 1
        assert data[0][3] == 2   # 3 - 1
        assert data[0][4] == 3   # 4 - 1

    def test_multiple_conect_records(self, tmp_path):
        pdb_content = (
            "CONECT    1    2    3\n"
            "CONECT    2    1    3    4\n"
            "CONECT    3    1    2\n"
            "END\n"
        )
        pdb_file = tmp_path / "test2.pdb"
        pdb_file.write_text(pdb_content)

        uc = udfcreate()
        data = uc.getconnectdata(str(pdb_file))

        # All three CONECT lines have >= 4 items (CONECT + at least 3 numbers)
        # Line 1: "CONECT 1 2 3" -> 4 items -> included
        # Line 2: "CONECT 2 1 3 4" -> 5 items -> included
        # Line 3: "CONECT 3 1 2" -> 4 items -> included
        assert len(data) == 3
        assert data[0] == ["CONECT", 0, 1, 2]
        assert data[1] == ["CONECT", 1, 0, 2, 3]
        assert data[2] == ["CONECT", 2, 0, 1]


class TestGetBatData:
    def test_bond_angle_torsion_extraction(self):
        """Test with a simple linear chain: 0-1-2-3."""
        # connect data as returned by getconnectdata (0-indexed, "CONECT" at [0])
        data = [
            ["CONECT", 0, 1],
            ["CONECT", 1, 0, 2],
            ["CONECT", 2, 1, 3],
        ]
        uc = udfcreate()
        bond, angle, torsion = uc.getbatdata(data)

        # Bonds: from data[i][1]-data[i][j] for j>=2
        # data[0]: bond 0-1
        # data[1]: bond 1-0, bond 1-2
        # data[2]: bond 2-1, bond 2-3
        assert [0, 1] in bond
        assert [1, 0] in bond
        assert [1, 2] in bond
        assert [2, 1] in bond
        assert [2, 3] in bond

        # Angles: for each data[i], pairs (j,k) with j<k from indices >=2
        # data[0]: only one element at j=2 -> no angle
        # data[1]: j=2(val=0), k=3(val=2) -> angle [0, 1, 2]
        # data[2]: j=2(val=1), k=3(val=3) -> angle [1, 2, 3]
        assert [0, 1, 2] in angle
        assert [1, 2, 3] in angle

    def test_no_data(self):
        uc = udfcreate()
        bond, angle, torsion = uc.getbatdata([])
        assert bond == []
        assert angle == []
        assert torsion == []

    def test_single_bond(self):
        """Single CONECT entry with one bond produces one bond, no angles, no torsions."""
        data = [
            ["CONECT", 0, 1],
        ]
        uc = udfcreate()
        bond, angle, torsion = uc.getbatdata(data)
        assert bond == [[0, 1]]
        assert angle == []
        assert torsion == []
