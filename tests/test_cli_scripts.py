"""CLIスクリプトの引数パーステスト。"""
import sys
import argparse

import pytest


# ---------------------------------------------------------------------------
# Helper: build an ArgumentParser that mirrors an inline script's parser
# so we can test argument parsing without executing the script body.
# ---------------------------------------------------------------------------


def _build_addsolvfrag_parser():
    """addsolvfrag.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--incoord', help='input frag info',
                        nargs='*', required=True)
    parser.add_argument('-temp', '--template', help='template file',
                        required=True)
    parser.add_argument('-solv', '--solvmol', help='output config file',
                        nargs='*', default=['HOH', 'WAT', 'NA'])
    parser.add_argument('-pb', '--solvation', action='store_true')
    parser.add_argument('-th', '--thicnv', type=float, default=1.0E-5)
    parser.add_argument('-arad', '--atmrad', default='delphi')
    parser.add_argument('-ajfv', '--ajfversion', default='rev23')
    parser.add_argument('-np', '--npro', type=int, default=4)
    parser.add_argument('-nopieda', '--nopieda', action='store_false')
    parser.add_argument('-cmm', '--cmm', action='store_true')
    parser.add_argument('-nocpf', '--nocpf', action='store_false')
    parser.add_argument('-cpfv', '--cpfver', default='10')
    parser.add_argument('-basis', '--basisset', default='6-31G*')
    parser.add_argument('-m', '--method', default='MP2')
    parser.add_argument('-ml', '--mldat', action='store_true')
    parser.add_argument('-mll', '--mllimit', type=int, default=0)
    parser.add_argument('-disp', '--disp', action='store_true')
    parser.add_argument('-dg', '--dgemm', action='store_true')
    parser.add_argument('-rp', '--resp', action='store_true')
    parser.add_argument('-nonbo', '--nonbo', action='store_false')
    parser.add_argument('-mem', '--memory', default='3000')
    parser.add_argument('-lc', '--ligandcharge', nargs=2, action='append')
    parser.add_argument('-rs', '--rsolv', nargs=2, action='append')
    parser.add_argument('-ma', '--manual', action='store_true')
    parser.add_argument('-bsse', '--bsse', action='store_true')
    return parser


def _build_ajf2config_parser():
    """ajf2config.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='*', required=True)
    parser.add_argument('-o', '--output', default='segment_data.dat')
    return parser


def _build_ajfserial_parser():
    """ajfserial.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-t', '--time', type=int, nargs=3)
    parser.add_argument('-str', '--string', default='xxx')
    return parser


def _build_generateajf_parser():
    """generateajf.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--incoord', required=True)
    parser.add_argument('-pb', '--solvation', action='store_true')
    parser.add_argument('-th', '--thicnv', type=float, default=1.0E-5)
    parser.add_argument('-arad', '--atmrad', default='delphi')
    parser.add_argument('-ajfv', '--ajfversion', default='rev23')
    parser.add_argument('-np', '--npro', type=int, default=4)
    parser.add_argument('-nopieda', '--nopieda', action='store_false')
    parser.add_argument('-cmm', '--cmm', action='store_true')
    parser.add_argument('-nocpf', '--nocpf', action='store_false')
    parser.add_argument('-cpfv', '--cpfver', default=23)
    parser.add_argument('-basis', '--basisset', default='6-31G*')
    parser.add_argument('-m', '--method', default='MP2')
    parser.add_argument('-ml', '--mldat', action='store_true')
    parser.add_argument('-mll', '--mllimit', type=int, default=0)
    parser.add_argument('-disp', '--disp', action='store_true')
    parser.add_argument('-dg', '--dgemm', action='store_true')
    parser.add_argument('-rp', '--resp', action='store_true')
    parser.add_argument('-nonbo', '--nonbo', action='store_false')
    parser.add_argument('-mem', '--memory', default='3000')
    parser.add_argument('-lc', '--ligandcharge', nargs=2, action='append')
    parser.add_argument('-rs', '--rsolv', nargs=2, action='append')
    parser.add_argument('-ma', '--manual', default=None)
    parser.add_argument('-ft', '--fragitype', default='ajf',
                        choices=['ajf', 'config'])
    parser.add_argument('-bsse', '--bsse', action='store_true')
    return parser


def _build_getcharge_parser():
    """getcharge.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--incoord', required=True)
    parser.add_argument('-t', '--type', required=True)
    parser.add_argument('-f', '--frag', nargs='*', required=True)
    return parser


def _build_pdb2fmo_parser():
    """pdb2fmo.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--incoord', nargs='*', required=True)
    parser.add_argument('-p', '--parameter', default='input_param')
    parser.add_argument('-xyz', '--xyz', action='store_true')
    parser.add_argument('-noreid', '--norefreshresid', action='store_false')
    parser.add_argument('-noreatm', '--norefreshatmtype', action='store_false')
    return parser


def _build_pdbmodify_parser():
    """pdbmodify.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='*', required=True)
    parser.add_argument('-move', '--move', action='store_true')
    parser.add_argument('-p', '--pos', nargs=3)
    parser.add_argument('-addc', '--addchain', nargs=3)
    parser.add_argument('-mol', '--mol', type=int)
    parser.add_argument('-into', '--into', action='store_true')
    parser.add_argument('-aresname', '--assignresname', action='store_true')
    parser.add_argument('-reid', '--refreshresid', action='store_true')
    parser.add_argument('-reatm', '--refreshatmtype', action='store_true')
    parser.add_argument('-mode', '--mode', default='resnum')
    parser.add_argument('-str', '--string', nargs='*')
    parser.add_argument('-s', '--sort', nargs=2)
    return parser


def _build_udf2fmo_parser():
    """udf2fmo.py のパーサーを再構築する。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--incoord', required=True)
    parser.add_argument('-p', '--parameter', default='input_param')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-s', '--solutes', nargs='*', type=int, default=None)
    parser.add_argument('-r', '--record', type=int, default=None)
    return parser


# ===================================================================
# 1. convertcpf
# ===================================================================
class TestConvertcpfArgs:
    """convertcpf.py の get_args() テスト。"""

    def test_required_input(self, monkeypatch):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog'])
        from abmptools.convertcpf import get_args
        with pytest.raises(SystemExit):
            get_args()

    def test_defaults(self, monkeypatch):
        """デフォルト値が正しく設定される。"""
        monkeypatch.setattr('sys.argv', ['prog', '-i', 'test.cpf'])
        from abmptools.convertcpf import get_args
        args = get_args()
        assert args.input == 'test.cpf'
        assert args.fragments == 0
        assert args.version == 23
        assert args.output is None

    def test_custom_args(self, monkeypatch):
        """カスタム引数が正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'in.cpf', '-f', '5', '-v', '10', '-o', 'out.cpf'
        ])
        from abmptools.convertcpf import get_args
        args = get_args()
        assert args.input == 'in.cpf'
        assert args.fragments == '5'
        assert args.version == '10'
        assert args.output == 'out.cpf'


# ===================================================================
# 2. cpf2ifielist
# ===================================================================
class TestCpf2ifielistArgs:
    """cpf2ifielist.py の getargs() テスト。"""

    @pytest.fixture(autouse=True)
    def _skip_if_no_pandas(self):
        """pandasが利用できない場合スキップする。"""
        pytest.importorskip('pandas')

    def test_required_input(self, monkeypatch):
        """必須引数-iと-fが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog'])
        from abmptools.cpf2ifielist import getargs
        with pytest.raises(SystemExit):
            getargs()

    def test_missing_fragment(self, monkeypatch):
        """必須引数-fが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog', '-i', 'test.cpf'])
        from abmptools.cpf2ifielist import getargs
        with pytest.raises(SystemExit):
            getargs()

    def test_defaults(self, monkeypatch):
        """デフォルト値が正しく設定される。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'test.cpf', '-f', '1', '10'
        ])
        from abmptools.cpf2ifielist import getargs
        args = getargs()
        assert args.input == 'test.cpf'
        assert args.fragment == [1, 10]
        assert args.output == ''
        assert args.andflag is False

    def test_custom_args(self, monkeypatch):
        """カスタム引数が正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'in.cpf', '-f', '5', '20',
            '-o', 'out.csv', '--and'
        ])
        from abmptools.cpf2ifielist import getargs
        args = getargs()
        assert args.input == 'in.cpf'
        assert args.fragment == [5, 20]
        assert args.output == 'out.csv'
        assert args.andflag is True


# ===================================================================
# 3. generate_difie
# ===================================================================
class TestGenerateDifieArgs:
    """generate_difie.py の get_args() テスト。"""

    @pytest.fixture(autouse=True)
    def _skip_if_no_pandas(self):
        """pandasが利用できない場合スキップする。"""
        pytest.importorskip('pandas')

    def test_required_input(self, monkeypatch):
        """必須引数-iと-tが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog'])
        from abmptools.generate_difie import get_args
        with pytest.raises(SystemExit):
            get_args()

    def test_missing_time(self, monkeypatch):
        """必須引数-tが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog', '-i', 'test.cpf'])
        from abmptools.generate_difie import get_args
        with pytest.raises(SystemExit):
            get_args()

    def test_defaults(self, monkeypatch):
        """デフォルト値が正しく設定される。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'test-xxx.cpf', '-t', '1', '10', '1'
        ])
        from abmptools.generate_difie import get_args
        args = get_args()
        assert args.input == 'test-xxx.cpf'
        assert args.time == [1, 10, 1]
        assert args.zero_padding == 1
        assert args.structure == 0
        assert args.fragments == 0
        assert args.version == 23
        assert args.np == 1

    def test_custom_args(self, monkeypatch):
        """カスタム引数が正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'data-xxx.cpf', '-t', '100', '500', '50',
            '-z', '4', '-s', '2', '-f', '10', '-v', '23', '-np', '8'
        ])
        from abmptools.generate_difie import get_args
        args = get_args()
        assert args.input == 'data-xxx.cpf'
        assert args.time == [100, 500, 50]
        assert args.zero_padding == 4
        assert args.structure == 2
        assert args.fragments == '10'
        assert args.version == '23'
        assert args.np == 8


# ===================================================================
# 4. getifiepieda
# ===================================================================
class TestGetifiepiedaArgs:
    """getifiepieda.py の get_args() テスト。"""

    @pytest.fixture(autouse=True)
    def _skip_if_no_pandas(self):
        """pandasが利用できない場合スキップする。"""
        pytest.importorskip('pandas')

    def test_all_defaults(self, monkeypatch):
        """引数なしでもパースできる（全てオプション）。"""
        monkeypatch.setattr('sys.argv', ['prog'])
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.input is None
        assert args.frag is None
        assert args.mol is None
        assert args.molname is None
        assert args.pynp == 1
        assert args.exclude == []
        assert args.nof90so is True
        assert args.noresinfo is True
        assert args.dimene == [2, 1]
        assert args.momene == [1]
        assert args.is_momdim is False
        assert args.zp == 0

    def test_with_frag_and_input(self, monkeypatch):
        """-fと-iを指定してパースできる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-f', '1-10', '101-200', '-i', 'test.log'
        ])
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.frag == ['1-10', '101-200']
        assert args.input == 'test.log'

    def test_custom_args(self, monkeypatch):
        """カスタム引数が正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-f', '10', '-d', '0.8', '-i', 'test.log',
            '-np', '4', '-ex', '102', '103', '-nof90',
            '-dimene', '3', '2', '-momene', '2', '-imd'
        ])
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.frag == ['10']
        assert args.dist == 0.8
        assert args.input == 'test.log'
        assert args.pynp == 4
        assert args.exclude == [102, 103]
        assert args.nof90so is False
        assert args.dimene == [3, 2]
        assert args.momene == [2]
        assert args.is_momdim is True

    def test_multi_mode(self, monkeypatch):
        """--multiモードの引数パースをテストする。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-mul', '1-100', '101-200',
            '-t', '100', '3100', '1000',
            '-i', '["prefix", "-suffix.log"]',
            '-np', '4'
        ])
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.multi == ['1-100', '101-200']
        assert args.time == ['100', '3100', '1000']
        assert args.pynp == 4

    def test_inputx(self, monkeypatch):
        """--inputxオプションが正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-ix', 'file-xxx-suffix.log'
        ])
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.inputx == 'file-xxx-suffix.log'
        assert args.input is None


# ===================================================================
# 5. log2config
# ===================================================================
class TestLog2configArgs:
    """log2config.py の get_args() テスト。"""

    def test_required_input(self, monkeypatch):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog'])
        from abmptools.log2config import get_args
        with pytest.raises(SystemExit):
            get_args()

    def test_defaults(self, monkeypatch):
        """デフォルト値が正しく設定される。"""
        monkeypatch.setattr('sys.argv', ['prog', '-i', 'test.log'])
        from abmptools.log2config import get_args
        args = get_args()
        assert args.input == 'test.log'
        assert args.pynp == 1
        assert args.output == 'segment_data.dat'

    def test_custom_args(self, monkeypatch):
        """カスタム引数が正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'calc.log', '-np', '8', '-o', 'custom.dat'
        ])
        from abmptools.log2config import get_args
        args = get_args()
        assert args.input == 'calc.log'
        assert args.pynp == 8
        assert args.output == 'custom.dat'


# ===================================================================
# 6. log2cpf
# ===================================================================
class TestLog2cpfArgs:
    """log2cpf.py の get_args() テスト。"""

    def test_required_input(self, monkeypatch):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        monkeypatch.setattr('sys.argv', ['prog'])
        from abmptools.log2cpf import get_args
        with pytest.raises(SystemExit):
            get_args()

    def test_defaults(self, monkeypatch):
        """デフォルト値が正しく設定される。"""
        monkeypatch.setattr('sys.argv', ['prog', '-i', 'test.log'])
        from abmptools.log2cpf import get_args
        args = get_args()
        assert args.input == 'test.log'
        assert args.frag is None
        assert args.pynp == 1
        assert args.output is None
        assert args.nof90so is True

    def test_custom_args(self, monkeypatch):
        """カスタム引数が正しくパースされる。"""
        monkeypatch.setattr('sys.argv', [
            'prog', '-i', 'calc.log', '-f', '1', '2', '3',
            '-np', '4', '-o', 'out.cpf', '-nof90'
        ])
        from abmptools.log2cpf import get_args
        args = get_args()
        assert args.input == 'calc.log'
        assert args.frag == ['1', '2', '3']
        assert args.pynp == 4
        assert args.output == 'out.cpf'
        assert args.nof90so is False


# ===================================================================
# 7. addsolvfrag (inline parser — tested via helper)
# ===================================================================
class TestAddsolvfragArgs:
    """addsolvfrag.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iと-tempが指定されない場合SystemExitを発生させる。"""
        parser = _build_addsolvfrag_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_missing_template(self):
        """必須引数-tempが指定されない場合SystemExitを発生させる。"""
        parser = _build_addsolvfrag_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['-i', 'test.pdb'])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_addsolvfrag_parser()
        args = parser.parse_args(['-i', 'test.pdb', '-temp', 'tmpl.ajf'])
        assert args.incoord == ['test.pdb']
        assert args.template == 'tmpl.ajf'
        assert args.solvmol == ['HOH', 'WAT', 'NA']
        assert args.solvation is False
        assert args.thicnv == 1.0E-5
        assert args.atmrad == 'delphi'
        assert args.ajfversion == 'rev23'
        assert args.npro == 4
        assert args.nopieda is True  # store_false default
        assert args.cmm is False
        assert args.nocpf is True  # store_false default
        assert args.cpfver == '10'
        assert args.basisset == '6-31G*'
        assert args.method == 'MP2'
        assert args.mldat is False
        assert args.mllimit == 0
        assert args.dgemm is False
        assert args.resp is False
        assert args.nonbo is True  # store_false default
        assert args.memory == '3000'
        assert args.ligandcharge is None
        assert args.rsolv is None
        assert args.manual is False
        assert args.bsse is False
        assert args.disp is False

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_addsolvfrag_parser()
        args = parser.parse_args([
            '-i', 'a.pdb', 'b.pdb',
            '-temp', 'template.ajf',
            '-solv', 'HOH', 'CL',
            '-pb', '-cmm', '-bsse',
            '-basis', 'cc-pVDZ', '-m', 'HF',
            '-np', '16', '-mem', '8000',
            '-lc', 'LIG', '-1',
        ])
        assert args.incoord == ['a.pdb', 'b.pdb']
        assert args.template == 'template.ajf'
        assert args.solvmol == ['HOH', 'CL']
        assert args.solvation is True
        assert args.cmm is True
        assert args.bsse is True
        assert args.basisset == 'cc-pVDZ'
        assert args.method == 'HF'
        assert args.npro == 16
        assert args.memory == '8000'
        assert args.ligandcharge == [['LIG', '-1']]


# ===================================================================
# 8. ajf2config (inline parser — tested via helper)
# ===================================================================
class TestAjf2configArgs:
    """ajf2config.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        parser = _build_ajf2config_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_ajf2config_parser()
        args = parser.parse_args(['-i', 'frag1.frag', 'frag2.frag'])
        assert args.input == ['frag1.frag', 'frag2.frag']
        assert args.output == 'segment_data.dat'

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_ajf2config_parser()
        args = parser.parse_args(['-i', 'a.frag', '-o', 'custom.dat'])
        assert args.input == ['a.frag']
        assert args.output == 'custom.dat'


# ===================================================================
# 9. ajfserial (inline parser — tested via helper)
# ===================================================================
class TestAjfserialArgs:
    """ajfserial.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        parser = _build_ajfserial_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_ajfserial_parser()
        args = parser.parse_args(['-i', 'template.ajf'])
        assert args.input == 'template.ajf'
        assert args.time is None
        assert args.string == 'xxx'

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_ajfserial_parser()
        args = parser.parse_args([
            '-i', 'tmpl.ajf', '-t', '1', '100', '1', '-str', 'YYY'
        ])
        assert args.input == 'tmpl.ajf'
        assert args.time == [1, 100, 1]
        assert args.string == 'YYY'


# ===================================================================
# 10. generateajf (inline parser — tested via helper)
# ===================================================================
class TestGenerateajfArgs:
    """generateajf.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        parser = _build_generateajf_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_generateajf_parser()
        args = parser.parse_args(['-i', 'mol.pdb'])
        assert args.incoord == 'mol.pdb'
        assert args.solvation is False
        assert args.thicnv == 1.0E-5
        assert args.atmrad == 'delphi'
        assert args.ajfversion == 'rev23'
        assert args.npro == 4
        assert args.nopieda is True
        assert args.cmm is False
        assert args.nocpf is True
        assert args.cpfver == 23
        assert args.basisset == '6-31G*'
        assert args.method == 'MP2'
        assert args.mldat is False
        assert args.mllimit == 0
        assert args.dgemm is False
        assert args.resp is False
        assert args.nonbo is True
        assert args.memory == '3000'
        assert args.ligandcharge is None
        assert args.rsolv is None
        assert args.manual is None
        assert args.fragitype == 'ajf'
        assert args.bsse is False
        assert args.disp is False

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_generateajf_parser()
        args = parser.parse_args([
            '-i', 'protein.pdb', '-pb', '-cmm', '-bsse',
            '-basis', '6-311G**', '-m', 'HF',
            '-np', '8', '-ft', 'config',
            '-ma', 'frags.ajf', '-nopieda', '-nocpf',
        ])
        assert args.incoord == 'protein.pdb'
        assert args.solvation is True
        assert args.cmm is True
        assert args.bsse is True
        assert args.basisset == '6-311G**'
        assert args.method == 'HF'
        assert args.npro == 8
        assert args.fragitype == 'config'
        assert args.manual == 'frags.ajf'
        assert args.nopieda is False
        assert args.nocpf is False


# ===================================================================
# 11. getcharge (inline parser — tested via helper)
# ===================================================================
class TestGetchargeArgs:
    """getcharge.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数が指定されない場合SystemExitを発生させる。"""
        parser = _build_getcharge_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_missing_type(self):
        """必須引数-tが指定されない場合SystemExitを発生させる。"""
        parser = _build_getcharge_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['-i', 'test.log', '-f', '1'])

    def test_missing_frag(self):
        """必須引数-fが指定されない場合SystemExitを発生させる。"""
        parser = _build_getcharge_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['-i', 'test.log', '-t', 'nbo'])

    def test_defaults(self):
        """全必須引数を指定してパースできる。"""
        parser = _build_getcharge_parser()
        args = parser.parse_args(['-i', 'test.log', '-t', 'nbo', '-f', '125'])
        assert args.incoord == 'test.log'
        assert args.type == 'nbo'
        assert args.frag == ['125']

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_getcharge_parser()
        args = parser.parse_args([
            '-i', 'calc.log', '-t', 'mulliken', '-f', '125', '126', '128'
        ])
        assert args.incoord == 'calc.log'
        assert args.type == 'mulliken'
        assert args.frag == ['125', '126', '128']


# ===================================================================
# 12. pdb2fmo (inline parser — tested via helper)
# ===================================================================
class TestPdb2fmoArgs:
    """pdb2fmo.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        parser = _build_pdb2fmo_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_pdb2fmo_parser()
        args = parser.parse_args(['-i', 'mol.pdb'])
        assert args.incoord == ['mol.pdb']
        assert args.parameter == 'input_param'
        assert args.xyz is False
        assert args.norefreshresid is True  # store_false default
        assert args.norefreshatmtype is True  # store_false default

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_pdb2fmo_parser()
        args = parser.parse_args([
            '-i', 'a.pdb', 'b.pdb',
            '-p', 'my_param',
            '-xyz', '-noreid', '-noreatm'
        ])
        assert args.incoord == ['a.pdb', 'b.pdb']
        assert args.parameter == 'my_param'
        assert args.xyz is True
        assert args.norefreshresid is False
        assert args.norefreshatmtype is False


# ===================================================================
# 13. pdbmodify (inline parser — tested via helper)
# ===================================================================
class TestPdbmodifyArgs:
    """pdbmodify.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        parser = _build_pdbmodify_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_pdbmodify_parser()
        args = parser.parse_args(['-i', 'mol.pdb'])
        assert args.input == ['mol.pdb']
        assert args.mode == 'resnum'
        assert args.move is False
        assert args.pos is None
        assert args.addchain is None
        assert args.mol is None
        assert args.into is False
        assert args.assignresname is False
        assert args.refreshresid is False
        assert args.refreshatmtype is False
        assert args.string is None
        assert args.sort is None

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_pdbmodify_parser()
        args = parser.parse_args([
            '-i', 'a.pdb', 'b.pdb',
            '-move', '-p', '10', '20', '30',
            '-addc', '307', '312', 'C',
            '-reid', '-reatm', '-aresname',
            '-mode', 'rename',
            '-str', '001', 'MRT', '002', 'CD7',
            '-s', 'WAT', 'NA',
        ])
        assert args.input == ['a.pdb', 'b.pdb']
        assert args.move is True
        assert args.pos == ['10', '20', '30']
        assert args.addchain == ['307', '312', 'C']
        assert args.refreshresid is True
        assert args.refreshatmtype is True
        assert args.assignresname is True
        assert args.mode == 'rename'
        assert args.string == ['001', 'MRT', '002', 'CD7']
        assert args.sort == ['WAT', 'NA']

    def test_mol_move(self):
        """-molオプションが正しくパースされる。"""
        parser = _build_pdbmodify_parser()
        args = parser.parse_args([
            '-i', 'test.pdb', '-move', '-mol', '3', '-into'
        ])
        assert args.mol == 3
        assert args.move is True
        assert args.into is True


# ===================================================================
# 14. udf2fmo (inline parser — tested via helper)
# ===================================================================
class TestUdf2fmoArgs:
    """udf2fmo.py の引数パーステスト（インラインパーサーの再構築版）。"""

    def test_required_input(self):
        """必須引数-iが指定されない場合SystemExitを発生させる。"""
        parser = _build_udf2fmo_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_defaults(self):
        """デフォルト値が正しく設定される。"""
        parser = _build_udf2fmo_parser()
        args = parser.parse_args(['-i', 'test.udf'])
        assert args.incoord == 'test.udf'
        assert args.parameter == 'input_param'
        assert args.output is None
        assert args.solutes is None
        assert args.record is None

    def test_custom_args(self):
        """カスタム引数が正しくパースされる。"""
        parser = _build_udf2fmo_parser()
        args = parser.parse_args([
            '-i', 'sim.udf', '-p', 'my_param',
            '-o', 'output_name', '-s', '1', '2', '3', '-r', '5'
        ])
        assert args.incoord == 'sim.udf'
        assert args.parameter == 'my_param'
        assert args.output == 'output_name'
        assert args.solutes == [1, 2, 3]
        assert args.record == 5
