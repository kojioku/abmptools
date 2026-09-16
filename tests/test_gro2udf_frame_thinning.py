"""frame の間引きと、gmx を呼ぶ条件。

**`--max-frames` は「合計何枚か」で、「何本に 1 本か」ではない。** 下流の
UDF の大きさが効くのは枚数なので、そちらを指定できるようにしてある
(``--frame-step`` が従来のストライド)。両方渡すのはエラー。

gro2udf は本来 gmx を呼ばないファイル変換器なので、**gmx が要るのは
``--edr`` と ``--prepare-nojump`` のときだけ**であることも押さえる。
"""
import pytest

from abmptools.gro2udf.top_exporter import _thin
from abmptools.trajectory.postprocess import skip_for_max_frames


class _F:
    """frame の代わり。中身は見ないので識別子だけ持つ。"""
    def __init__(self, i):
        self.i = i
    def __repr__(self):
        return f"F{self.i}"
    def __eq__(self, other):
        return isinstance(other, _F) and other.i == self.i


def _frames(n):
    return [_F(i) for i in range(n)]


@pytest.mark.parametrize("total,want,expect_n", [
    (101, 10, 10),     # 割り切れない
    (101, 25, 21),     # 割り切れない
    (101, 101, 101),   # ちょうど
    (100, 100, 100),
    (10, 100, 10),     # 元がすでに少ない -- 増えない
    (101, 1, 1),
])
def test_max_frames_never_exceeds_the_request(total, want, expect_n):
    got = _thin(_frames(total), max_frames=want)
    assert len(got) == expect_n
    assert len(got) <= want or total <= want


def test_max_frames_keeps_the_right_frames():
    """間引いた先が元の該当 frame であること (順序も保つ)。"""
    got = _thin(_frames(101), max_frames=10)
    step = -(-101 // 10)
    assert got == [_F(i * step) for i in range(len(got))]


def test_frame_step_is_a_stride():
    got = _thin(_frames(10), frame_step=3)
    assert got == [_F(0), _F(3), _F(6), _F(9)]


def test_step_one_returns_everything():
    frames = _frames(5)
    assert _thin(frames, frame_step=1) == frames
    assert _thin(frames) == frames


@pytest.mark.parametrize("bad", [0, -1])
def test_rejects_nonsense(bad):
    with pytest.raises(ValueError):
        _thin(_frames(5), max_frames=bad)
    with pytest.raises(ValueError):
        _thin(_frames(5), frame_step=bad)


@pytest.mark.parametrize("total,want,skip", [
    (101, 100, 2),    # **1 枚超えただけで skip 2 = 約半分になる**
    (101, 10, 11),
    (10000, 100, 100),
    (99, 100, 1),     # 元が少なければ間引かない
])
def test_skip_for_max_frames(total, want, skip):
    assert skip_for_max_frames(total, want) == skip
    assert -(-total // skip) <= max(want, 1) or total <= want


def test_skip_math_matches_the_thinner():
    """gmx 経路 (skip) と読み込み経路 (_thin) が同じ枚数になること。

    別々に実装されているので、片方だけ直して食い違うのを防ぐ。
    """
    for total in (7, 10, 101, 1000):
        for want in (1, 3, 10, 100):
            n_gmx = -(-total // skip_for_max_frames(total, want))
            n_read = len(_thin(_frames(total), max_frames=want))
            assert n_gmx == n_read, (total, want, n_gmx, n_read)
