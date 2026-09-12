# -*- coding: utf-8 -*-
"""
abmptools.gro2udf
-----------------
GROMACS の .gro / .top を COGNAC の UDF に変換する。

**入口が 2 つある。何を持っているかで選ぶ。**

=========================  ==========================  ======================
持っているもの             使うクラス                  出るもの
=========================  ==========================  ======================
UDF (力場つき) + .gro      :class:`Exporter`           座標だけ差し替えた UDF
.top + .gro                :class:`TopExporter`        力場ごと組んだ UDF
=========================  ==========================  ======================

:class:`Exporter` は**テンプレートの UDF が力場を持っている**ことが前提で、
座標と セルだけを書き換える。COGNAC で組んだ系を GROMACS に出して戻す往復に
使う。:class:`TopExporter` は ``.top`` からトポロジと力場を組み立てるので、
**元になる UDF が要らない**。GROMACS だけで作った系を UDF にするのはこちら。

座標を差し替えるだけ (UDF がテンプレート)::

    from abmptools.gro2udf import Exporter
    Exporter().export("test.udf", "output.gro")
    # writes: test_groout.udf  in current directory

.top から組み立てる::

    from abmptools.gro2udf import TopExporter
    TopExporter().export("system.top", "output.gro",
                         template_path="template.udf",
                         out_path="result.udf")

CLI usage (via __main__.py or standalone script)::

    # 座標を差し替えるだけ
    python -m abmptools.gro2udf test.udf output.gro

    # .top から組み立てる
    python -m abmptools.gro2udf --from-top system.top output.gro \\
        [--template template.udf] [--out result.udf]

"""
from .exporter import Exporter
from .top_exporter import TopExporter

__all__ = ["Exporter", "TopExporter"]
