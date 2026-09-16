# -*- coding: utf-8 -*-
"""
Allows running as:  python -m abmptools.gro2udf <input.udf> <input.gro>
"""
import sys
import traceback
from .cli import main

if __name__ == "__main__":
    try:
        # **返り値を捨てないこと。** udf-and-gro モードの Exporter.export()
        # は失敗を例外ではなく 1 で返すので、 捨てると
        # `No module named 'UDFManager'` と表示しながら **RC 0 で終わる**。
        # 呼んだ側 (sample/gro2udf/run.sh や後続のパイプライン) は成功と
        # 読み、 UDF が書かれていないことに後から気付く。
        # --from-top モードは None を返すので or 0 で 0 に落とす。
        sys.exit(main(sys.argv) or 0)
    except SystemExit:
        raise
    except Exception as exc:
        # Print the *full* error message (UDFExportError carries section
        # context + template path + underlying error + hint) plus the
        # Python traceback so users can pinpoint where conversion failed.
        # Earlier versions only printed "ERROR: gro2udf failed." which
        # swallowed all the diagnostic information added by
        # top_exporter.UDFExportError.
        print(f"ERROR: gro2udf failed: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        print("", file=sys.stderr)
        print("--- traceback ---", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
