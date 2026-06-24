"""abmptools.udfcharge — OCTA/COGNAC UDF への per-atom 電荷割り当て。

単分子 UDF (電荷あり) → バルク UDF (同名分子、 電荷なし) へ電荷を転写する。

>>> from abmptools.udfcharge import read_molecule_charges, assign_charges_to_bulk
>>> tmpl = read_molecule_charges("mol.udf")
>>> assign_charges_to_bulk("bulk.udf", tmpl, "bulk_charged.udf")

CLI: ``python -m abmptools.udfcharge --template mol.udf --bulk bulk.udf --out bulk_charged.udf``
"""

from .core import (
    CHARGE_UNIT,
    POINT_CHARGE,
    AssignResult,
    MoleculeChargeTemplate,
    RestoreResult,
    assign_charges_to_bulk,
    read_molecule_charges,
    restore_formal_charge,
)

__all__ = [
    "CHARGE_UNIT",
    "POINT_CHARGE",
    "AssignResult",
    "MoleculeChargeTemplate",
    "RestoreResult",
    "assign_charges_to_bulk",
    "read_molecule_charges",
    "restore_formal_charge",
]
