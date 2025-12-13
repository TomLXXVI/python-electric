from dataclasses import dataclass
from enum import StrEnum

__all__ = [
    "Insulation",
    "InsulationMaterials",
    "INSULATION_MATERIALS"
]


@dataclass(frozen=True)
class Insulation:
    type: str
    T_max_cont: float   # degC
    T_max_short: float  # degC


class InsulationMaterials(StrEnum):
    PVC = "PVC"
    RUBBER = "RUBBER"
    XLPE = "XLPE"
    EPR = "EPR"
    PRC = "PRC"
    B = "B"


INSULATION_MATERIALS: dict[str, Insulation] = {
    "PVC": Insulation("PVC", 70.0, 160.0),
    "RUBBER": Insulation("RUBBER", 60.0, 200.0),
    "XLPE": Insulation("XLPE", 90.0, 250.0),
    "EPR": Insulation("EPR", 90.0, 250.0)
}
