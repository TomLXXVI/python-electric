from dataclasses import dataclass
from enum import StrEnum

__all__ = [
    "InsulationProperties",
    "InsulationMaterial",
    "INSULATION_PROPS"
]


@dataclass(frozen=True)
class InsulationProperties:
    type: str
    T_max_cont: float   # degC
    T_max_short: float  # degC


class InsulationMaterial(StrEnum):
    PVC = "PVC"
    RUBBER = "RUBBER"
    XLPE = "XLPE"
    EPR = "EPR"
    PRC = "PRC"
    B = "B"


INSULATION_PROPS: dict[str, InsulationProperties] = {
    "PVC": InsulationProperties("PVC", 70.0, 160.0),
    "RUBBER": InsulationProperties("RUBBER", 60.0, 200.0),
    "XLPE": InsulationProperties("XLPE", 90.0, 250.0),
    "EPR": InsulationProperties("EPR", 90.0, 250.0)
}
