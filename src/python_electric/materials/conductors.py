from dataclasses import dataclass
from enum import StrEnum

__all__ = [
    "ConductorProperties",
    "ConductorMaterial",
    "CONDUCTOR_PROPS"
]


@dataclass(frozen=True)
class ConductorProperties:
    density: float        # kg/m³
    specific_heat: float  # J/(kg.K)
    resistivity20: float  # Ohm.mm²/m
    temp_coeff20: float   # 1/K
    type: str

    def resistivity(self, T: float) -> float:
        """
        Returns the resistivity of the conductor material at the given
        temperature.

        Parameters
        ----------
        T: float
            Temperature of the conductor in degrees Celsius.

        Returns
        -------
        float
        """
        return self.resistivity20 * (1.0 + self.temp_coeff20 * (T - 20.0))


class ConductorMaterial(StrEnum):
    COPPER = "copper"
    ALUMINIUM = "aluminium"
    STEEL = "steel"
    LEAD = "lead"


CONDUCTOR_PROPS: dict[str, ConductorProperties] = {
    "copper": ConductorProperties(8940.0, 385.0, 1.78e-2, 0.0039, "copper"),
    "aluminium": ConductorProperties(2700.0, 897.0, 2.86e-2, 0.0040, "aluminium"),
    "steel": ConductorProperties(7860.0, 450.0, 13.8e-2, 0.0050, "steel"),
    "lead": ConductorProperties(11340.0, 127.0, 20.8e-2, 0.0037, "lead")
}
