import math
from enum import StrEnum, Enum

__all__ = [
    "VoltReference",
    "PhaseSystem"
]


class VoltReference(float, Enum):
    """
    Enum-class to indicate what kind of voltage is meant. The member values
    are the coefficients to be used in the equation for calculating voltage
    drop.
    """
    PH3_LINE_TO_LINE = math.sqrt(3.0)
    PH3_GROUND_TO_LINE = 1.0
    PH1 = 2.0
    DC = 2.0


class PhaseSystem(StrEnum):
    THREE_PHASE = "three_phase"
    SINGLE_PHASE = "single_phase"

    @property
    def three_phase(self) -> bool:
        return self == PhaseSystem.THREE_PHASE

    @property
    def single_phase(self) -> bool:
        return self == PhaseSystem.SINGLE_PHASE

    def cP(self, volt_ref: VoltReference | None = None) -> float:
        """
        Returns the coefficient in the power equation `P = cP * U * I * cos_phi`
        depending on the phase system and the voltage reference.

        Parameters
        ----------
        volt_ref:
            Indicates if the power equation uses either phase or line voltage
            and current. This parameter is only relevant with three-phase
            systems and can be left as None in the case of single-phase systems.

        Returns
        -------
        float
        """
        if self == PhaseSystem.THREE_PHASE:
            if volt_ref == VoltReference.PH3_LINE_TO_LINE or volt_ref is None:
                return math.sqrt(3)
            elif volt_ref == VoltReference.PH3_GROUND_TO_LINE:
                return 3.0
            else:
                raise ValueError(
                    f"Voltage reference {volt_ref} is not "
                    f"valid in case of a three-phase cable."
                )
        elif self == PhaseSystem.SINGLE_PHASE:
            return 1.0
        else:
            raise ValueError("Value of `phase system` is not recognized.")
