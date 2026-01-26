from dataclasses import dataclass, asdict

from .. import Quantity, Q_
from ..materials import ConductorMaterial, InsulationMaterial


__all__ = ["SCCalcConfig", "PEConductorConfig"]


@dataclass
class SCCalcConfig:
    """
    Configuration for building sequence short-circuit networks.
    """
    # Which short-circuit case we are preparing impedances for (affects
    # IEC-like factors).
    sc_case: str = "MAX"  # "MAX" or "MIN"

    # Global base apparent power for PU conversion.
    S_base: Quantity = Q_(100.0, "MVA")

    # Voltage factors (tune to your standard/model if needed).
    volt_factor_max: float = 1.10
    volt_factor_min: float = 0.95

    # Conductor temperatures for impedance calculation where required.
    busbar_T_max: Quantity = Q_(20.0, "degC")
    busbar_T_min: Quantity = Q_(150.0, "degC")

    cable_T_max: Quantity = Q_(20.0, "degC")
    cable_T_min: Quantity = Q_(150.0, "degC")

    # If True: ignore induction motors in series impedance of branches
    # (recommended default).
    ignore_induction_motors_in_branch_Z: bool = False

    # If True: add motor *source* branches ground->bus for induction motors
    # above threshold (sequence=1 only).
    include_induction_motor_sources: bool = False
    induction_motor_min_Pn: Quantity = Q_(0.0, "kW")


@dataclass
class PEConductorConfig:
    """
    Configuration settings for sizing PE-conductors.
    """
    cond_mat: ConductorMaterial = ConductorMaterial.COPPER
    insul_mat: InsulationMaterial = InsulationMaterial.PVC
    mech_protected: bool = False  # e.g. PE-conductor in tube or not
    separated: bool = False  # PE-conductor separate from cable or not
    t_interrupt: Quantity = Q_(200, 'ms')  # assumed interruption time of current-protective device = duration of fault current

    def __str__(self) -> str:
        d = asdict(self)
        s_list = []
        for k, v in d.items():
            if isinstance(v, Quantity):
                s_list.append(f"{k}: {v:~P.0f}")
            else:
                s_list.append(f"{k}: {v}")
        return "\n".join(s_list)
