from dataclasses import dataclass, asdict

from .. import Quantity, Q_
from ..materials import ConductorMaterial, InsulationMaterial


__all__ = ["SCCalcConfig", "PEConductorConfig"]


@dataclass
class SCCalcConfig:
    """
    Configuration for building sequence short-circuit networks.

    Attributes
    ----------
    S_base: Quantity, default Q_(100.0, "MVA")
        Global base apparent power for PU conversion.
    volt_factor_max: float, default 1.10
        Voltage factor used for calculating maximum short-circuit currents.
    volt_factor_min: float, default 0.95
        Voltage factor used for calculating minimum short-circuit currents.
    busbar_T_max: Quantity, default Q_(20.0, "degC")
        Temperature of bus bars used for calculating maximum short-circuit
        currents.
    busbar_T_min: Quantity, default Q_(150.0, "degC")
        Temperature of bus bars used for calculating minimum short-circuit
        currents.
    cable_T_max: Quantity, default Q_(20.0, "degC")
        Temperature of cables used for calculating maximum short-circuit
        currents.
    cable_T_min: Quantity, default Q_(150.0, "degC")
        Temperature of cables used for calculating minimum short-circuit
        currents.
    ignore_induction_motors_in_branch_Z: bool, default False
        Indicates whether the impedance of induction motors must be included
        in the short-circuit calculations.
    include_induction_motor_sources: bool, default False
        Indicates whether induction motors must be included as supply sources
        in the short-circuit calculations.
    induction_motor_min_Pn: Quantity, default Q_(0.0, "kW")
        Specifies the minimum nominal power an induction motor must have to be
        included as a source in the short-circuit calculations.
    """
    # Used internally by the Sequence Network Builder to see whether minimum or
    # maximum short-circuit currents are to be calculated.
    sc_case: str = "MAX"

    S_base: Quantity = Q_(100.0, "MVA")
    volt_factor_max: float = 1.10
    volt_factor_min: float = 0.95
    busbar_T_max: Quantity = Q_(20.0, "degC")
    busbar_T_min: Quantity = Q_(150.0, "degC")
    cable_T_max: Quantity = Q_(20.0, "degC")
    cable_T_min: Quantity = Q_(150.0, "degC")
    ignore_induction_motors_in_branch_Z: bool = False
    include_induction_motor_sources: bool = False
    induction_motor_min_Pn: Quantity = Q_(0.0, "kW")

@dataclass
class PEConductorConfig:
    """
    Configuration settings for sizing PE-conductors.

    Attributes
    ----------
    cond_mat: ConductorMaterial, default ConductorMaterial.COPPER
        Material the core of the PE-conductor(s) is made of.
    insul_mat: InsulationMaterial, default InsulationMaterial.PVC
        Insulation material around the PE-conductors.
    mech_protected: bool, default False
        Indicates if the PE-conductors are mechanically protected, e.g.
        PE-conductors inside tubes.
    separated: bool, default False
        Indicates whether PE-conductors are separate from the cables or not.
    t_interrupt: Quantity, default Q_(200, 'ms')
        Assumed interruption time of current-protective devices = duration of
        the fault current to earth.
    """
    cond_mat: ConductorMaterial = ConductorMaterial.COPPER
    insul_mat: InsulationMaterial = InsulationMaterial.PVC
    mech_protected: bool = False
    separated: bool = True
    t_interrupt: Quantity = Q_(200, 'ms')

    def __str__(self) -> str:
        d = asdict(self)
        s_list = []
        for k, v in d.items():
            if isinstance(v, Quantity):
                s_list.append(f"{k}: {v:~P.0f}")
            else:
                s_list.append(f"{k}: {v}")
        return "\n".join(s_list)
