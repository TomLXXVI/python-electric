from enum import StrEnum
from dataclasses import dataclass, field

from python_electric import Quantity, Q_


__all__ = [
    "EarthingSystem",
    "IndirectContactProtResult",
    "get_max_allow_disconnect_time"
]


class EarthingSystem(StrEnum):
    TT = "TT"
    TN = "TN"
    TN_C = "TN-C"
    TN_S = "TN-S"
    TN_CS = "TN-CS"
    IT = "IT"

    def is_TN(self) -> bool:
        return self.value.startswith("TN")

    def is_TT(self) -> bool:
        return self.value.startswith("TT")

    def is_IT(self) -> bool:
        return self.value.startswith("IT")

    def is_TN_C(self) -> bool:
        return self.value == "TN-C"


def get_max_allow_disconnect_time(
    U_phase: Quantity,
    earthing_system: EarthingSystem
) -> Quantity:
    """
    Returns the maximum allowable disconnection time for final circuits with a
    nominal current not exceeding 32 A in low-voltage systems with a TT- or
    TN-earthing scheme according to IEC 60364-4-41.

    Parameters
    ----------
    U_phase: Quantity
        Phase (line-to-neutral/ground) system voltage. Must be higher than 50 V.
    earthing_system: EarthingSystem
        Either a TT- or one of the TN-earthing systems.

    Returns
    -------
    Quantity
    """
    U_phase = U_phase.to('V').m
    if earthing_system.is_TT() or earthing_system.is_TN():
        if 50 < U_phase <= 120:
            return Q_(300, 'ms') if earthing_system.is_TT() else Q_(800, 'ms')
        elif 120 < U_phase <= 230:
            return Q_(200, 'ms') if earthing_system.is_TT() else Q_(400, 'ms')
        elif 230 < U_phase <= 400:
            return Q_(70, 'ms') if earthing_system.is_TT() else Q_(200, 'ms')
        elif U_phase > 400:
            return Q_(40, 'ms') if earthing_system.is_TT() else Q_(100, 'ms')
        else:
            raise ValueError("No result for `U_phase` equal to 50 V or less.")
    else:
        raise ValueError(
            "Function only usable for TT- and TN-earthing schemes."
        )


@dataclass
class IndirectContactProtResult:
    """
    Holds the results returned from methods of `AbstractEarthingSystem`-derived
    classes that check protection against electric shock due to indirect
    contact (insulation faults).

    Attributes
    ----------
    I_f: Quantity, optional
        Fault current due to a line-to-ground fault.
    U_f: Quantity, optional
        Fault voltage (i.e. the touch potential rise) between earth (0 V) and an
        exposed conductive part of the low-voltage network.
    L_max: Quantity, optional
        Maximum allowable cable length between the distribution board and
        the appliance (or sub-distribution board), for which the circuit
        breaker's maximum short-circuit tripping current I_m_max becomes equal
        to the short-circuit current that flows when a line-to-ground fault
        occurs at the end of the cable.
    t_contact_max: Quantity, optional
        Maximum allowable contact duration according to the safety curve.
    R_e_max: Quantity, optional
        Maximum allowable earth spreading resistance of the low-voltage
        distribution network (TN), or of the consumer installation (IT).
    R_pe_max: Quantity, optional
        Maximum allowable resistance of the PE-conductor(s) between the main
        earthing terminal and the exposed conductive part at the end of the
        cable, given the fault voltage U_f and the circuit breaker's maximum
        short-circuit tripping current I_m_max. The actual, measured resistance
        must be smaller than this R_pe_max.
    """
    I_f: Quantity | None = None
    U_f: Quantity | None = None
    L_max: Quantity | None = None
    t_contact_max: Quantity | None = None
    R_e_max: Quantity | None = None
    R_pe_max: Quantity | None = None

    passed: bool = field(init=False, default=False)

    def __str__(self) -> str:
        s = [f"\tpassed: {self.passed}"]
        if self.I_f is not None:
            s.append(f"\tfault current: {self.I_f.to('A'):~P.1f}")
        if self.U_f is not None:
            s.append(f"\tfault voltage: {self.U_f.to('V'):~P.1f}")
        if self.t_contact_max is not None:
            s.append(
                f"\tmaximum allowable fault duration: "
                f"{self.t_contact_max.to('ms'):~P.0f}"
            )
        if self.L_max is not None:
            s.append(
                f"\tmaximum allowable cable length: "
                f"{self.L_max.to('m'):~P.0f}"
            )
        if self.R_e_max is not None:
            s.append(
                f"\tmaximum allowable earth-spreading resistance: "
                f"{self.R_e_max.to('ohm'):~P.0f}"
            )
        if self.R_pe_max is not None:
            s.append(
                f"\tmaximum allowable resistance of the PE-conductor(s): "
                f"{self.R_pe_max.to('mohm'):~P.0f}"
            )
        return "\n".join(s)
