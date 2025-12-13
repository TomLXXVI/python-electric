from enum import StrEnum
from dataclasses import dataclass
from abc import ABC, abstractmethod

from python_electric import Quantity

Q_ = Quantity

__all__ = [
    "EarthingSystem",
    "IndirectContactProtectionResult"
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
        Line-to-neutral system voltage. Must be higher than 50 V.
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
            "Function only usable for TT- and IT-earthing schemes."
        )


@dataclass
class IndirectContactProtectionResult:
    """
    Holds the results returned from methods of `AbstractEarthingSystem`-derived
    classes that check protection against electric shock due to indirect
    contact (insulation faults).

    Attributes
    ----------
    I_f: Quantity | None
        Fault current due to a line-to-ground fault.
    U_f: Quantity | None
        Voltage (touch potential rise) between earth (0 V) and an exposed
        conductive part of the low-voltage distribution network.
    L_max: Quantity | None
        Maximum allowable length of the cable between the distribution board and
        the appliance/sub-distribution board, so that the minimum short-circuit
        current calculated at the end of the cable would equal the maximum limit
        of the magnetic trip current of the circuit breaker.
    t_c_max: Quantity | None
        Maximum allowable contact duration with touch voltage according to the
        applicable safety curve.
    R_b_max: Quantity | None
        Maximum allowable earth spreading resistance of the low-voltage
        distribution network.
    R_pe_max: Quantity | None
        Maximum allowable resistance of the PE-conductor between the main
        earthing terminal and the most distant exposed conductive part in the
        distribution network, so that the fault current would equal the maximum
        limit of the magnetic trip current of the circuit breaker just upstream
        of the most distant exposed conductive part.
    """
    I_f: Quantity | None = None
    U_f: Quantity | None = None
    L_max: Quantity | None = None
    t_c_max: Quantity | None = None
    R_b_max: Quantity | None = None
    R_pe_max: Quantity | None = None


class AbstractEarthingSystem(ABC):

    @staticmethod
    @abstractmethod
    def check_indirect_contact(*args, **kwargs) -> IndirectContactProtectionResult:
        ...
