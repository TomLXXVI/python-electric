from typing import Union
import math
import numpy.typing as npt

from ... import Quantity

__all__ = [
    "PerUnitSystem",
    "convert_per_unit_impedance"
]


class PerUnitSystem:

    def __init__(self, S_base: Quantity, U_base: Quantity):
        """
        Creates a `PerUnitSystem` object.
        
        Parameters
        ----------
        S_base:
            Base power of the per-unit system.
        U_base:
            Base voltage of the per-unit system (line-to-line voltage).
        """
        self.S_base = S_base
        self.U_base = U_base
        self.I_base: Quantity = S_base / (math.sqrt(3) * U_base)
        self.Z_base: Quantity = U_base ** 2 / S_base

    def get_per_unit_impedance(self, Z_act: Quantity) -> Union[float, complex]:
        """Returns the per-unit value of the actual impedance `Z_act`."""
        Z_pu: Quantity = Z_act / self.Z_base
        return Z_pu.to('ohm / ohm').magnitude

    def get_per_unit_current(self, I_act: Quantity) -> Union[float, complex]:
        """Returns the per-unit value of the actual current `I_act`."""
        I_pu: Quantity = I_act / self.I_base
        return I_pu.to('A / A').magnitude

    def get_per_unit_voltage(self, U_act: Quantity) -> Union[float, complex]:
        """Returns the per-unit value of the actual voltage `U_act`."""
        U_pu: Quantity = U_act / self.U_base
        return U_pu.to('V / V').magnitude

    def get_per_unit_power(self, S_act: Quantity) -> Union[float, complex]:
        """Returns the per-unit value of the actual power `S_act`."""
        S_pu: Quantity = S_act / self.S_base
        return S_pu.to('VA / VA').magnitude

    def get_actual_impedance(self, Z_pu: float) -> Quantity:
        Z_act = Z_pu * self.Z_base
        return Z_act

    # noinspection PyTypeHints
    def get_actual_current(self, I_pu: float | npt.ArrayLike) -> Quantity:
        I_act = I_pu * self.I_base
        return I_act

    # noinspection PyTypeHints
    def get_actual_voltage(self, U_pu: float | npt.ArrayLike) -> Quantity:
        U_act = U_pu * self.U_base
        return U_act

    def get_actual_power(self, S_pu: float) -> Quantity:
        S_act = S_pu * self.S_base
        return S_act


def convert_per_unit_impedance(
    Z_pu_source: Union[float, complex],
    source: PerUnitSystem,
    dest: PerUnitSystem
) -> Union[float, complex]:
    """
    Converts a per-unit impedance to another per-unit system.

    Parameters
    ----------
    Z_pu_source:
        Per-unit impedance.
    source:
        Per-unit system of `Z_pu`.
    dest:
        Per-unit system the per-unit impedance `Z_pu` is to be converted to.
    """
    Z_pu_dest = Z_pu_source * source.Z_base / dest.Z_base
    return Z_pu_dest.magnitude
