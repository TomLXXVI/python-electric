import math

from .. import Quantity
from ..general import VoltReference

Q_ = Quantity

__all__ = [
    "get_3ph_line_current",
    "get_1ph_line_current",
    "voltage_drop"
]


def get_3ph_line_current(
    P: float | Quantity,
    U_line: float | Quantity,
    cos_phi: float = 1.0
) -> float:
    if isinstance(P, Quantity):
        P = P.to('W').m
    if isinstance(U_line, Quantity):
        U_line = U_line.to('V').m
    I = P / (math.sqrt(3) * U_line * cos_phi)
    return I


def get_1ph_line_current(
    P: float | Quantity,
    U_phase: float | Quantity,
    cos_phi: float = 1.0
) -> float:
    if isinstance(P, Quantity):
        P = P.to('W').m
    if isinstance(U_phase, Quantity):
        U_phase = U_phase.to('V').m
    I = P / (U_phase * cos_phi)
    return I


def voltage_drop(
    R: float | Quantity,
    X: float | Quantity,
    I: float | Quantity,
    volt_ref: VoltReference,
    cos_phi: float = 0.8
) -> float:
    """
    Returns the voltage drop in volts.

    Parameters
    ----------
    R: float | Quantity
        Resistance of cable conductor in ohms (if float).
    X: float | Quantity
        Reactance of cable conductor in ohms (if float).
    I: float | Quantity
        Current through the cable in amperes (if float).
    volt_ref: VoltReference
        Indicates the reference against which the voltage is measured. See enum
        VoltageRef.
    cos_phi: float, default 0.8
        Power factor of the current.

    Returns
    -------
    float
    """
    if isinstance(R, Quantity):
        R = R.to('ohm').m
    if isinstance(X, Quantity):
        X = X.to('ohm').m
    if isinstance(I, Quantity):
        I = I.to('A').m
    sin_phi = math.sqrt(1.0 - cos_phi ** 2)
    U_drop = volt_ref.value * I * (R * cos_phi + X * sin_phi)
    return U_drop


def delta_to_star(
    Z_a: Quantity,
    Z_b: Quantity,
    Z_c: Quantity
) -> tuple[Quantity, Quantity, Quantity]:
    """
    Delta-to-star transformation.
    """
    Z_a = Z_a.to('ohm').m
    Z_b = Z_b.to('ohm').m
    Z_c = Z_c.to('ohm').m
    den = Z_a + Z_b + Z_c
    Z_alpha = Z_b * Z_c / den
    Z_beta = Z_a * Z_c / den
    Z_gamma = Z_a * Z_b / den
    return Q_(Z_alpha, 'ohm'), Q_(Z_beta, 'ohm'), Q_(Z_gamma, 'ohm')


def star_to_delta(
    Z_alpha: Quantity,
    Z_beta: Quantity,
    Z_gamma: Quantity
) -> tuple[Quantity, Quantity, Quantity]:
    """
    Star-to-delta transformation.
    """
    Z_alpha = Z_alpha.to('ohm').m
    Z_beta = Z_beta.to('ohm').m
    Z_gamma = Z_gamma.to('ohm').m
    num = Z_alpha * Z_beta + Z_beta * Z_gamma + Z_alpha * Z_gamma
    Z_a = num / Z_alpha
    Z_b = num / Z_beta
    Z_c = num / Z_gamma
    return Q_(Z_a, 'ohm'), Q_(Z_b, 'ohm'), Q_(Z_c, 'ohm')
