import math

from .. import Quantity, Q_, VoltReference, PhaseSystem


__all__ = [
    "load_current",
    "voltage_drop",
    "delta_to_star",
    "star_to_delta",
    "R_cond"
]


def load_current(
    U_l: Quantity,
    cos_phi: float,
    P_e: Quantity | None = None,
    P_m: Quantity | None = None,
    eta: Quantity = Q_(100, 'pct'),
    phase_system: PhaseSystem = PhaseSystem.PH3
) -> Quantity:
    """
    Returns the load current when load power is given (either electrical power
    or mechanical power).

    Parameters
    ----------
    U_l: Quantity
        Line-to-line voltage in case of a three-phase network, or neutral-to-line
        voltage in case of a single-phase network.
    cos_phi: float
        Power factor.
    P_e: Quantity, optional
        Electrical power.
    P_m: Quantity, optional
        Mechanical power.
    eta: Quantity, optional
        Electromechanical efficiency.
    phase_system:
        Indicates if the electrical network is either three- or a single-phase.

    Returns
    -------
    Quantity
    """
    cP = phase_system.cP()
    if P_e is not None:
        I_b = P_e / (cP * U_l * cos_phi)
        return I_b.to('A')
    elif P_m is not None:
        P_e = P_m / eta
        I_b = P_e / (cP * U_l * cos_phi)
        return I_b.to('A')
    else:
        raise ValueError(f"Either 'P_e', or 'P_m' and 'eta' must be specified.")


def voltage_drop(
    R: Quantity,
    X: Quantity,
    I: Quantity,
    volt_ref: VoltReference,
    cos_phi: float = 0.8
) -> Quantity:
    """
    Returns voltage drop.

    Parameters
    ----------
    R: Quantity
        Resistance of cable conductors.
    X: Quantity
        Reactance of cable conductors.
    I: Quantity
        Load current through cable.
    volt_ref: VoltReference
        Reference of voltage measurement. See enum `VoltReference`.
    cos_phi: float, default 0.8
        Power factor of the load current.

    Returns
    -------
    Quantity
    """
    R = R.to('ohm').m
    X = X.to('ohm').m
    I = I.to('A').m
    sin_phi = math.sqrt(1.0 - cos_phi ** 2)
    U_drop = volt_ref.value * I * (R * cos_phi + X * sin_phi)
    return Q_(U_drop, 'V')


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


def R_cond(rho: Quantity, L: Quantity, S: Quantity) -> Quantity:
    """
    Returns the ohmic resistance of a conductor.

    Parameters
    ----------
    rho: Quantity
        Resistivity of conductor material.
    L: Quantity
        Conductor length.
    S: Quantity
        Cross-sectional area of conductor.

    Returns
    -------
    Quantity
    """
    R = rho * L / S
    return R.to('ohm')
