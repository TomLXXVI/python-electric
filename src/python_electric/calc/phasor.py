import math

__all__ = [
    "phasor",
    "phasor_rad",
    "polar",
    "phasor_from_cos",
    "cosphi"
]


def phasor(magnitude: float, angle_deg: float) -> complex:
    """Return a complex number for a phasor (angle in degrees)."""
    rad = math.radians(angle_deg)
    return magnitude * (math.cos(rad) + 1j * math.sin(rad))


def phasor_rad(magnitude: float, angle_rad: float) -> complex:
    """Return a complex number for a phasor (angle in radians)."""
    return magnitude * (math.cos(angle_rad) + 1j * math.sin(angle_rad))


def polar(z: complex) -> tuple[float, float]:
    """Return (magnitude, angle_deg)."""
    magnitude = abs(z)
    angle_deg = math.degrees(math.atan2(z.imag, z.real))
    return magnitude, angle_deg


def phasor_from_cos(magnitude: float, cos_phi: float, *, lagging=True) -> complex:
    """
    Create a phasor when only magnitude and cos(phi) are given.
    lagging=True  negative angle (inductive, current behind voltage)
    lagging=False positive angle (capacitive, current ahead of voltage)
    """
    phi = math.acos(cos_phi)
    angle = -phi if lagging else +phi
    return phasor_rad(magnitude, angle)


def cosphi(z: complex) -> float:
    """Return cos(phi) for a complex phasor."""
    return math.cos(math.atan2(z.imag, z.real))
