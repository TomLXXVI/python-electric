"""
Implements the sizing of electrical cables, i.e. the routines to determine the
cross-sectional area, current-carrying capacity, and nominal current of the
loaded conductors.

References
----------
Knockaert, J., Sergeant, P., Corne, B., Debruyne, C., Dereyne, S.,
Descheemaeker, J., Hespel, L., Vansteenberge, C., & Verhelst, B. (2015).
Laagspanningsinstallaties: technologie en ontwerp.
"""
import math
from enum import StrEnum, IntEnum
from dataclasses import dataclass, field

from ...utils.lookup_table import LookupTable
from ...materials import (
    ConductorMaterial,
    CONDUCTOR_PROPS,
    InsulationMaterial,
    INSULATION_PROPS
)
from .ampacity import unburied
from .ampacity import buried

__all__ = [
    "Ambient",
    "InstallMethod",
    "CableMounting",
    "CableArrangement",
    "CABLE_UNIT_REACTANCES",
    "CableData",
    "get_size",
    "set_size",
    "get_nominal_current",
]


def _create_unburied_group_correction_table() -> LookupTable:
    """
    Reduction factors for unburied groups of more than one circuit or multicore
    cable in a single layer or bundle, touching each other.
    """
    data = [
        [1.00, 0.80, 0.70, 0.65, 0.60, 0.57, 0.54, 0.52, 0.50, 0.45, 0.41, 0.38],
        [1.00, 0.85, 0.79, 0.75, 0.73, 0.72, 0.72, 0.71, 0.70, 0.70, 0.70, 0.70],
        [0.95, 0.81, 0.72, 0.68, 0.66, 0.64, 0.63, 0.62, 0.61, 0.61, 0.61, 0.61],
        [1.00, 0.88, 0.82, 0.77, 0.75, 0.73, 0.73, 0.72, 0.72, 0.72, 0.72, 0.72],
        [1.00, 0.87, 0.82, 0.80, 0.80, 0.79, 0.79, 0.78, 0.78, 0.78, 0.78, 0.78]
    ]
    col_header = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 20]  # number of circuits
    row_header = {  # installation method
        1: "Bunched in air, on a surface, embedded or enclosed",
        2: "Single layer on wall, floor, or unperforated tray",
        3: "Single layer fixed directly under a wooden ceiling",
        4: "Single layer on a perforated horizontal or vertical tray",
        5: "Single layer on ladder support or cleats, etc."
    }
    description = (
        "Reduction factors for unburied groups of more than one circuit or "
        "multicore cable in a single layer or bundle, touching each other."
    )
    cols_description = "number of circuits or multicore cables"
    rows_description = "method of installation"
    return LookupTable.create(
        row_header,
        col_header,
        data,
        description,
        rows_description,
        cols_description
    )

tbl_unburied_group_correction = _create_unburied_group_correction_table()


class CableArrangement(StrEnum):
    BUS_BAR = "bus_bar"
    MULTICORE = "multicore"
    SINGLE_SPACED = "single_spaced"
    SINGLE_SPACED_2R = "single_spaced_2r"
    SINGLE_SPACED_4R = "single_spaced_4r"
    SINGLE_TREFOIL = "single_trefoil"
    SINGLE_FLAT = "single_flat"


CABLE_UNIT_REACTANCES: dict[str, float] = {  # ohm / m
    'bus_bar': 0.15e-3,
    'multicore': 0.08e-3,
    'single_spaced': 0.15e-3,
    'single_spaced_2r': 0.145e-3,
    'single_spaced_4r': 0.19e-3,
    'single_trefoil': 0.085e-3,
    'single_flat': 0.095e-3
}


class InstallMethod(StrEnum):
    A1 = "A1"
    A2 = "A2"
    B1 = "B1"
    B2 = "B2"
    C = "C"
    D = "D"
    E = "E"
    F = "F"

    @property
    def unburied(self) -> bool:
        if self is not InstallMethod.D:
            return True
        return False
    
    @property
    def buried(self) -> bool:
        return not self.unburied


class Ambient(StrEnum):
    AIR = "air"
    GROUND = "ground"


class CableMounting(IntEnum):
    BUNCHED = 1  # Bunched in air, on a surface, embedded or enclosed
    UNPERFORATED_TRAY = 2  # Single layer on wall, floor, or unperforated tray
    CEILING = 3  # Single layer fixed directly under a wooden ceiling
    PERFORATED_TRAY = 4  # Single layer on a perforated horizontal or vertical tray
    LADDER = 5  # Single layer on ladder support or cleats, etc.


def get_nominal_current(I_load: float) -> float:
    nom_sizes = [1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0, 5.0, 6.3, 8.0]
    if 0 < I_load <= 10:
        I_load = round(I_load, 2)
        c = 1
    elif 10 < I_load <= 100:
        I_load = round(I_load / 10, 2)
        c = 10
    elif 100 < I_load <= 1000:
        I_load = round(I_load / 100, 2)
        c = 100
    elif 1000 < I_load <= 8_000:
        I_load = round(I_load / 1000, 2)
        c = 1000
    else:
        raise ValueError("Rated load current cannot exceed 8 kA.")
    nom_sizes.append(I_load)
    nom_sizes.sort()
    i = nom_sizes.index(I_load)
    if i == 0:
        I_nom = nom_sizes[1] * c
    elif i == len(nom_sizes) - 1:
        if c < 1000:
            I_nom = nom_sizes[0] * c * 10
        else:
            raise ValueError("Nominal current out of range.")
    else:
        I_nom = nom_sizes[i + 1] * c
    return I_nom


def _lookup_T_cmax(insulation: str, short_time: bool = False) -> float:
    """
    Returns the maximum allowable conductor temperature in °C depending on the
    given type of conductor insulation.

    Parameters
    ----------
    insulation: str, {"PVC", "rubber", "XLPE", "EPR"}
        Type of conductor insulation
    short_time: bool, default False
        Indicates short-time adiabatic heating of the conductor, i.e. at
        short-circuit conditions.

    Returns
    -------
    float

    Raises
    ------
    ValueError if the insulation type is unknown.
    """
    ins = INSULATION_PROPS.get(insulation.upper())
    if ins is None:
        raise ValueError(f"Unknown insulation type: {insulation}.")
    if short_time:
        return ins.T_max_short
    return ins.T_max_cont


def _temperature_correction(insulation: str, ambient: str, T_amb: float) -> float:
    """
    Returns the temperature correction factor with which the standard maximum
    continuous current must be multiplied for the given ambient temperature.

    Parameters
    ----------
    insulation: str, {"PVC", "rubber", "XLPE", "EPR"}
        Type of conductor insulation
    ambient: str, {"air", "ground"}
        Indicates whether the electrical cable is above or in the ground.
    T_amb: float
        Ambient temperature, degC or K

    Returns
    -------
    float
    """
    T_cmax = _lookup_T_cmax(insulation)
    T_amb_std = None
    match ambient.lower():
        case "air":
            T_amb_std = 30.0
        case "ground":
            T_amb_std = 20.0
    if T_amb_std is None:
        raise ValueError("Unknown ambient; must be either 'air' or 'ground'.")
    if T_amb >= T_cmax:
        raise ValueError(
            f"Ambient temperature cannot exceed the maximum "
            f"allowable conductor temperature {T_cmax} °C of "
            f"insulation type {insulation}."
        )
    k_T = math.sqrt((T_cmax - T_amb) / (T_cmax - T_amb_std))
    return round(k_T, 2)


def _harmonics_correction(harmonic3_content: float) -> float:
    """
    Returns the harmonic content correction factor for 3rd order harmonics in
    4-core and 5-core cables.
    """
    if not (0.0 <= harmonic3_content <= 1.0):
        raise ValueError("3rd harmonic fraction must be between 0.0 and 1.0")
    if 0.0 <= harmonic3_content < 0.15:
        k_h = 1.0
    elif 0.15 <= harmonic3_content < 0.33:
        k_h = 0.86
    elif 0.33 <= harmonic3_content < 0.45:
        # size selection is based on neutral current = 3x 3rd-harmonic phase current
        k_h = 0.86 / (3 * harmonic3_content)
    else:
        k_h = 1.0 / (3 * harmonic3_content)
    return k_h


def _joule_integral(conductor: str, insulation: str, S: float) -> float:
    """
    Returns the "Joule integral" in A².s of a conductor, i.e., the total thermal
    energy dissipated as heat during a short circuit, which is proportional to
    the square of the current and the time the fault persists. It is a critical
    value used to assess a conductor's ability to withstand the heat from a
    fault before its insulation gets damaged.

    Parameters
    ----------
    conductor: str, {"copper", "aluminium", "steel", "lead"}
        Indicates the type of conductor material
    insulation: str, {"PVC", "XLPE", "EPR", "rubber"}
        Indicates the type of insulation
    S:
        Cross-sectional area of 1 conductor, mm²

    Returns
    -------
    float
    """
    T_20 = 20.0  # °C
    cond = CONDUCTOR_PROPS[conductor]
    T_1 = _lookup_T_cmax(insulation)  # °C
    T_2 = _lookup_T_cmax(insulation, short_time=True)  # °C
    cp = cond.specific_heat  # J/(kg.K)
    d = cond.density  # kg/m³
    rho20 = cond.resistivity20 * 1e-6  # Ohm.m²/m = Ohm.m
    alp20 = cond.temp_coeff20  # 1/K
    k1 = cp * d / (rho20 * alp20)

    def f(T: float) -> float:
        return 1 + alp20 * (T - T_20)

    k2 = math.log(f(T_2) / f(T_1))
    k = math.sqrt(k1 * k2) * 1e-6
    ji = (k * S) ** 2
    return ji  # A².s


@dataclass
class CableData:
    """
    Dataclass used for sizing an electrical cable.

    Parameters
    ----------
    I_b: float
        The rated current through the cable in amperes.
    L: float
        Cable length in meters.
    conductor_type: ConductorMaterial
        Conductor material. See enum ConductorMaterial.
    insulation_type: InsulationMaterial
        Insulation material. See enum InsulationMaterial.
    ambient: Ambient
        The ambient where the cable is installed. See enum Ambient.
    T_amb: float.
        Steady-state operating temperature of the cable in degrees Celsius.
    install_method: InstallMethod
        The installation method of the cable. See enum InstallationMethod.
    cable_mounting: CableMounting
        Specifies in more detail the way cables are mounted. See enum
        CableMounting.
    cable_arrangement: CableArrangement
        Specifies in more detail how cables are arranged with respect to each
        other. See enum CableArrangement.
    num_circuits: int, default 1
        The number of circuits or multicore cables in the proximity of the
        cable.
    num_loaded_conductors: int, {2, 3 (default)}
        Number of loaded conductors in the cable/ single-core cables in the
        conduit.
    h3_content: float, default 0.0
        Fraction of third-order harmonics present in the total current,
        a number between 0.0 and 1.0.

    Attributes
    ----------
    S: float, default 0.0
        Cross-sectional area of the cable conductors in square millimeters.
    I_z: float, default 0.0
        Current-carrying capacity of the conductors in amperes at actual
        installation conditions.
    I_z0: float, default 0.0
        Current-carrying capacity of the conductors in amperes at standard
        installation conditions.
    I_n: float, default 0.0
        Standardized nominal current in amperes which is just greater than the
        rated load current of the cable, and that can be used for the
        current-protection device.
    k_z: float, default 0.0
        Total reduction factor to be applied to the current-carrying capacity of
        the cable at standard installation conditions.
    I2t: float, default 0.0
        Joule integral of the cable in A².s.
    """
    I_b: float # A
    L: float  # m
    conductor_type: ConductorMaterial
    insulation_type: InsulationMaterial
    ambient: Ambient
    T_amb: float  # degC
    install_method: InstallMethod
    cable_mounting: CableMounting
    cable_arrangement: CableArrangement
    num_circuits: int = 1
    num_loaded_conductors: int = 3
    h3_content: float = 0.0  # fraction 0..1

    S: float = field(init=False, default=0.0)  # mm²
    I_z: float = field(init=False, default=0.0)  # A
    I_z0: float = field(init=False, default=0.0)  # A
    I_n: float = field(init=False, default=0.0)  # A
    k_z: float = field(init=False, default=0.0)
    I2t: float = field(init=False, default=0.0)  # A².s

    def __post_init__(self):
        self.cond_props = CONDUCTOR_PROPS[self.conductor_type]
        self.insul_props = INSULATION_PROPS[self.insulation_type]
    
    def update(self, k, I_z0, S, I_n) -> None:
        I_z = k * I_z0  # current-carrying capacity at actual conditions
        self.k_z = k
        self.S = S
        self.I_z = I_z
        self.I_z0 = I_z0
        self.I_n = I_n
        self.I2t = _joule_integral(self.cond_props.type, self.insul_props.type, S)


def get_size(
    cable_data: CableData,
    based_on_I_nom: bool = True
) -> CableData:
    """
    Determines the required minimal cross-sectional area S of the cable
    conductors in square millimeters (mm²).

    Parameters
    ----------
    cable_data:
        Instance of `CableSizingData` holding data about the cable needed to
        determine the required minimal cross-sectional area of the cable
        conductors.
    based_on_I_nom: bool, default True
        Indicates that the cross-sectional area of the cable conductors needs to
        be determined based on the next standardized nominal current which is
        just greater than the rated load current. If False, conductors are sized
        based on the rated load current of the cable.
    """
    k_T = _temperature_correction(
        cable_data.insul_props.type,
        cable_data.ambient,
        cable_data.T_amb
    )
    k_G = tbl_unburied_group_correction.data_value(cable_data.cable_mounting, cable_data.num_circuits)
    k_H = _harmonics_correction(cable_data.h3_content)
    k = k_T * k_G * k_H
    I_b = cable_data.I_b
    I_b_corr = I_b / k
    I_n = get_nominal_current(I_b)
    I_n_corr = I_n / k  # refer I_n to standard conditions
    if cable_data.install_method.unburied:
        S = unburied.get_cross_sectional_area(
            conductor_props=cable_data.cond_props,
            insulation_props=cable_data.insul_props,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            ref_method=cable_data.install_method,
            current=I_n_corr if based_on_I_nom else I_b_corr
        )
        I_z0 = unburied.get_ampacity(
            conductor_props=cable_data.cond_props,
            insulation_props=cable_data.insul_props,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            ref_method=cable_data.install_method,
            cross_sectional_area=S
        )
    elif cable_data.install_method.buried:
        S = buried.get_cross_sectional_area(
            conductor_props=cable_data.cond_props,
            insulation_props=cable_data.insul_props,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            current=I_n_corr if based_on_I_nom else I_b_corr
        )
        I_z0 = buried.get_ampacity(
            conductor_props=cable_data.cond_props,
            insulation_props=cable_data.insul_props,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            cross_sectional_area=S
        )
    else:
        raise ValueError(f"Unknown installation method.")
    cable_data.update(k, I_z0, S, I_n)
    return cable_data


def set_size(cable_data: CableData, S: float) -> CableData:
    """
    Determines the nominal current, the ampacity, and the Joule-integral of the
    cable when the cross-sectional area of the cable conductors is given.

    Parameters
    ----------
    cable_data:
        Instance of `CableSizingData` holding the data about the conductor
        material, insulation material, and the installation method of the cable.
    S:
        Cross-sectional area of the cable conductors in square millimeters.

    Returns
    -------
    CableData
    """
    k_T = _temperature_correction(
        cable_data.insul_props.type,
        cable_data.ambient,
        cable_data.T_amb
    )
    k_G = tbl_unburied_group_correction.data_value(cable_data.cable_mounting, cable_data.num_circuits)
    k_H = _harmonics_correction(cable_data.h3_content)
    k = k_T * k_G * k_H
    I_b = cable_data.I_b
    I_n = get_nominal_current(I_b)
    if cable_data.install_method.unburied:
        I_z0 = unburied.get_ampacity(
            conductor_props=cable_data.cond_props,
            insulation_props=cable_data.insul_props,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            ref_method=cable_data.install_method,
            cross_sectional_area=S
        )
    elif cable_data.install_method.buried:
        I_z0 = buried.get_ampacity(
            conductor_props=cable_data.cond_props,
            insulation_props=cable_data.insul_props,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            cross_sectional_area=S
        )
    else:
        raise ValueError(f"Unknown installation method.")
    cable_data.update(k, I_z0, S, I_n)
    return cable_data
