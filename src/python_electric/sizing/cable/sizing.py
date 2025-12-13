"""
Implementation of the sizing routine for electrical cables, i.e. to determine 
the cross-sectional area, current-carrying capacity, and nominal current of the 
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
    ConductorMaterials,
    CONDUCTOR_MATERIALS,
    InsulationMaterials,
    INSULATION_MATERIALS
)
from .ampacity import unburied
from .ampacity import buried

__all__ = [
    "Ambient",
    "InstallationMethods",
    "CableMounting",
    "CableArrangement",
    "CABLE_UNIT_REACTANCES",
    "CableSizingData",
    "get_cable_sizing",
    "set_cable_sizing",
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


class InstallationMethods(StrEnum):
    A1 = "A1"
    A2 = "A2"
    B1 = "B1"
    B2 = "B2"
    C = "C"
    D = "D"
    E = "E"
    F = "F"


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
    ins = INSULATION_MATERIALS.get(insulation.upper())
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
    cond = CONDUCTOR_MATERIALS[conductor]
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
class CableSizingData:
    """
    Class for sizing electrical cables in low-voltage systems.

    After instantiation, call method `size()` to size the cable. After a
    cable has been sized, the voltage drop across the cable can be calculated
    with method `voltage_drop()`.

    Parameters
    ----------
    rated_load_current: float
        The rated current through the cable in amperes.
    length: float
        Cable length in meters.
    conductor_material: ConductorMaterials
        Conductor material. See enum ConductorMaterials.
    insulation_material: InsulationMaterials
        Insulation material. See enum InsulationMaterials.
    ambient: Ambient
        The ambient where the cable is installed. See enum Ambient.
    amb_temperature: float.
        Steady-state operating temperature of the cable in degrees Celsius.
    install_method: InstallationMethods
        The installation method of the cable. See enum InstallationMethods.
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
    harmonic3_content: float, default 0.0
        Fraction of third-order harmonics present in the total current,
        a number between 0.0 and 1.0.

    Attributes
    ----------
    cross_section_area: float, default 0.0
        Cross-sectional area of the conductors in square millimeters, calculated
        by calling function size().
    current_capacity: float, default 0.0
        Current-carrying capacity of the conductors in amperes at actual
        installation conditions, calculated by calling function size().
    current_capacity_ref: float, default 0.0
        Current-carrying capacity of the conductors in amperes at standard
        installation conditions, calculated by calling function size().
    nom_current: float, default 0.0
        Standardized nominal current in amperes which is just greater than the
        rated load current of the cable, and that can be used for the
        current-protection device. Calculated by calling function size().
    sizing_factor: float, default 0.0
        Total reduction factor to be applied to the current-carrying capacity of
        the cable at standard installation conditions, determined by calling
        function size().
    joule_integral: float, default 0.0
        Joule integral of the cable in A².s, determined by calling function
        size()
    """
    rated_load_current: float # A
    length: float  # m
    conductor_material: ConductorMaterials
    insulation_material: InsulationMaterials
    ambient: Ambient
    amb_temperature: float  # degC
    install_method: InstallationMethods
    cable_mounting: CableMounting
    cable_arrangement: CableArrangement
    num_circuits: int = 1
    num_loaded_conductors: int = 3
    harmonic3_content: float = 0.0  # fraction 0..1

    cross_section_area: float = field(init=False, default=0.0)  # mm²
    current_capacity: float = field(init=False, default=0.0)  # A
    current_capacity_ref: float = field(init=False, default=0.0)  # A
    nom_current: float = field(init=False, default=0.0)  # A
    sizing_factor: float = field(init=False, default=0.0)
    joule_integral: float = field(init=False, default=0.0)  # A².s

    def __post_init__(self):
        self.conductor = CONDUCTOR_MATERIALS[self.conductor_material]
        self.insulation = INSULATION_MATERIALS[self.insulation_material]
        self._install_methods_unburied = (
            InstallationMethods.A1,
            InstallationMethods.A2,
            InstallationMethods.B1,
            InstallationMethods.B2,
            InstallationMethods.C,
            InstallationMethods.E,
            InstallationMethods.F
        )
        self._install_methods_buried = (
            InstallationMethods.D,
        )


def get_cable_sizing(cable_data: CableSizingData, based_on_I_nom: bool = True) -> CableSizingData:
    """
    Sizes the cable. Returns the required cross-sectional area of the
    conductors in square millimeters.

    Parameters
    ----------
    cable_data:
        Instance of `CableData` holding the data that is needed to size the
        cable.
    based_on_I_nom: bool, default True
        Indicates that the determination of the cross-sectional area of the
        conductors should be based on the next standardized nominal current
        greater than the rated load current. If False, conductors are sized
        based on the rated load current.
    """
    k_T = _temperature_correction(
        cable_data.insulation.type,
        cable_data.ambient,
        cable_data.amb_temperature
    )
    k_G = tbl_unburied_group_correction.data_value(cable_data.cable_mounting, cable_data.num_circuits)
    k_H = _harmonics_correction(cable_data.harmonic3_content)
    k = k_T * k_G * k_H
    I_b = cable_data.rated_load_current
    I_b_corr = I_b / k
    I_n = get_nominal_current(I_b)
    I_n_corr = I_n / k  # refer I_n to standard conditions
    if cable_data.install_method in cable_data._install_methods_unburied:
        S = unburied.get_cross_sectional_area(
            conductor_material=cable_data.conductor,
            insulation=cable_data.insulation,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            ref_method=cable_data.install_method,
            current=I_n_corr if based_on_I_nom else I_b_corr
        )
        I_z0 = unburied.get_ampacity(
            conductor_material=cable_data.conductor,
            insulation=cable_data.insulation,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            ref_method=cable_data.install_method,
            cross_sectional_area=S
        )
    elif cable_data.install_method in cable_data._install_methods_buried:
        S = buried.get_cross_sectional_area(
            conductor_material=cable_data.conductor,
            insulation=cable_data.insulation,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            current=I_n_corr if based_on_I_nom else I_b_corr
        )
        I_z0 = buried.get_ampacity(
            conductor_material=cable_data.conductor,
            insulation=cable_data.insulation,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            cross_sectional_area=S
        )
    else:
        raise ValueError(f"Unknown installation method.")
    I_z = k * I_z0  # current-carrying capacity at actual conditions

    cable_data.sizing_factor = k
    cable_data.cross_section_area = S
    cable_data.current_capacity = I_z
    cable_data.current_capacity_ref = I_z0
    cable_data.nom_current = I_n
    cable_data.joule_integral = _joule_integral(
        cable_data.conductor.type,
        cable_data.insulation.type,
        S
    )
    return cable_data


def set_cable_sizing(cable_data: CableSizingData, S: float) -> CableSizingData:
    """
    Determines the nominal current, ampacity, and Joule-integral of the cable
    whose conductor cross-sectional area is given.

    Parameters
    ----------
    cable_data:
        Instance of `CableData` holding data about the conductor material,
        the insulation material, and the installation of the cable.
    S:
        Cross-sectional area of the cable conductors in square millimeters.

    Returns
    -------
    CableSizingData
    """
    k_T = _temperature_correction(
        cable_data.insulation.type,
        cable_data.ambient,
        cable_data.amb_temperature
    )
    k_G = tbl_unburied_group_correction.data_value(cable_data.cable_mounting, cable_data.num_circuits)
    k_H = _harmonics_correction(cable_data.harmonic3_content)
    k = k_T * k_G * k_H
    I_b = cable_data.rated_load_current
    I_n = get_nominal_current(I_b)
    if cable_data.install_method in cable_data._install_methods_unburied:
        I_z0 = unburied.get_ampacity(
            conductor_material=cable_data.conductor,
            insulation=cable_data.insulation,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            ref_method=cable_data.install_method,
            cross_sectional_area=S
        )
    elif cable_data.install_method in cable_data._install_methods_buried:
        I_z0 = buried.get_ampacity(
            conductor_material=cable_data.conductor,
            insulation=cable_data.insulation,
            num_loaded_conductors=cable_data.num_loaded_conductors,
            cross_sectional_area=S
        )
    else:
        raise ValueError(f"Unknown installation method.")
    I_z = k * I_z0  # current-carrying capacity at actual conditions

    cable_data.sizing_factor = k
    cable_data.cross_section_area = S
    cable_data.current_capacity = I_z
    cable_data.current_capacity_ref = I_z0
    cable_data.nom_current = I_n
    cable_data.joule_integral = _joule_integral(
        cable_data.conductor.type,
        cable_data.insulation.type,
        S
    )
    return cable_data
