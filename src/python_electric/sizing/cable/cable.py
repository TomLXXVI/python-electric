"""
Cable sizing and impedance calculation.

Integrates the sizing routine and impedance calculation of a cable inside a
single class `Cable`.

References
----------
Knockaert, J., Sergeant, P., Corne, B., Debruyne, C., Dereyne, S.,
Descheemaeker, J., Hespel, L., Vansteenberge, C., & Verhelst, B. (2015).
Laagspanningsinstallaties: technologie en ontwerp.
"""
import math
from copy import deepcopy
from dataclasses import dataclass, field
from typing import ClassVar
from abc import ABC

import numpy as np

from ... import Quantity
from ...utils.charts import LineChart
from ...general import (
    PhaseSystem,
    VoltReference
)
from ...materials import (
    ConductorMaterials,
    InsulationMaterials
)
from ...calc import voltage_drop
from ...protection import (
    CircuitBreaker,
    check_current_based_selectivity,
)
from ...protection import earthing_system
from ...protection.earthing_system import (
    EarthingSystem,
    IndirectContactProtectionResult,
)
from . import sizing
from .sizing import (
    Ambient,
    InstallationMethods,
    CableMounting,
    CableArrangement,
    CableSizingData
)

__all__ = [
    "ThreePhaseCable",
    "SinglePhaseCable",
    "AbstractCable",
    "plot_cable_characteristic",
    "check_selectivity"
]

Q_ = Quantity


@dataclass
class AbstractCable(ABC):
    """
    Class for sizing electrical cables in low-voltage three-phase systems.

    The cable object is sized immediately on instantiation of the class.
    If a cross-sectional area is already specified by the user, it will be
    checked if the conductors can carry the specified load current. If not the
    case, a `ValueError` exception will be raised.

    Parameters
    ----------
    P_elec: Quantity | None = None
        Rated active electric power that flows through the cable.
    P_mech: Quantity | None = None
        Rated mechanical power of the machine fed through the cable.
    eta: float = 1.0
        Rated electromechanical conversion efficiency of the connected machine.
    I_load: Quantity | None = None
        Rated load current through the cable.
    cos_phi: float = 0.8
        Rated power factor, i.e. cosine of the phase shift between voltage and
        rated load current.
    U_line: Quantity = Q_(400, 'V')
        Line-to-line voltage in case of a three-phase system or neutral-to-line
        voltage in case of a single-phase system.
    k_simul: float = 1.0
        Simultaneity factor to be applied for (industrial) distribution boards.
            Number of circuits  k_simul
            2 to 3              0.9
            4 to 5              0.8
            6 to 9              0.7
            >= 10               0.6
    k_ext: float = 1.0
        Expansion factor that takes into account possible future expansions.
        For new, industrial installations, a guideline value is 1.2. 
    L: Quantity = Q_(10, 'm')
        Cable length.
    conductor_material: ConductorMaterials = ConductorMaterials.COPPER
        Material the cable conductors are made of.
    insulation_material: InsulationMaterials = InsulationMaterials.XLPE
        Insulation material around the conductors of the cable.
    ambient: Ambient = Ambient.AIR
        Ambient where the cable is installed, either in air or buried in the
        ground.
    T_amb: Quantity = Q_(35, 'degC')
        Ambient temperature of the cable.
    install_method: InstallationMethod = InstallationMethod.E
        Installation method according to enum InstallationMethods.
    cable_mounting: CableMounting = CableMounting.PERFORATED_TRAY
        Specifies in more detail the way cables are mounted. See enum
        CableMounting.
    cable_arrangement: CableArrangement = CableArrangement.MULTICORE
        Specifies in more detail how cables are arranged with respect to each
        other. See enum CableArrangement.
    num_circuits: int = 1
        The number of circuits or multicore cables in the proximity of the
        cable.
    harmonic3_content: float = 0.0
        Fraction of third-order harmonics present in the total load current;
        a number between 0.0 and 1.0.
    sizing_based_on_I_nom: bool = False
        Indicates that the determination of the cross-sectional area of the
        conductors should be based on the standardized nominal current which is
        just greater than the rated load current. If False, conductors are sized
        based on rated load current.
    name: str, optional
        To uniquely identify the cable object.
    S: Quantity, optional
        Cross-sectional area of the loaded cable conductors. By default, this
        cross-sectional area is set to None, meaning that a valid
        cross-sectional area is to be determined on instantiation of the cable
        object.
    S_pe: Quantity, optional
        Cross-sectional area of the protective earth conductor. By default, this
        cross-sectional area is set equal to the cross-sectional area of the
        loaded cable conductors.
    earthing_system: EarthingScheme, default E
        Indicates the earthing system of the low-voltage system. By default,
        this is set to the TN-earthing system.

    Attributes
    ----------
    I_z: Quantity
        Current-carrying capacity (ampacity) of the conductors at actual
        installation conditions.
    I_nom: Quantity
        Standardized nominal current which is just greater than the rated load
        current of the cable. This current can be used to select the
        current-protective device for the cable.
    I_load_corr: Quantity
        Load current corrected by multiplication with the simultaneity factor
        and expansion factor. This current is used to size the cable.
    Z: dict[str, dict[int, Quantity]]
        Dictionary that contains the zero-sequence (key 0), positive sequence
        (key 1), and negative sequence (key 2) impedances of the cable
        conductors at three different temperatures:
        -   at 20°C (key "T_20") used for calculating maximum short-circuit
            currents.
        -   at the maximum continuous temperature of the insulation
            (key "T_nom").
        -   at 150°C (key "T_150") used for calculating minimum short-circuit
            currents.
        E.g. to get Z1 at 150 °C, use Z["T_150"][1].
    U_drop: Quantity
        Absolute drop in voltage across the cable in volts. Available after
        method `voltage_drop()` has been called.
    U_drop_pct: Quantity
        Relative drop in voltage in percent. Available after method
        `voltage_drop()` has been called.
    joule_integral: Quantity
        Joule integral of the cable.
    circuit_breaker: CircuitBreaker
        Circuit breaker connected to the cable. Available after method
        `connect_circuit_breaker()` has been called.
    I_sc_max: Quantity
        Calculated maximum short-circuit current at the location of the
        circuit breaker (most often caused by a three-phase fault).
    I_sc_min: Quantity
        Calculated minimum short-circuit current at the end of the cable
        (usually caused by a single line-to-ground fault).

    Methods
    -------
    voltage_drop:
        Calculates the voltage drop across the cable.
    connect_circuit_breaker:
        Connects a circuit breaker to the cable.
    check_TN_indirect_contact_protection:
        Can be used for low-voltage systems having a TN-earthing scheme.
        Returns the maximum allowable interruption time in order to protect
        against indirect contact in case of an insulation fault. If a circuit
        breaker is connected to the cable, also returns the maximum length of
        the cable that can be protected by this circuit breaker.
    """
    phase_system: ClassVar[PhaseSystem]
    P_elec: Quantity | None = None
    P_mech: Quantity | None = None
    eta: float = 1.0
    I_load: Quantity | None = None
    cos_phi: float = 0.8
    U_line: Quantity = Q_(400, 'V')
    k_simul: float = 1.0
    k_ext: float = 1.0
    L: Quantity = Q_(10, 'm')
    conductor_material: ConductorMaterials = ConductorMaterials.COPPER
    insulation_material: InsulationMaterials = InsulationMaterials.XLPE
    ambient: Ambient = Ambient.AIR
    T_amb: Quantity = Q_(30, 'degC')
    install_method: InstallationMethods = InstallationMethods.E
    cable_mounting: CableMounting = CableMounting.PERFORATED_TRAY
    cable_arrangement: CableArrangement = CableArrangement.MULTICORE
    num_circuits: int = 1
    harmonic3_content: float = 0.0
    sizing_based_on_I_nom: bool = False
    name: str = "cable"
    S: Quantity | None = None
    S_pe: Quantity | None = None
    earthing_system: EarthingSystem = EarthingSystem.TN

    U_drop: Quantity = field(init=False)
    U_drop_pct: Quantity = field(init=False)
    circuit_breaker: CircuitBreaker = field(init=False, default=None)
    I_sc_max: Quantity = field(init=False, default=None)
    I_sc_min: Quantity = field(init=False, default=None)
    U_phase: Quantity = field(init=False, default=None)

    def __post_init__(self):
        self._calc_load_current()

        if self.S is None:
            self._calculate_cable_sizing()
        else:
            S_user = deepcopy(self.S)
            self._calculate_cable_sizing()
            if self.S.to('mm ** 2') > S_user.to('mm ** 2'):
                raise ValueError(
                    f"The given csa {S_user.to('mm ** 2'):~P.2f} is too small."
                    f"A csa of at least {self.S.to('mm ** 2'):~P.2f} is "
                    f"required."
                )
            else:
                self._recalculate_cable_sizing(S_user)

        if self.S_pe is None:
            self.S_pe = self.S

        if self.phase_system.three_phase:
            self.U_phase = self.U_line / math.sqrt(3)
        else:
            self.U_phase = self.U_line

    def _calc_load_current(self) -> None:
        cP = self.phase_system.cP()
        if self.P_elec is not None:
            self.I_load = self.P_elec / (cP * self.U_line * self.cos_phi)
        elif self.P_mech is not None:
            self.P_elec = self.P_mech / self.eta
            self.I_load = self.P_elec / (cP * self.U_line * self.cos_phi)
        elif self.I_load is not None:
            pass
        else:
            raise ValueError(
                "Either `P_elec`, `P_mech` (and `eta`) or `I_b` must be "
                "specified."
            )
        self.I_load_corr = self.I_load * self.k_simul * self.k_ext

    def _calculate_cable_sizing(self) -> None:
        if self.phase_system.three_phase:
            self.num_loaded_conductors: int = 3
        else:
            self.num_loaded_conductors: int = 2
        sizing_data = CableSizingData(
            rated_load_current=self.I_load_corr.to('A').m,
            length=self.L.to('m').m,
            conductor_material=self.conductor_material,
            insulation_material=self.insulation_material,
            ambient=self.ambient,
            amb_temperature=self.T_amb.to('degC').m,
            install_method=self.install_method,
            cable_mounting=self.cable_mounting,
            cable_arrangement=self.cable_arrangement,
            num_circuits=self.num_circuits,
            num_loaded_conductors=self.num_loaded_conductors,
            harmonic3_content=self.harmonic3_content,
        )
        sizing_data = sizing.get_cable_sizing(sizing_data, self.sizing_based_on_I_nom)
        self._set_attributes(sizing_data)
        return None

    def _recalculate_cable_sizing(self, S: Quantity) -> None:
        data = sizing.set_cable_sizing(self.sizing_data, S.to('mm ** 2').m)
        self._set_attributes(data)
        return None

    def _set_attributes(self, sizing_data: CableSizingData) -> None:
        self.sizing_data = sizing_data
        self.S = Q_(self.sizing_data.cross_section_area, 'mm ** 2')
        self.I_z = Q_(self.sizing_data.current_capacity, 'A')
        self.I_nom = Q_(self.sizing_data.nom_current, 'A')
        self.joule_integral = Q_(self.sizing_data.joule_integral, 'A ** 2 * s')
        T_nom = self.sizing_data.insulation.T_max_cont
        self.Z: dict[str, dict[int, Quantity]] = {
            "T_20": self._calc_impedance(Q_(20, 'degC')),
            "T_nom": self._calc_impedance(Q_(T_nom, 'degC')),
            "T_150": self._calc_impedance(Q_(150, 'degC'))
        }
        return None

    def _calc_impedance(self, T: Quantity) -> dict[int, Quantity]:
        from python_electric.equipment import Cable
        cable = Cable(
            length=self.L,
            cross_section_area=self.S,
            conductor_material=ConductorMaterials.COPPER,
            cable_arrangement=CableArrangement.MULTICORE,
            temperature=T
        )
        return {1: cable.Z1, 2: cable.Z2, 0: cable.Z0}

    @property
    def has_circuit_breaker(self) -> bool:
        return isinstance(self.circuit_breaker, CircuitBreaker)

    def voltage_drop(
        self,
        U_line: Quantity | None = None,
        I_load: Quantity | None = None,
        cos_phi: float | None = None,
        volt_ref: VoltReference | None = None
    ) -> None:
        """
        Calculates the absolute voltage drop in volts and the relative voltage
        drop in percent of the given input voltage. Results are stored in
        attributes `U_drop` and `U_drop_pct`.

        Parameters
        ----------
        U_line: Quantity, optional
            Line-to-line voltage at the entry in case of a three-phase system,
            or ground-to-line voltage at the entry in case of a single-phase
            system.
        I_load: Quantity, optional
            Load current through the cable.
        cos_phi: float, optional
            Power factor of the current, i.e. the cosine of the phase shift
            between voltage and current.
        volt_ref: VoltReference, optional
            In case of a three-phase cable, if `volt_ref` is set to
            `VoltReference.PH3_LINE_TO_LINE`, the absolute voltage drop is
            measured between two phase conductors, and the relative voltage drop
            is referred to the line-to-line voltage. If `volt_ref` is set to
            `VoltReference.PH3_GROUND_TO_LINE`, the voltage drop is measured
            across a single line conductor with reference to ground, and the
            relative voltage drop is referred to the ground-to-line voltage.
            If `volt_ref` is left as None, `VoltReference.PH3_GROUND_TO_LINE`
            will be assumed.
            In case of a single-phase cable, parameter `volt_ref` has no real
            use and can be left as None.

        Returns
        -------
        None
        """
        if U_line is None:
            U_line = self.U_line
        if I_load is None:
            I_load = self.I_load
        if cos_phi is None:
            cos_phi = self.cos_phi

        if self.phase_system.three_phase and volt_ref is None:
            volt_ref = VoltReference.PH3_GROUND_TO_LINE
        elif self.phase_system.single_phase and volt_ref is None:
            volt_ref = VoltReference.PH1

        R = self.Z["T_nom"][1].real
        X = self.Z["T_nom"][1].imag

        # noinspection PyTypeChecker
        U_drop = voltage_drop(R, X, I_load, volt_ref, cos_phi)

        U = U_line
        if self.phase_system.three_phase and volt_ref == VoltReference.PH3_GROUND_TO_LINE:
            U = U_line / math.sqrt(3)
        U_drop_pct = 100 * U_drop / U.to('V').m

        self.U_drop = Q_(U_drop, 'V')
        self.U_drop_pct = Q_(U_drop_pct, 'pct')
        return None

    def connect_circuit_breaker(
        self,
        standard: CircuitBreaker.Standard,
        category: CircuitBreaker.Category,
        ultim_break_capacity: Quantity,
        I_sc_max: Quantity,
        I_sc_min: Quantity,
        k_magn_trip: float | None = None,
        E_through: Quantity | None = None,
        t_m_lim: Quantity | None = None
    ) -> None:
        """
        Connects a circuit breaker with the bus bars.

        Parameters
        ----------
        standard: CircuitBreaker.Standard
            Either the residential standard (IEC 60898-1) or the industrial
            standard (IEC 60947-2), which specifies the minimal performance
            requirements of the circuit breaker. See enum
            CircuitBreaker.Standard.
        category: CircuitBreaker.Category
            Specifies the category of the circuit breaker depending on its
            magnetic trip threshold. Either category B, C, D, or AJUSTABLE. See
            enum CircuitBreaker.Category. Note that category ADJUSTABLE also
            demands that the circuit breaker is of the industrial type.
        ultim_break_capacity: Quantity
            Ultimate breaking capacity of the circuit breaker.
        I_sc_max: Quantity
            Calculated maximum short-circuit current at the location of the
            circuit breaker (most often caused by a three-phase fault).
        I_sc_min: Quantity
            Calculated minimum short-circuit current at the end of the cable
            (usually caused by a single line-to-ground fault).
        k_magn_trip: float, optional
            Multiplication factor that determines the rated magnetic trip
            current as a multiple of the thermal current setting if the circuit
            breaker is of the industrial type and adjustable.
        E_through: Quantity, optional
            Thermal energy (I²t) let through by the circuit breaker at the
            calculated maximum short-circuit current.
        t_m_lim: Quantity, optional
            Upper limit of instantaneous magnetic tripping time with regard to
            short-circuits. By default, this time limit is set to 100 ms
            according to the residential standard IEC 60898-1.
        """
        cb = CircuitBreaker(
            standard=standard,
            category=category,
            load_current=self.I_load,
            nom_current=self.I_nom,
            ampacity=self.I_z,
            joule_integral=self.joule_integral,
            ultim_break_capacity=ultim_break_capacity,
            let_through_energy=E_through,
            k_magn_trip=k_magn_trip,
            t_m_lim=t_m_lim
        )
        cb.check_overload_protection()
        cb.check_shortcircuit_protection(I_sc_max, I_sc_min)
        self.circuit_breaker = cb
        self.I_sc_max = I_sc_max
        self.I_sc_min = I_sc_min
        return None

    def check_indirect_contact_protection(
        self,
        S_pe: Quantity | None = None,
        skin_condition: str = "BB2",
        final_circuit: bool = True
    ) -> IndirectContactProtectionResult | None:
        """
        Returns the requirements so that the circuit breaker to which the
        cable is connected also protects against indirect contact in a
        low-voltage system.

        Parameters
        ----------
        S_pe: Quantity, optional
            Cross-sectional area of the earth protective conductor. If given, it
            overrides the value of instance attribute `self.S_pe`.
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the condition of the human skin: either dry
            (BB1) or wet (BB2).
        final_circuit: bool, default True
            Indicates that the cable directly feeds a consumer (of which the
            nominal current does not exceed 32 A).

        Returns
        -------
        IndirectContactProtectionResult | None
        """
        if self.earthing_system.is_TN():
            return earthing_system.TN.check_indirect_contact(
                U_phase=self.U_phase,
                L_cable=self.L,
                S_phase=self.S,
                conductor_material=self.conductor_material,
                I_m_cb= self.circuit_breaker.I_m_max if self.circuit_breaker else None,
                S_pe = self.S_pe if S_pe is None else S_pe,
                skin_condition = skin_condition,
                final_circuit = final_circuit
            )
        return None


class ThreePhaseCable(AbstractCable):
    phase_system = PhaseSystem.THREE_PHASE


class SinglePhaseCable(AbstractCable):
    phase_system = PhaseSystem.SINGLE_PHASE


def plot_cable_characteristic(cable: AbstractCable):
    t_min, t_max = 1e-6, 1e6
    t_5s = 5.0
    ji = cable.joule_integral.to('A ** 2 * s').m
    I_z = cable.I_z.to('A').m
    I_min = 1e-3

    t_adiabatic_ar1 = np.linspace(t_min, t_5s)
    I_adiabatic_ar1 = np.sqrt(ji / t_adiabatic_ar1)
    I_max = np.max(I_adiabatic_ar1)

    t_ampacity_ini = ji / (I_z ** 2)
    t_adiabatic_ar_2 = np.linspace(t_5s, t_ampacity_ini)
    I_adiabatic_ar_2 = np.sqrt(ji / t_adiabatic_ar_2)

    t_ampacity_ar = np.linspace(t_ampacity_ini, t_max)
    I_ampacity_ar = np.array([I_z for _ in range(len(t_ampacity_ar))])

    chart = LineChart()
    chart.x1.set_log_scale()
    chart.y1.set_log_scale()

    chart.add_xy_data(
        label="joule-integral-1",
        x1_values=I_adiabatic_ar1,
        y1_values=t_adiabatic_ar1,
        style_props={"color": "tab:green"}
    )
    chart.add_xy_data(
        label="joule-integral-2",
        x1_values=I_adiabatic_ar_2,
        y1_values=t_adiabatic_ar_2,
        style_props={"linestyle": "--", "color": "tab:green"}
    )
    chart.add_xy_data(
        label="ampacity",
        x1_values=I_ampacity_ar,
        y1_values=t_ampacity_ar,
        style_props={"color": "tab:green"}
    )

    if cable.circuit_breaker:
        I_nf = cable.circuit_breaker.I_nf.to('A').m
        I_f = cable.circuit_breaker.I_f.to('A').m
        t_conv = cable.circuit_breaker.t_conv.to('s').m
        I_m_min = cable.circuit_breaker.I_m_min.to('A').m
        I_m_max = cable.circuit_breaker.I_m_max.to('A').m
        t_m_lim = cable.circuit_breaker.t_m_lim.to('s').m
        I_cu = cable.circuit_breaker.I_cu.to('A').m
        I_sc_min = cable.I_sc_min.to('A').m
        I_sc_max = cable.I_sc_max.to('A').m

        color_overload = "tab:orange"
        color_sc = "tab:red"
        alpha_min = 0.7
        alpha_max = 1.0
        alpha_time = 1.0

        cb_lines = [
            ("t_conv", (I_min, I_max), (t_conv, t_conv), color_overload, alpha_time),
            ("I_nf", (I_nf, I_nf), (t_min, t_max), color_overload, alpha_min),
            ("I_f", (I_f, I_f), (t_min, t_max), color_overload, alpha_max),
            ("t_m_lim", (I_min, I_max), (t_m_lim, t_m_lim), color_sc, alpha_time),
            ("I_m_min", (I_m_min, I_m_min), (t_min, t_max), color_sc, alpha_min),
            ("I_m_max", (I_m_max, I_m_max), (t_min, t_max), color_sc, alpha_max),
            ("I_cu", (I_cu, I_cu), (t_min, t_max), color_sc, alpha_max)
        ]
        for line in cb_lines:
            chart.add_xy_data(
                label=line[0],
                x1_values=line[1],
                y1_values=line[2],
                style_props={
                    "linestyle": "-" if line[0] == "I_cu" else "--",
                    "color": line[3],
                    "alpha": line[4],
                    "linewidth": 3.0 if line[0] == "I_cu" else 1.5
                }
            )

        sc_lines = [
            ("I_sc_min", (I_sc_min, I_sc_min), (t_min, t_max), "tab:purple"),
            ("I_sc_max", (I_sc_max, I_sc_max), (t_min, t_max), "tab:purple")
        ]
        for line in sc_lines:
            chart.add_xy_data(
                label=line[0],
                x1_values=line[1],
                y1_values=line[2],
                style_props={"color": line[3]}
            )

    chart.x1.scale(I_min, I_max)
    chart.y1.scale(t_min, t_max)

    return chart


@dataclass
class SelectivityResult:
    total: bool
    exists: bool
    t_trip_max: Quantity | None = None
    t_margin: Quantity | None = None


def check_selectivity(
    cable_up: AbstractCable,
    cable_down: AbstractCable,
) -> SelectivityResult:
    """
    Checks the selectivity between the circuit breaker of the downstream cable
    and the circuit breaker of the nearest upstream cable.

    Parameters
    ----------
    cable_up: AbstractCable
        Upstream cable with circuit breaker connected.
    cable_down: AbstractCable
        Downstream cable with circuit breaker connected.

    Returns
    -------
    SelectivityResult
        `exists` is True if current based selectivity exists between the
        downstream and upstream circuit breaker, i.e., if the time-current
        characteristic curve of the upstream circuit breaker is completely to
        the right of the curve of the downstream circuit breaker.
        `total` is True if current based selectivity exists, and the calculated
        maximum short-circuit current in the downstream cable is smaller than
        the minimum magnetic tripping current of the upstream circuit breaker.
        `t_trip_max` and `t_margin` are set when the upstream circuit breaker is
        of the industrial type and adjustable, current based selectivity
        exists, but selectivity is not total. `t_trip_max` is the maximum
        allowable tripping time of the upstream circuit breaker taking account
        of the Joule-integral of the upstream cable, and also the maximum
        allowable fault duration to ensure protection against indirect contact
        (assuming a TN-earthing system). `t_margin` is the margin between
        `t_trip_max` and `t_m_lim`, the instantaneous magnetic tripping time
        limit of the circuit breaker.

    Notes
    -----
    - For residential circuit breakers (IEC 60898-1), the magnetic tripping
      time is not adjustable. The time value that could be computed from the
      Joule integral and the safety curve is therefore not returned and cannot
      be used as a setting; only current-based and total selectivity are
      evaluated.

    - For industrial circuit breakers (IEC 60947-2) with an adjustable delay,
      the returned time (if any) represents the *maximum allowable tripping
      delay* for the upstream breaker, limited by:
        * the thermal withstand (Joule integral) of the upstream cable, and
        * the maximum fault duration allowed for protection against indirect
          contact (default BB2 – wet skin).
      The actual setting must be chosen below or equal to this value.
    """
    if cable_up.has_circuit_breaker and cable_down.has_circuit_breaker:
        cb_up = cable_up.circuit_breaker
        cb_down = cable_down.circuit_breaker
        res = check_current_based_selectivity(cb_down, cb_up)
        if res:
            I_sc_max_down = cable_down.I_sc_max
            if I_sc_max_down < cb_up.I_m_min:
                return SelectivityResult(total=True, exists=True)
            elif not cb_up.has_adjustable_delay:
                return SelectivityResult(total=False, exists=True)
            t_cable_max = cable_up.joule_integral / I_sc_max_down ** 2
            I_sc_min_up = cable_up.I_sc_min
            t_cable_min_up = cable_up.joule_integral / I_sc_min_up ** 2
            res = cable_up.check_indirect_contact_protection(final_circuit=False)
            t_safety_max = res.t_c_max
            t_trip_max = min(t_cable_max, t_cable_min_up, t_safety_max)
            t_margin = t_trip_max - cb_up.t_m_lim
            return SelectivityResult(
                total=True,
                exists=True,
                t_trip_max=t_trip_max.to('ms'),
                t_margin=t_margin.to('ms')
            )
        else:
            return SelectivityResult(total=False, exists=False)
    else:
        raise ValueError(
            "To check selectivity circuit breakers must be "
            "connected to both upstream and downstream cables."
        )
