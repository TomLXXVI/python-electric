"""
Cable: sizing - impedance - voltage drop - indirect contact protection

References
----------
Knockaert, J., Sergeant, P., Corne, B., Debruyne, C., Dereyne, S.,
Descheemaeker, J., Hespel, L., Vansteenberge, C., & Verhelst, B. (2015).
Laagspanningsinstallaties: technologie en ontwerp.
"""
import math
from dataclasses import dataclass, field
from copy import deepcopy

import numpy as np

from ... import Quantity, Q_, VoltReference, PhaseSystem
from ...materials import ConductorMaterial, InsulationMaterial, INSULATION_PROPS
from ...sizing.cable import *
from ...protection.earthing_system import EarthingSystem, ICPResult
from ...protection import CircuitBreaker, earthing_system, check_current_based_selectivity
from ...calc import voltage_drop
from ...utils.charts import LineChart

from ..network import Component

__all__ = [
    "Cable",
    "ConductorMaterial",
    "InsulationMaterial",
    "Ambient",
    "InstallMethod",
    "CableMounting",
    "CableArrangement",
    "plot_cable_characteristic",
    "check_selectivity"
]


@dataclass
class Cable(Component):
    """
    Represents an electrical cable in a low-voltage network.

    The cable is sized on instantiation, unless the conductor cross-sectional
    area S is already specified by the user. In that case, it is checked if the
    cable conductors are able to carry the specified load current. Should this
    be not the case, a `ValueError` exception is raised.

    Parameters
    ----------
    name: str
        Uniquely identifies the cable in the network.
    I_b: Quantity
        Rated load current through the cable.
    U_l: Quantity = Q_(400, 'V')
        Line-to-line voltage in case of a three-phase system or neutral-to-line
        voltage in case of a single-phase system.
    cos_phi: float = 0.8
        Rated power factor, i.e. cosine of the phase shift between voltage and
        rated load current.
    L: Quantity = Q_(1, 'm')
        Cable length.
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
    conductor_material: ConductorMaterial = ConductorMaterial.COPPER
        Material the cable conductors are made of.
    insulation_material: InsulationMaterial = InsulationMaterial.XLPE
        Insulation material around the conductors of the cable.
    ambient: Ambient = Ambient.AIR
        Ambient where the cable is installed, either in air or buried in the
        ground.
    T_amb: Quantity = Q_(30, 'degC')
        Ambient temperature of the cable.
    install_method: InstallMethod = InstallMethod.E
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
    h3_fraction: float = 0.0
        Fraction of third-order harmonics present in the total load current;
        a number between 0.0 and 1.0.
    sizing_based_on_I_nom: bool = False
        Indicates that the determination of the cross-sectional area of the
        conductors should be based on the standardized nominal current which is
        just greater than the rated load current. If False, conductors are sized
        based on rated load current.
    S: Quantity, optional
        Cross-sectional area of the loaded cable conductors. By default, this
        cross-sectional area is set to None, meaning that a valid
        cross-sectional area is to be determined on instantiation of the cable
        object.
    S_pe: Quantity, optional
        Cross-sectional area of the protective earth conductor. By default, this
        cross-sectional area is set equal to the cross-sectional area of the
        loaded cable conductors.
    earthing_system: EarthingScheme, default EarthingScheme.TN
        Indicates the earthing system of the low-voltage system. By default,
        this is set to the TN-earthing system.
    phase_system: PhaseSystem, default PhaseSystem.3PH
        Indicates whether the cable is three-phase or single phase
        (PhaseSystem.1PH).

    Attributes
    ----------
    I_z: Quantity
        Current-carrying capacity (ampacity) of the conductors at actual
        installation conditions.
    I_z0: Quantity
        Current-carrying capacity (ampacity) of the conductors at standard
        installation conditions.
    I_n: Quantity
        Standardized nominal current which is just greater than the rated load
        current of the cable. This current can be used to select the
        current-protective device for the cable.
    I_bc: Quantity
        Load current corrected by multiplication with the simultaneity factor
        and expansion factor. This current is used to size the cable.
    Z_dict: dict[str, dict[int, Quantity]]
        Dictionary containing the zero-sequence (key 0), positive-sequence
        (key 1), and negative-sequence (key 2) impedances of the cable
        conductors at three different temperatures:
        -   at 20°C (key "T20") used for calculating maximum short-circuit
            currents.
        -   at the maximum continuous temperature of the insulation
            (key "T_n").
        -   at 150°C (key "T150") used for calculating minimum short-circuit
            currents.
        E.g. to get Z1 at 150 °C, use Z["T150"][1].
    dU: Quantity
        Absolute voltage drop across the cable. In case of a three-phase cable,
        this voltage drop is specified with respect to ground/neutral.
    dU_rel: Quantity
        Relative voltage drop.
    I2t: Quantity
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
    U_ph: Quantity
        In case of a three-phase system, the ground-to-line voltage. Otherwise,
        equal to attribute `U_l`.

    Methods
    -------
    get_cable_impedance:
        Calculates the impedance of the cable at a given temperature and returns
        a dict with three values: 1: the positive sequence, 2: the negative
        sequence, and 0: the zero-sequence impedance.
    get_voltage_drop:
        Calculates the voltage drop across the cable. Returns a tuple with the
        absolute and relative voltage drop.
    connect_circuit_breaker:
        Connects a circuit breaker to the cable.
    check_indirect_contact_protection:
        Can be used for low-voltage systems having a TN- or IT-earthing scheme.
        Returns the maximum allowable interruption time in order to protect
        against indirect contact in case of an insulation fault. If a circuit
        breaker is connected to the cable, also returns the maximum length of
        the cable that can be protected by this circuit breaker.
    """
    name: str
    I_b: Quantity
    U_l: Quantity = Q_(400, 'V')
    cos_phi: float = 0.8
    L: Quantity = Q_(1.0, 'm')
    k_simul: float = 1.0
    k_ext: float = 1.0
    conductor_material: ConductorMaterial = ConductorMaterial.COPPER
    insulation_material: InsulationMaterial = InsulationMaterial.XLPE
    ambient: Ambient = Ambient.AIR
    T_amb: Quantity = Q_(30, 'degC')
    install_method: InstallMethod = InstallMethod.E
    cable_mounting: CableMounting = CableMounting.PERFORATED_TRAY
    cable_arrangement: CableArrangement = CableArrangement.MULTICORE
    num_circuits: int = 1
    h3_fraction: float = 0.0
    sizing_based_on_I_nom: bool = False
    earthing_system: EarthingSystem = EarthingSystem.TN
    phase_system: PhaseSystem = PhaseSystem.PH3

    S: Quantity | None = None
    S_pe: Quantity | None = None

    # Calculated attributes:
    dU: Quantity = field(init=False)
    dU_rel: Quantity = field(init=False)

    circuit_breaker: CircuitBreaker = field(init=False, default=None)

    I_bc: Quantity = field(init=False, default=None)
    I_z: Quantity = field(init=False, default=None)
    I_z0: Quantity = field(init=False, default=None)
    I_n: Quantity = field(init=False, default=None)
    I2t: Quantity = field(init=False, default=None)
    I_sc_max: Quantity = field(init=False, default=None)
    I_sc_min: Quantity = field(init=False, default=None)

    U_ph: Quantity = field(init=False, default=None)

    Z_dict: dict[str, dict[int, Quantity[float | complex]]] = field(init=False, default_factory=dict)

    def __post_init__(self):
        super().__init__(self.name)
        self.I_bc = self._get_corr_load_current()
        self.U_ph = self._get_phase_voltage()
        self._size_cable()
        self.Z_dict = self._calc_impedance_dict()
        self.dU, self.dU_rel = self.get_voltage_drop(self.U_l, self.I_b, self.cos_phi)

    def _get_corr_load_current(self) -> Quantity:
        return self.I_b * self.k_simul * self.k_ext

    def _get_phase_voltage(self) -> Quantity:
        if self.phase_system.is_ph3:
            return self.U_l / math.sqrt(3)
        else:
            return self.U_l

    def _size_cable(self):
        cd = self.__collect_data()
        if self.S is None:
            # If conductor csa is not specified, calculate minimal required csa
            self.__size_cable(cd)
        else:
            # User specified the csa of cable conductors: check if ok.
            S_user = deepcopy(self.S)
            self.__size_cable(cd)
            if self.S.to('mm ** 2') > S_user.to('mm ** 2'):
                raise ValueError(
                    f"The specified csa {S_user.to('mm ** 2'):~P.2f} is too "
                    f"small. A csa of at least {self.S.to('mm ** 2'):~P.2f} "
                    f"is required."
                ) from None
            else:
                self.__update_cable_params(cd, S_user)

    def __collect_data(self) -> CableData:
        if self.phase_system.is_ph3:
            num_loaded_conductors: int = 3
        else:
            num_loaded_conductors: int = 2

        cd = CableData(
            I_b=self.I_bc.to('A').m,
            L=self.L.to('m').m,
            conductor_type=self.conductor_material,
            insulation_type=self.insulation_material,
            ambient=self.ambient,
            T_amb=self.T_amb.to('degC').m,
            install_method=self.install_method,
            cable_mounting=self.cable_mounting,
            cable_arrangement=self.cable_arrangement,
            num_circuits=self.num_circuits,
            num_loaded_conductors=num_loaded_conductors,
            h3_content=self.h3_fraction,
        )
        return cd

    def __size_cable(self, cd: CableData) -> None:
        cd = get_size(cd, self.sizing_based_on_I_nom)
        self.S = Q_(cd.S, 'mm**2')
        self.I_z = Q_(cd.I_z, 'A')
        self.I_z0 = Q_(cd.I_z0, 'A')
        self.I_n = Q_(cd.I_n, 'A')
        self.I2t = Q_(cd.I2t, 'A**2*s')
        return None

    def __update_cable_params(self, cd: CableData, S_user: Quantity) -> None:
        cd = set_size(cd, S_user.to('mm**2').m)
        self.S = Q_(cd.S, 'mm**2')
        self.I_z = Q_(cd.I_z, 'A')
        self.I_z0 = Q_(cd.I_z0, 'A')
        self.I_n = Q_(cd.I_n, 'A')
        self.I2t = Q_(cd.I2t, 'A**2*s')
        return None

    def get_impedance(
        self,
        T: Quantity,
        z0_r_factor: float = 3.0,
        z0_x_factor: float = 3.0
    ) -> dict[int, Quantity[complex | float]]:
        """
        Returns the positive, negative, and zero-sequence impedance of the cable
        at the given temperature `T` of the conductors.

        Parameters
        ----------
        T: Quantity
            Temperature of the conductors.
        z0_r_factor: float, default 3.0
            Scaling factor applied to the positive-sequence resistance of the
            cable to determine its zero-sequence resistance.
        z0_x_factor: float, default 3.0
            Scaling factor applied to the negative-sequence reactance of the
            cable to determine its zero-sequence reactance.

        Returns
        -------
        dict[int, Quantity[complex | float]]
            A dictionary of which the keys 1, 2, 0 map to the positive,
            negative, and zero-sequence impedance of the cable at the given
            temperature `T`.
        """
        from python_electric.equipment import Cable
        c = Cable(
            length=self.L,
            cross_section_area=self.S,
            conductor_material=self.conductor_material,
            cable_arrangement=self.cable_arrangement,
            temperature=T,
            z0_r_factor=z0_r_factor,
            z0_x_factor=z0_x_factor,
            name=self.name
        )
        return {1: c.Z1, 2: c.Z2, 0: c.Z0}

    def _calc_impedance_dict(self) -> dict[str, dict[int, Quantity]]:
        insprops = INSULATION_PROPS[self.insulation_material]
        T_nom = insprops.T_max_cont
        T_rng = {"T20": Q_(20, 'degC'), "T_n": Q_(T_nom, 'degC'), "T150": Q_(150, 'degC')}
        return {k: self.get_impedance(v) for k, v in T_rng.items()}

    def get_voltage_drop(
        self,
        U_l: Quantity,
        I_b: Quantity,
        cos_phi: float,
        volt_ref: VoltReference | None = None
    ) -> tuple[Quantity, Quantity]:
        """
        Returns the voltage drop across the cable at the given load conditions.

        Parameters
        ----------
        U_l: Quantity
            Line-to-line voltage in case of a three-phase cable. Otherwise,
            in case of a single-phase cable, the neutral/ground-to-line voltage
            (phase voltage).
        I_b: Quantity
            Load current.
        cos_phi: float
            Power factor of the load.
        volt_ref: VoltReference | None = None
            In case of a three-phase cable, it indicates the reference of the
            voltage measurement: either the voltage drop is measured between
            ground (neutral) and a line, or it is measured between two lines.
            If `None`, voltage drop is taken between ground/neutral and a line.
        """
        if volt_ref is None:
            if self.phase_system.is_ph3:
                volt_ref = VoltReference.PH3_GROUND_TO_LINE
            if self.phase_system.is_ph1:
                volt_ref = VoltReference.PH1

        R = Q_(self.Z_dict["T_n"][1].m.real, 'ohm')
        X = Q_(self.Z_dict["T_n"][1].m.imag, 'ohm')
        dU = voltage_drop(R, X, I_b, volt_ref, cos_phi)

        U = U_l.to('V')
        if volt_ref == VoltReference.PH3_GROUND_TO_LINE:
            U = U / math.sqrt(3)
        dU_rel = dU / U
        return dU, dU_rel.to('pct')

    def connect_circuit_breaker(
        self,
        standard: CircuitBreaker.Standard,
        category: CircuitBreaker.Category,
        I_cu: Quantity,
        I_sc_max: Quantity,
        I_sc_min: Quantity,
        k_m: float | None = None,
        E_t: Quantity | None = None,
        t_m: Quantity | None = None
    ) -> None:
        cb = CircuitBreaker(
            standard=standard,
            category=category,
            I_b=self.I_b,
            I_n=self.I_n,
            I_z=self.I_z,
            I2t=self.I2t,
            I_cu=I_cu,
            E_t=E_t,
            k_m=k_m,
            t_m=t_m
        )
        c1 = cb.check_overload_protection()
        c2 = cb.check_shortcircuit_protection(I_sc_max, I_sc_min)
        if not (c1 and c2):
            raise ValueError("Circuit-breaker protection is inadequate.")
        self.circuit_breaker = cb
        self.I_sc_max = I_sc_max
        self.I_sc_min = I_sc_min
        return None

    @property
    def has_circuit_breaker(self) -> bool:
        return isinstance(self.circuit_breaker, CircuitBreaker)

    def check_indirect_contact_protection(
        self,
        S_pe: Quantity | None = None,
        skin_condition: str = "BB2",
        final_circuit: bool = True,
        neutral_distributed: bool = True,
        R_e: Quantity | None = None
    ) -> ICPResult | None:
        if self.earthing_system.is_TN():
            return earthing_system.TN.check_indirect_contact(
                U_phase=self.U_ph,
                L_cable=self.L,
                S_phase=self.S,
                conductor_material=self.conductor_material,
                I_m_cb= self.circuit_breaker.I_m_max if self.circuit_breaker else None,
                S_pe = self.S_pe if S_pe is None else S_pe,
                skin_condition = skin_condition,
                final_circuit = final_circuit
            )
        elif self.earthing_system.is_IT():
            return earthing_system.IT.check_indirect_contact(
                U_phase=self.U_ph,
                L_cable=self.L,
                S_phase=self.S,
                conductor_material=self.conductor_material,
                I_m_cb=self.circuit_breaker.I_m_max if self.circuit_breaker else None,
                S_pe=self.S_pe if S_pe is None else S_pe,
                skin_condition=skin_condition,
                final_circuit=final_circuit,
                neutral_distributed=neutral_distributed,
                R_e=R_e
            )
        return None


def plot_cable_characteristic(cable: Cable):
    t_min, t_max = 1e-6, 1e6
    t_5s = 5.0
    I2t = cable.I2t.to('A ** 2 * s').m
    I_z = cable.I_z.to('A').m
    I_min = 1e-3

    t_adiabatic_ar1 = np.linspace(t_min, t_5s)
    I_adiabatic_ar1 = np.sqrt(I2t / t_adiabatic_ar1)
    I_max = np.max(I_adiabatic_ar1)

    t_ampacity_ini = I2t / (I_z ** 2)
    t_adiabatic_ar_2 = np.linspace(t_5s, t_ampacity_ini)
    I_adiabatic_ar_2 = np.sqrt(I2t / t_adiabatic_ar_2)

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
        t_m_lim = cable.circuit_breaker.t_m.to('s').m
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
    cable_up: Cable,
    cable_down: Cable,
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
        exists: bool
            Is True if current based selectivity exists between the downstream
            and upstream circuit breaker, i.e., if the time-current
            characteristic curve of the upstream circuit breaker is completely
            to the right of the curve of the downstream circuit breaker.
        total: bool
            Is True if current based selectivity exists, and the calculated
            maximum short-circuit current in the downstream cable is smaller
            than the minimum magnetic tripping current of the upstream circuit
            breaker.
        t_trip_max: Quantity, optional
        t_margin: Quantity, optional
            Are set when the upstream circuit breaker is of the industrial type
            and adjustable, current based selectivity exists, but selectivity is
            not total.
            `t_trip_max` is the maximum allowable tripping time of the upstream
            circuit breaker taking account of the Joule-integral of the upstream
            cable, and also the maximum allowable fault duration to ensure
            protection against indirect contact (assuming a TN-earthing system).
            `t_margin` is the margin between `t_trip_max` and `t_m_lim`, the
            current instantaneous magnetic tripping time limit of the circuit
            breaker.

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
            t_cable_max = cable_up.I2t / I_sc_max_down ** 2
            I_sc_min_up = cable_up.I_sc_min
            t_cable_min_up = cable_up.I2t / I_sc_min_up ** 2
            res = cable_up.check_indirect_contact_protection(final_circuit=False)
            t_safety_max = res.t_c_max
            t_trip_max = min(t_cable_max, t_cable_min_up, t_safety_max)
            t_margin = t_trip_max - cb_up.t_m
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
