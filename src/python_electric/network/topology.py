import typing
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from enum import StrEnum
import warnings

from .. import Quantity, Q_
from ..utils.charts import LineChart
from ..protection import PEConductor, CircuitBreaker
from .config import SCCalcConfig, PEConductorConfig
from .network import Network, Connection, Component, TComponent
from .short_circuit_calc import ShortCircuitCalc
from .components import (
    Load,
    BusBar,
    Cable,
    Grid,
    Transformer,
    InductionMotor
)
from .components.cable import (
    ConductorMaterial,
    InsulationMaterial,
    Ambient,
    InstallMethod,
    CableMounting,
    CableArrangement,
    EarthingSystem,
    PhaseSystem,
    plot_cable,
    check_selectivity,
    SelectivityResult
)


__all__ = [
    "NetworkTopology",
    "BusBarInput",
    "CableInput",
    "GridInput",
    "InductionMotorInput",
    "TransformerInput",
    "LoadInput"
]


class ComponentInput(ABC):

    def __init__(self, name: str) -> None:
        self.name = name

    @abstractmethod
    def create_component(self, *args, **kwargs) -> Component:
        ...


@dataclass
class BusBarInput(ComponentInput):
    name: str
    S: Quantity
    I_z: Quantity
    L: Quantity
    conductor_material: ConductorMaterial = ConductorMaterial.COPPER
    T_max_cont: Quantity = Q_(100, 'degC')
    z0_r_factor: float = 1.0
    z0_x_factor: float = 1.0

    load: Load = field(init=False, default=None)

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self, load: Load) -> BusBar:
        self.load = load
        busbar = BusBar(
            name=self.name,
            S=self.S,
            I_z=self.I_z,
            L=self.L,
            U_l=self.load.U_l,
            I_b=self.load.I_b,
            cos_phi=self.load.cos_phi,
            conductor_material=self.conductor_material,
            T_max_cont=self.T_max_cont,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor
        )
        return busbar


@dataclass
class CableInput(ComponentInput):
    name: str
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
    z0_r_factor: float = 3.0
    z0_x_factor: float = 3.0

    load: Load = field(init=False, default=None)

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self, load: Load) -> Component:
        self.load = load
        cable = Cable(
            name=self.name,
            I_b=self.load.I_b,
            U_l=self.load.U_l,
            cos_phi=self.load.cos_phi,
            L=self.L,
            k_simul=self.k_simul,
            k_ext=self.k_ext,
            conductor_material=self.conductor_material,
            insulation_material=self.insulation_material,
            ambient=self.ambient,
            T_amb=self.T_amb,
            install_method=self.install_method,
            cable_mounting=self.cable_mounting,
            cable_arrangement=self.cable_arrangement,
            num_circuits=self.num_circuits,
            h3_fraction=self.h3_fraction,
            sizing_based_on_I_nom=self.sizing_based_on_I_nom,
            earthing_system=self.earthing_system,
            phase_system=self.phase_system,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor
        )
        return cable


@dataclass
class GridInput(ComponentInput):
    name: str
    U_l: Quantity
    S_sc: Quantity
    R_to_X: float = 0.1
    z0_r_factor: float = 1.0
    z0_x_factor: float = 3.0

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self) -> Component:
        grid = Grid(
            name=self.name,
            U_l=self.U_l,
            S_sc=self.S_sc,
            R_to_X=self.R_to_X,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor
        )
        return grid


@dataclass
class InductionMotorInput(ComponentInput):
    name: str
    U_n: Quantity
    I_n: Quantity
    P_m: Quantity
    eta: Quantity
    cos_phi: float
    k_start: float = 6.0
    R_to_X: float = 0.42
    z2_factor: float = 1.0
    z0_factor: float = 10.0

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self) -> Component:
        motor = InductionMotor(
            name=self.name,
            U_n=self.U_n,
            I_n=self.I_n,
            P_m=self.P_m,
            eta=self.eta,
            cos_phi=self.cos_phi,
            k_start=self.k_start,
            R_to_X=self.R_to_X,
            z2_factor=self.z2_factor,
            z0_factor=self.z0_factor
        )
        return motor


@dataclass
class TransformerInput(ComponentInput):
    WindingConn = Transformer.WindingConn

    name: str
    S_n: Quantity
    U_lp: Quantity
    U_ls: Quantity
    u_cc: Quantity
    P_Cu: Quantity
    pri_conn: WindingConn | None = None
    sec_conn: WindingConn | None = None
    Zn_pri: Quantity = Q_(0, 'ohm')
    Zn_sec: Quantity = Q_(0, 'ohm')
    z0_r_factor: float = 1.0
    z0_x_factor: float = 1.0

    load: Load = field(init=False, default=None)

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self, load: Load) -> Component:
        self.load = load
        transfo = Transformer(
            name=self.name,
            S_n=self.S_n,
            U_lp=self.U_lp,
            U_ls=self.U_ls,
            u_cc=self.u_cc,
            P_Cu=self.P_Cu,
            I_b=load.I_b,
            cos_phi=load.cos_phi,
            pri_conn=self.pri_conn,
            sec_conn=self.sec_conn,
            Zn_pri=self.Zn_pri,
            Zn_sec=self.Zn_sec,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor
        )
        return transfo


@dataclass
class LoadInput:
    U_l: Quantity | None = None
    cos_phi: float = 0.8
    P_e: Quantity | None = None
    P_m: Quantity | None = None
    eta: Quantity = Q_(100, 'pct')
    I_b: Quantity | None = None

    def to_dict(self) -> dict[str, Quantity | float | None]:
        return asdict(self)


class NetworkTopology:
    """
    Central class for designing a low-voltage network.
    """
    GROUND_ID = Network.GROUND_ID

    class ShortCircuitCase(StrEnum):
        MAX = "MAX"
        MIN = "MIN"

    @dataclass
    class ShortCircuitResult:
        max: Quantity
        min: Quantity

    def __init__(
        self,
        name: str,
        U_lv: Quantity,
        sc_config: SCCalcConfig | None = None,
        pe_config: PEConductorConfig | None = None
    ) -> None:
        """
        Initializes a `NetworkTopology` object.

        Parameters
        ----------
        name: str
            Name for the network.
        U_lv: Quantity
            Nominal network voltage.
        sc_config: SCCalcConfig, optional
            Global configuration settings for the maximum/minimum short-circuit
            calculations. If None, the default configuration is used (see
            python_electric/network/short_circuit_calc/sequence_network_builder.py).
        pe_config: PEConductorConfig, optional
            Global configuration settings for sizing PE-conductors. If None, the
            default configuration is used (see python_electric/protection/earthing.py).
        """
        self.name = name
        self.U_lv = U_lv

        self._network = Network(name)

        # Register of all the components in the network
        self._comp_register: dict[str, str] = {}  # component.name, conn_id

        # Configuration settings for the min/max short-circuit calculation
        self._sc_glob_config = sc_config

        # Global configuration settings applied to all PE-conductors
        if pe_config is None:
            self._pe_conductor_glob_cfg = PEConductorConfig()
        else:
            self._pe_conductor_glob_cfg = pe_config

        self._loads: dict[str, LoadInput | None] = {}  # loads are associated with connections -> see add_connection()

        self._sc_calc: ShortCircuitCalc | None = None
        self._sc_dict: dict[str, NetworkTopology.ShortCircuitResult] = {}  # bus_id, min/max short-circuit current

    def add_connection(
        self,
        conn_id: str,
        start_id: str = GROUND_ID,
        end_id: str = GROUND_ID,
        load_data: LoadInput | None = None
    ) -> None:
        """
        Adds a new connection to the low-voltage network.

        Parameters
        ----------
        conn_id: str
            Identifies the connection.
        start_id: str
            Identifies the network bus where the connection leaves.
        end_id: str
            Identifies the network bus where the connection arrives.
        load_data: LoadInput, optional
            Contains the user-input data about the load that flows through this
            connection.

        Returns
        -------
        None
        """
        if start_id == end_id:
            raise ValueError("A connection cannot have equal start_id and end_id.")
        self._network.add_connection(conn_id, start_id, end_id)
        if load_data is not None:
            load_data.U_l = load_data.U_l if load_data.U_l is not None else self.U_lv
        self._loads[conn_id] = load_data

    def add_component(
        self,
        conn_id: str,
        comp_data: ComponentInput
    ) -> None:
        """
        Adds a single component to a network connection.

        Parameters
        ----------
        conn_id: str
            Identifies the connection.
        comp_data: ComponentInput
            Contains the user-input data about the component.

        Returns
        -------
        None
        """
        try:
            load_data = self._loads[conn_id]
        except KeyError:
            raise KeyError(f"Connection '{conn_id}' not found.")

        comp = None
        if isinstance(comp_data, (BusBarInput, CableInput)):
            comp = comp_data.create_component(Load(**load_data.to_dict()))
        if isinstance(comp_data, GridInput):
            comp = comp_data.create_component()
        if isinstance(comp_data, InductionMotorInput):
            comp = comp_data.create_component()
        if isinstance(comp_data, TransformerInput):
            comp = comp_data.create_component(Load(**load_data.to_dict()))

        if comp is None and comp_data.name:
            raise ValueError(f"Failed to add component '{comp_data.name}' to network.")
        if comp is None:
            raise ValueError("Failed to add component.")

        conn_id_ret = self._comp_register.setdefault(comp.name, conn_id)
        if conn_id_ret != conn_id:
            raise KeyError(
                f"A component with ID '{comp.name}' has already been "
                f"added to the network."
            )
        self._network.add_component(conn_id, comp)

    def get_component(self, comp_id: str) -> TComponent:
        """
        Returns the network component with the specified ID.
        """
        try:
            conn_id = self._comp_register[comp_id]
        except KeyError:
            raise KeyError(f"Component '{comp_id}' not found.")

        return self._network.get_component(conn_id, comp_id)

    def get_cable(self, cable_id: str) -> tuple[Cable, Connection]:
        """
        Returns the Cable or BusBar object with the specified ID and the
        Connection object to which it belongs.
        """
        try:
            conn_id = self._comp_register[cable_id]
        except KeyError:
            raise KeyError(f"Cable '{cable_id}' not found.")

        cable = self._network.get_component(conn_id, cable_id)

        if not isinstance(cable, (Cable, BusBar)):
            cls = cable.__class__.__name__
            raise ValueError(f"Component '{cable_id}' is not of type Cable or BusBar, but {cls}.")

        conn = self._network.get_connection(conn_id)
        return cable, conn

    def iter_cables(self) -> typing.Iterator[Cable]:
        """
        Returns an iterator over the cables in the network.
        """
        cables = self._network.find_components(Cable)
        for cable, *_ in cables:
            yield cable

    @property
    def global_config_pe_conductors(self) -> PEConductorConfig:
        """Returns the global settings for sizing PE-conductors."""
        return self._pe_conductor_glob_cfg

    def get_load(self, conn_id: str) -> Load:
        """Returns the Load object associated with the specified connection."""
        try:
            load_data = self._loads[conn_id]
        except KeyError:
            raise KeyError(f"No load data available with connection '{conn_id}'.")

        load = Load(**load_data.to_dict())
        return load

    def get_short_circuit_current(
        self,
        bus_id: str,
        sc_case: ShortCircuitCase,
        sc_config: SCCalcConfig | None = None
    ) -> Quantity:
        """
        Returns the short-circuit current at the specified bus of the network.

        Parameters
        ----------
        bus_id: str
            Specifies the bus of the network.
        sc_case: ShortCircuitCase
            Indicates whether the maximum or minimum short-circuit current
            is to be returned.
        sc_config: SCCalcConfig, optional
            Configuration settings to be used in the short-circuit calculation.
            Overwrites the global short-circuit configuration of the network.

        Returns
        -------
        Quantity
            The maximum/minimum short-circuit current.
        """
        if len(self._network.connections) < 1:
            raise ValueError("The network has no connections yet.")

        if isinstance(sc_config, SCCalcConfig):
            self._sc_glob_config = sc_config

        if self._sc_calc is None or sc_config is not None:
            self._sc_calc = ShortCircuitCalc(self._network, self._sc_glob_config)

        if sc_case == self.ShortCircuitCase.MAX:
            try:
                return self._sc_calc.max(bus_id)
            except:
                raise ValueError(
                    f"Short-circuit calculation failed on bus '{bus_id}'."
                )
        elif sc_case == self.ShortCircuitCase.MIN:
            try:
                return self._sc_calc.min(bus_id)
            except:
                raise ValueError(
                    f"Short-circuit calculation failed on bus '{bus_id}'."
                )
        else:
            raise ValueError(f"Unknown short-circuit case: '{sc_case}'.")

    def run_short_circuit_calculation(
        self,
        sc_config: SCCalcConfig | None = None
    ) -> dict[str, ShortCircuitResult]:
        """
        Calculates the maximum and minimum short-circuit current at each bus of
        the network.

        Parameters
        ----------
        sc_config: SCCalcConfig, optional
            Configuration settings to be used in the short-circuit calculations.
            Overrides the global short-circuit configuration of the network.

        Returns
        -------
        dict[str, ShortCircuitResult]
            A dictionary of which the keys are the bus IDs that map to the
            short-circuit calculation results (see dataclass ShortCircuitResult).
        """
        # If a short-circuit config is passed, replace the previous one.
        if isinstance(sc_config, SCCalcConfig):
            self._sc_glob_config = sc_config

        # Empty the short-circuit dict each time run_short_circuit_calculation()
        # is called.
        self._sc_dict = {}

        # Calculate maximum/minimum short-circuit currents at the busses.
        for bus_id in self._network.busses:
            if bus_id != Network.GROUND_ID:
                I_sc_max = self.get_short_circuit_current(bus_id, self.ShortCircuitCase.MAX)
                I_sc_min = self.get_short_circuit_current(bus_id, self.ShortCircuitCase.MIN)
                self._sc_dict[bus_id] = self.ShortCircuitResult(I_sc_max, I_sc_min)

        # Assign the maximum/minimum short-circuit currents to cables/busbars.
        if self._sc_dict:
            cables = self._network.find_components(Cable)
            for cable, _, start_id, end_id in cables:
                cable.I_sc_max = self._sc_dict[start_id].max
                cable.I_sc_min = self._sc_dict[end_id].min

            busbars = self._network.find_components(BusBar)
            for busbar, _, start_id, end_id in busbars:
                busbar.I_sc_max = self._sc_dict[start_id].max
                busbar.I_sc_min = self._sc_dict[end_id].min

        return self._sc_dict

    def size_pe_conductor(
        self,
        cable_id: str,
        pe_conductor_cfg: PEConductorConfig | None = None
    ) -> Quantity:
        """
        Sizes the single PE-conductor that belongs to the specified cable.

        Parameters
        ----------
        cable_id: str
            Identifies the cable in the network.
        pe_conductor_cfg: PEConductorConfig, optional
            Configuration settings to size the PE-conductor. Overrides the
            global configuration settings assigned to the network.

        Returns
        -------
        Quantity
            Calculated cross-sectional area of the PE-conductor.

        Notes
        -----
        Internally the cross-sectional area of the PE-conductor is also set on
        the Cable instance. If the PE-conductor is also the neutral conductor
        (TN-C), it may be that the calculated cross-sectional is overwritten due
        to other requirements that apply to PEN-conductors.
        """
        pe_conductor_cfg = self._pe_conductor_glob_cfg if pe_conductor_cfg is None else pe_conductor_cfg

        # Get minimum short-circuit current
        cable, conn = self.get_cable(cable_id)
        I_sc_min = self.get_short_circuit_current(conn.end.name, self.ShortCircuitCase.MIN)

        pe_conductor = PEConductor(
            pe_conductor_cfg.cond_mat,
            pe_conductor_cfg.insul_mat,
            pe_conductor_cfg.separated
        )
        S_pe = pe_conductor.cross_section_area(
            I_sc_min,
            pe_conductor_cfg.t_u,
            pe_conductor_cfg.mech_protected
        )
        cable.S_pe = S_pe
        return S_pe

    def connect_circuit_breaker(
        self,
        cable_id: str,
        *,
        standard: CircuitBreaker.Standard,
        category: CircuitBreaker.Category,
        I_cu: Quantity,
        k_m: float | None = None,
        E_t: Quantity | None = None,
        t_m: Quantity | None = None
    ) -> CircuitBreaker:
        """
        Connects a circuit breaker to the specified cable.

        Parameters
        ----------
        cable_id: str
            Identifies the cable in the network.
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
        I_cu: Quantity
            Ultimate breaking capacity of the circuit breaker.
        k_m: float, optional
            Multiplication factor that determines the rated magnetic trip
            current as a multiple of the thermal current setting if the circuit
            breaker is of the industrial type and adjustable.
        E_t: Quantity, optional
            Thermal energy (I²t) let through by the circuit breaker at the
            calculated maximum short-circuit current during the interruption
            time.
        t_m: Quantity, optional
            Upper limit of instantaneous magnetic tripping time with regard to
            short-circuits. By default, this time limit is set to 100 ms
            according to the residential standard IEC 60898-1.

        Returns
        -------
        CircuitBreaker
        """
        if not self._sc_dict:
            if self._sc_glob_config is None:
                warnings.warn(
                    "Short-circuit currents have not been calculated yet. "
                    "These are calculated now with a default short-circuit "
                    "configuration.", category=UserWarning
                )
            self.run_short_circuit_calculation()

        cable, _ = self.get_cable(cable_id)

        cable.connect_circuit_breaker(
            standard=standard,
            category=category,
            I_cu=I_cu,
            k_m=k_m,
            E_t=E_t,
            t_m=t_m
        )
        return cable.circuit_breaker

    def plot_cable(self, cable_id: str) -> LineChart:
        cable, _ = self.get_cable(cable_id)
        plt = plot_cable(cable)
        return plt

    def check_selectivity(self, up: str, down: str) -> SelectivityResult:
        cable_up, _ = self.get_cable(up)
        cable_down, _ = self.get_cable(down)
        res = check_selectivity(cable_up, cable_down)
        return res

    def size_all_pe_conductors(self):
        """
        Sizes the PE-conductor of all cables in the network using the global
        configuration settings in self._pe_glob_config.
        """
        if not self._sc_dict:
            if self._sc_glob_config is None:
                warnings.warn(
                    "Short-circuit currents have not been calculated yet. "
                    "These are calculated now with a default short-circuit "
                    "configuration.", category=UserWarning
                )
            self.run_short_circuit_calculation()

        cables = self._network.find_components(Cable)

        for cable, *_, end_id in cables:
            pe_conductor = PEConductor(
                self._pe_conductor_glob_cfg.cond_mat,
                self._pe_conductor_glob_cfg.insul_mat,
                self._pe_conductor_glob_cfg.separated
            )
            S_pe = pe_conductor.cross_section_area(
                cable.I_sc_min,
                self._pe_conductor_glob_cfg.t_u,
                self._pe_conductor_glob_cfg.mech_protected
            )
            cable.S_pe = S_pe

    # def add_circuit_breakers(self):
    #     """
    #     Adds a circuit-breaker to all cables in the network at once.
    #     """
    #     if not self._sc_dict:
    #         if self._sc_glob_config is None:
    #             warnings.warn(
    #                 "Short-circuit currents have not been calculated yet. "
    #                 "These are calculated now with a default short-circuit "
    #                 "configuration.", category=UserWarning
    #             )
    #         self.run_short_circuit_calculation()
    #
    #     cables = self._network.find_components(Cable)
    #
    #     for cable, conn_id, *_ in cables:
    #         cable = typing.cast(Cable, cable)
    #         cable.connect_circuit_breaker()

    def get_voltage_drop(self, conn_id: str) -> tuple[Quantity, Quantity]:
        conn = self._network.get_connection(conn_id)
        dU = Q_(0.0, 'V')
        dU_rel = Q_(0.0, 'pct')
        for comp in conn.components.values():
            try:
                comp.dU
            except AttributeError:
                continue
            dU += comp.dU
            dU_rel += comp.dU_rel
        return dU, dU_rel
