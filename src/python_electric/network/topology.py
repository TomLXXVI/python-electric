from typing import Iterator, Sequence
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import StrEnum
import warnings

from .. import Quantity, Q_
from ..utils.charts import LineChart
from ..protection import (
    PEConductor,
    CircuitBreaker,
    CircuitBreakerAdvisor,
    CircuitBreakerSuggestion
)
from ..protection.earthing_system import IndirectContactProtectionResult
from .config import SCCalcConfig, PEConductorConfig
from .graph import NetworkGraph, Connection, Component, TComponent
from .components import (
    Load,
    BusBar,
    Cable,
    Grid,
    Transformer,
    InductionMotor,
    Generator,
    SynchronousMotor
)
from .components.cable import (
    ConductorMaterial,
    InsulationMaterial,
    InstallMethod,
    CableMounting,
    CableArrangement,
    EarthingSystem,
    PhaseSystem,
    plot_cable,
    SelectivityResult
)

__all__ = [
    "NetworkTopology",
    "BusBarInput",
    "CableInput",
    "GridInput",
    "InductionMotorInput",
    "TransformerInput",
    "ShortCircuitResult"
]


@dataclass
class ShortCircuitResult:
    max: Quantity | None = None
    min: Quantity | None = None
    fault_type_min: str = ""
    fault_type_max: str = ""


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
    T_amb: Quantity = Q_(30, 'degC')
    install_method: InstallMethod = InstallMethod.E
    cable_mounting: CableMounting = CableMounting.PERFORATED_TRAY
    cable_arrangement: CableArrangement = CableArrangement.MULTICORE
    num_other_circuits: int = 0
    h3_fraction: float = 0.0
    sizing_based_on_I_nom: bool = False
    earthing_system: EarthingSystem | None = None
    phase_system: PhaseSystem = PhaseSystem.PH3
    z0_r_factor: float = 3.0
    z0_x_factor: float = 3.0
    n_phase: int = 1

    load: Load = field(init=False, default=None)

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self, load: Load) -> Component:
        self.load = load
        cable = Cable(
            name=self.name,
            I_b_ph=self.load.I_b,
            U_l=self.load.U_l,
            cos_phi=self.load.cos_phi,
            L=self.L,
            k_simul=self.k_simul,
            k_ext=self.k_ext,
            conductor_material=self.conductor_material,
            insulation_material=self.insulation_material,
            T_amb=self.T_amb,
            install_method=self.install_method,
            cable_mounting=self.cable_mounting,
            cable_arrangement=self.cable_arrangement,
            num_other_circuits=self.num_other_circuits,
            h3_fraction=self.h3_fraction,
            sizing_based_on_I_nom=self.sizing_based_on_I_nom,
            earthing_system=self.earthing_system,
            phase_system=self.phase_system,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor,
            n_phase=self.n_phase
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
class GeneratorInput(ComponentInput):
    name: str
    U_n: Quantity
    S_n: Quantity
    x: float
    cos_phi: float = 0.8
    R_to_X: float = 0.15
    z0_factor: float = 1.0

    def __post_init__(self):
        super().__init__(self.name)

    def create_component(self) -> Component:
        gen = Generator(
            name=self.name,
            U_n=self.U_n,
            S_n=self.S_n,
            x=self.x,
            cos_phi=self.cos_phi,
            R_to_X=self.R_to_X,
            z0_factor=self.z0_factor
        )
        return gen


@dataclass
class SynchronousMotorInput(GeneratorInput):

    def create_component(self) -> Component:
        motor = SynchronousMotor(
            name=self.name,
            U_n=self.U_n,
            S_n=self.S_n,
            x=self.x,
            cos_phi=self.cos_phi,
            R_to_X=self.R_to_X,
            z0_factor=self.z0_factor
        )
        return motor


class NetworkTopology:
    """
    Central class for designing a low-voltage network.
    """
    GROUND_ID = NetworkGraph.GROUND_ID

    class ShortCircuitCase(StrEnum):
        MAX = "MAX"
        MIN = "MIN"

    def __init__(
        self,
        name: str,
        U_n: Quantity,
        sc_config: SCCalcConfig | None = None,
        pe_config: PEConductorConfig | None = None,
        earthing_system: EarthingSystem = EarthingSystem.TN,
        cb_standard: CircuitBreaker.Standard = CircuitBreaker.Standard.INDUSTRIAL,
        skin_condition: str = "BB2",
        neutral_distributed: bool = True,
        R_e: Quantity | None = None
    ) -> None:
        """
        Initializes a `NetworkTopology` object.

        Parameters
        ----------
        name: str
            Name for the network.
        U_n: Quantity
            Nominal network voltage.
        sc_config: SCCalcConfig, optional
            Global configuration settings for the maximum/minimum short-circuit
            calculations. If None, the default configuration is used (see
            python_electric/network/short_circuit_calc/sequence_network_builder.py).
        pe_config: PEConductorConfig, optional
            Global configuration settings for sizing PE-conductors. If None, the
            default configuration is used (see python_electric/protection/earthing.py).
        cb_standard: CircuitBreaker.Standard, default CircuitBreaker.Standard.INDUSTRIAL
            IEC-standard to which the circuit breakers in the LV-network apply.
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the condition of the human skin: either dry
            (BB1) or wet (BB2).
        neutral_distributed: bool, default True
            Indicates whether the neutral conductor is also distributed in the
            network or not. Only used in an IT-earthing system.
        R_e: Quantity, optional
            Earth-spreading resistance of the consumer installation. Only
            used in an IT-earthing system.
        """
        self.name = name
        self.U_n = U_n
        self.earthing_system = earthing_system
        self.cb_standard = cb_standard
        self.skin_condition = skin_condition
        self.neutral_distributed = neutral_distributed
        self.R_e = R_e

        # Configuration settings for the min/max short-circuit calculation
        self._glob_sc_config = sc_config
        self.sccalc = None
        self._sc_dict: dict[str, ShortCircuitResult] = {}  # bus_id -> min/max short-circuit current

        # Global configuration settings applied to all PE-conductors
        if pe_config is None:
            self._glob_pe_cfg = PEConductorConfig()
        else:
            self._glob_pe_cfg = pe_config

        # Register of all the components in the network
        # component.name -> conn_id
        self._comp_register: dict[str, str] = {}

        # Loads are associated with connections -> see add_connection()
        self._loads: dict[str, Load | None] = {}

        # Core network graph object.
        self._nw_graph = NetworkGraph(name)

    @property
    def graph(self) -> NetworkGraph:
        return self._nw_graph

    def add_source_grid_connection(
        self,
        conn_id: str,
        end_id: str,
        U_l: Quantity,
        S_sc: Quantity,
        R_to_X: float = 0.1,
        z0_r_factor: float = 0.0,
        z0_x_factor: float = 0.0
    ) -> None:
        grid = GridInput(
            name="grid",
            U_l=U_l,
            S_sc=S_sc,
            R_to_X=R_to_X,
            z0_r_factor=z0_r_factor,
            z0_x_factor=z0_x_factor
        )
        self.add_connection(conn_id, self.GROUND_ID, end_id)
        self.add_component(conn_id, grid)

    def add_connection(
        self,
        conn_id: str,
        start_id: str = GROUND_ID,
        end_id: str = GROUND_ID,
        load: Load | None = None
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
        load: LoadInput, optional
            Contains the user-input data about the load that flows through this
            connection.

        Returns
        -------
        None
        """
        if start_id == end_id:
            raise ValueError(
                "A connection cannot have equal start_id and end_id."
            )
        self._nw_graph.add_connection(conn_id, start_id, end_id)
        if load is not None:
            load.U_l = (
                load.U_l if load.U_l is not None
                else self.U_n
            )
        self._loads[conn_id] = load

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
            load = self._loads[conn_id]
        except KeyError:
            raise KeyError(f"Connection '{conn_id}' not found.")

        comp = None

        if isinstance(comp_data, (BusBarInput, CableInput, TransformerInput)):
            if isinstance(comp_data, CableInput):
                if comp_data.earthing_system is None:
                    comp_data.earthing_system = self.earthing_system
            comp = comp_data.create_component(load)
        elif isinstance(comp_data, ComponentInput):
            comp = comp_data.create_component()

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
        self._nw_graph.add_component(conn_id, comp)

    def get_component(self, comp_id: str) -> TComponent:
        """
        Returns the network component with the specified ID.
        """
        try:
            conn_id = self._comp_register[comp_id]
        except KeyError:
            raise KeyError(f"Component '{comp_id}' not found.")

        return self._nw_graph.get_component(conn_id, comp_id)

    def get_cable(self, cable_id: str) -> tuple[Cable, Connection]:
        """
        Returns the Cable or BusBar object with the specified ID and the
        Connection object to which it belongs.
        """
        try:
            conn_id = self._comp_register[cable_id]
        except KeyError:
            raise KeyError(f"Cable '{cable_id}' not found.")

        cable = self._nw_graph.get_component(conn_id, cable_id)

        if not isinstance(cable, (Cable, BusBar)):
            cls = cable.__class__.__name__
            raise ValueError(f"Component '{cable_id}' is not of type Cable or BusBar, but {cls}.")

        conn = self._nw_graph.get_connection(conn_id)
        return cable, conn

    def iter_all_cables(self) -> Iterator[Cable]:
        """
        Returns an iterator over the cables in the network.
        """
        cables = self._nw_graph.find_components(Cable)
        for cable, *_ in cables:
            yield cable

    @property
    def glob_pe_config(self) -> PEConductorConfig:
        """Returns the global settings for sizing PE-conductors."""
        return self._glob_pe_cfg

    @glob_pe_config.setter
    def glob_pe_config(self, v: PEConductorConfig) -> None:
        """Sets the global settings for sizing PE-conductors."""
        self._glob_pe_cfg = v

    def get_load(self, conn_id: str) -> Load:
        """Returns the Load object associated with the specified connection."""
        try:
            load = self._loads[conn_id]
        except KeyError:
            raise KeyError(f"No load data available with connection '{conn_id}'.")
        return load

    def get_shortcircuit_current(
        self,
        bus_id: str,
        sc_case: ShortCircuitCase,
        sc_config: SCCalcConfig | None = None
    ) -> ShortCircuitResult:
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
        ShortCircuitResult
            The maximum/minimum short-circuit current.
        """
        from python_electric.network.short_circuit_calc import ShortCircuitCalc

        if len(self._nw_graph.connections) < 1:
            raise ValueError("The network has no connections yet.")

        if isinstance(sc_config, SCCalcConfig):
            self._glob_sc_config = sc_config

        if self.sccalc is None or sc_config is not None:
            self.sccalc = ShortCircuitCalc(self, self._glob_sc_config)

        if sc_case == self.ShortCircuitCase.MAX:
            try:
                return self.sccalc.max(bus_id)
            except Exception as err:
                raise ValueError(
                    f"Short-circuit calculation failed on bus '{bus_id}'."
                    f"Got exception: {type(err).__name__} - {err}"
                )
        elif sc_case == self.ShortCircuitCase.MIN:
            try:
                return self.sccalc.min(bus_id)
            except Exception as err:
                raise ValueError(
                    f"Short-circuit calculation failed on bus '{bus_id}'. "
                    f"Got exception: {type(err).__name__} - {err}"
                )
        else:
            raise ValueError(f"Unknown short-circuit case: '{sc_case}'.")

    def run_shortcircuit_calc(
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
            self._glob_sc_config = sc_config

        # Empty the short-circuit dict each time run_short_circuit_calculation()
        # is called.
        self._sc_dict: dict[str, ShortCircuitResult] = {}

        # Calculate maximum/minimum short-circuit currents at the busses.
        for bus_id in self._nw_graph.busses:
            if bus_id != NetworkGraph.GROUND_ID:
                sc_res1 = self.get_shortcircuit_current(bus_id, self.ShortCircuitCase.MAX)
                sc_res2 = self.get_shortcircuit_current(bus_id, self.ShortCircuitCase.MIN)
                self._sc_dict[bus_id] = ShortCircuitResult(
                    max=sc_res1.max,
                    min=sc_res2.min,
                    fault_type_max=sc_res1.fault_type_max,
                    fault_type_min=sc_res2.fault_type_min
                )

        # Assign the maximum/minimum short-circuit currents to cables/busbars.
        if self._sc_dict:
            cables = self._nw_graph.find_components(Cable)
            for cable, _, start_id, end_id in cables:
                cable.I_sc_max = self._sc_dict[start_id].max
                cable.I_sc_min = self._sc_dict[end_id].min

            busbars = self._nw_graph.find_components(BusBar)
            for busbar, _, start_id, end_id in busbars:
                busbar.I_sc_max = self._sc_dict[start_id].max
                busbar.I_sc_min = self._sc_dict[end_id].min

        return self._sc_dict

    def suggest_circuit_breaker(
        self,
        cable_id: str,
        *,
        safety_factor_Icu: float = 1.0,
        icu_series: Sequence[Quantity | float] | None = None,
        km_grid: Sequence[float] | None = None,
        prefer_adjustable: bool = False,
    ) -> CircuitBreakerSuggestion:
        cable, _ = self.get_cable(cable_id)

        advisor = CircuitBreakerAdvisor(
            cable=cable,
            standard=self.cb_standard,
            safety_factor_Icu=safety_factor_Icu,
            icu_series=icu_series,
            km_grid=km_grid,
            prefer_adjustable=prefer_adjustable
        )

        best = advisor.suggest(top=1)
        if not best:
            all_results = advisor.suggest_all()
            return all_results[0]
        else:
            return best[0]

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
            if self._glob_sc_config is None:
                warnings.warn(
                    "Short-circuit currents have not been calculated yet. "
                    "These are calculated now with a default short-circuit "
                    "configuration.", category=UserWarning
                )
            self.run_shortcircuit_calc()

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

    def check_selectivity(self, up_id: str, down_id: str) -> SelectivityResult:
        from python_electric.network.components import cable

        cable_up, _ = self.get_cable(up_id)
        cable_down, _ = self.get_cable(down_id)

        res = cable.check_selectivity(
            cable_up, cable_down,
            neutral_distributed=self.neutral_distributed
        )
        return res

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
        pe_conductor_cfg = self._glob_pe_cfg if pe_conductor_cfg is None else pe_conductor_cfg

        # Get minimum short-circuit current
        cable, conn = self.get_cable(cable_id)
        sc_res = self.get_shortcircuit_current(conn.end.name, self.ShortCircuitCase.MIN)

        pe_conductor = PEConductor(
            conductor_material=pe_conductor_cfg.cond_mat,
            insulation_material=pe_conductor_cfg.insul_mat,
            separated=pe_conductor_cfg.separated
        )
        S_pe = pe_conductor.cross_section_area(
            If=sc_res.min,
            t_interrupt=pe_conductor_cfg.t_interrupt,
            mech_protected=pe_conductor_cfg.mech_protected
        )
        if not pe_conductor.separated and S_pe < cable.S:
            cable.S_pe = cable.S
        else:
            cable.S_pe = S_pe
        return S_pe

    def size_all_pe_conductors(self):
        """
        Sizes the PE-conductor of all cables in the network using the global
        configuration settings in self._pe_glob_config.
        """
        if not self._sc_dict:
            if self._glob_sc_config is None:
                warnings.warn(
                    "Short-circuit currents have not been calculated yet. "
                    "These are calculated now with a default short-circuit "
                    "configuration.", category=UserWarning
                )
            self.run_shortcircuit_calc()

        cables = self._nw_graph.find_components(Cable)
        for cable, *_, end_id in cables:
            pe_conductor = PEConductor(
                self._glob_pe_cfg.cond_mat,
                self._glob_pe_cfg.insul_mat,
                self._glob_pe_cfg.separated
            )
            S_pe = pe_conductor.cross_section_area(
                cable.I_sc_min,
                self._glob_pe_cfg.t_interrupt,
                self._glob_pe_cfg.mech_protected
            )
            if pe_conductor.separated:
                if cable.S <= Q_(16, 'mm ** 2'):
                    cable.S_pe = max(S_pe, cable.S)
                elif Q_(16, 'mm ** 2') < cable.S <= Q_(35, 'mm ** 2'):
                    cable.S_pe = max(S_pe, Q_(16, 'mm ** 2'))
                elif cable.S > Q_(35, 'mm ** 2'):
                    cable.S_pe = max(S_pe, PEConductor.get_std_csa(cable.S / 2))
            else:
                cable.S_pe = max(S_pe, cable.S)

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

    def get_voltage_drop_of_conn(self, conn_id: str) -> tuple[Quantity, Quantity]:
        """
        Returns the absolute and relative voltage drop across the specified
        connection.
        """
        conn = self._nw_graph.get_connection(conn_id)
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

    def get_voltage_drop_at_bus(
        self,
        bus_id: str,
        *,
        start_bus_id: str | None = None,
    ) -> tuple[Quantity, Quantity]:
        """
        Return the cumulative voltage drop (absolute and relative) from a start
        bus to a target bus.

        If `start_bus_id` is None, the method tries to infer the source bus as
        follows:
        -   Find all connections that start from `ground`.
        -   If exactly one such connection exists, the inferred start bus is the
            end bus of that connection (i.e. the first energized bus downstream
            of ground).
        -   Otherwise, the caller must provide `start_bus_id`.

        Parameters
        ----------
        bus_id:
            Target bus ID where the cumulative voltage drop is requested.
        start_bus_id:
            Start bus ID. If None, the method attempts to infer the source bus
            of the network to be used as the start bus.

        Returns
        -------
        (dU, dU_rel):
            Absolute voltage drop in volts and relative voltage drop in percent.

        Raises
        ------
        KeyError:
            If the bus does not exist.
        ValueError:
            If no directed path exists, or if `start_bus_id` cannot be inferred.
        """
        # Validate target bus exists (Network will also check, but this gives a
        # clearer error here)
        if bus_id not in self._nw_graph.busses:
            raise KeyError(f"Bus '{bus_id}' not found.")

        # Infer start bus if not provided
        if start_bus_id is None:
            start_bus_ids = self._nw_graph.get_source_bus_ids()
            if len(start_bus_ids) != 1:
                raise ValueError(
                    "Cannot infer start_bus_id automatically: "
                    f"found {len(start_bus_ids)} connection(s) starting "
                    f"from '{self.GROUND_ID}'. "
                    f"Please provide start_bus_id explicitly."
                )
            start_bus_id = start_bus_ids[0]

        # Path as ordered list of connection IDs from start_bus_id -> bus_id
        conn_path = self._nw_graph.find_connection_path(start_bus_id, bus_id)

        # Sum voltage drops over the connections on the path
        dU = Q_(0.0, "V")
        dU_rel = Q_(0.0, "pct")
        for conn_id in conn_path:
            dU_i, dU_rel_i = self.get_voltage_drop_of_conn(conn_id)
            dU += dU_i
            dU_rel += dU_rel_i

        return dU, dU_rel

    def check_indirect_contact_protection(
        self,
        cable_id: str,
        *,
        skin_condition: str = "BB2",
        final_circuit: bool = True,
        neutral_distributed: bool = True,
        R_e: Quantity | None = None
    ) -> IndirectContactProtectionResult:
        """
        Checks whether the specified cable is protected against indirect
        contact.

        Parameters
        ----------
        cable_id: str
            Identifies the cable in the network.
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the condition of the human skin: either dry
            (BB1) or wet (BB2).
        final_circuit: bool, default True
            Indicates that the cable is feeding an electrical consumer of which
            the nominal current does not exceed 32 A.
        neutral_distributed: bool, default True
            Only used with IT-earthing. Indicates whether the neutral conductor
            is also distributed in the network or not.
        R_e: Quantity, optional
            Only used with IT-earthing. Indicates the spreading resistance of
            the earth electrodes in the consumer installation. Must be set when
            two exposed conductive parts in the fault loop are both connected to
            a different, separate earth electrode. It is assumed that all earth
            electrodes in the consumer installation exhibit the same spreading
            resistance.

        Returns
        -------
        IndirectContactProtectionResult
        """
        cable, _ = self.get_cable(cable_id)
        res = cable.check_indirect_contact_protection(
            S_pe=cable.S_pe,
            skin_condition=skin_condition,
            final_circuit=final_circuit,
            neutral_distributed=neutral_distributed,
            R_e=R_e
        )
        return res

    def check_final_circuits(
        self,
        source_bus_id: str | None = None
    ) -> dict[str, IndirectContactProtectionResult]:
        """
        Checks the protection against indirect contact of all the final
        circuits in the network.

        Parameters
        ----------
        source_bus_id: str, optional
            ID of the network source bus. If None, a source bus is searched for.

        Returns
        -------
        dict[str, IndirectContactProtectionResult]

        Raises
        ------
        ValueError
            If a network source bus cannot be found.
        """
        if source_bus_id is None:
            # Get a source bus ID (it doesn't matter which one, should there be
            # several ones)
            source_bus_ids = self._nw_graph.get_source_bus_ids()
            if not source_bus_ids:
                raise ValueError(
                    "Cannot find a network source bus."
                    "Please provide source_bus_id explicitly."
                )
            source_bus_id = source_bus_ids[0]

        # Get final circuits
        final_bus_ids = [
            final_bus.name for final_bus in
            self._nw_graph.find_reachable_final_busses(source_bus_id)
        ]
        final_conn_ids = [
            self._nw_graph.find_connection_path(source_bus_id, bus_id)[-1]
            for bus_id in final_bus_ids
        ]
        final_circuits = [
            self._nw_graph.get_connection(conn_id)
            for conn_id in final_conn_ids
        ]

        # Get cables from final circuits
        cables = [
            component
            for final_circuit in final_circuits
            for component in final_circuit.components.values()
            if isinstance(component, Cable)
        ]

        results = {
            cable.name: self.check_indirect_contact_protection(
                cable.name,
                skin_condition=self.skin_condition,
                final_circuit=True,
                neutral_distributed=self.neutral_distributed,
                R_e=self.R_e
            )
            for cable in cables
        }
        return results
