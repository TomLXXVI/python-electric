from ... import Quantity, Q_
from ...short_circuit.network import Network as SCNetwork
from ...short_circuit.network import PerUnitSystem
from ...short_circuit.faults import ThreePhaseFault, LineToGroundFault, LineToLineFault
from ..config import SCCalcConfig
from ..topology import NetworkTopology
from .sequence_network_builder import build_sequence_network, MissingConnectionError


__all__ = ["ShortCircuitCalc"]


class ShortCircuitCalc:

    class Max:

        def __init__(self, network: NetworkTopology, cfg: SCCalcConfig | None) -> None:
            self.nw = network
            if isinstance(cfg, SCCalcConfig):
                self.config.sc_case = "MAX"
            else:
                self.config = SCCalcConfig("MAX")
            self.nw1: SCNetwork | None = None
            self._fault_type = ThreePhaseFault
            self._fault = None

            self._create_sequence_network()
            self._create_fault()

        def _create_sequence_network(self) -> None:
            self.nw1 = build_sequence_network(
                self.nw.graph,
                sequence=1,
                config=self.config
            )

        def _create_fault(self) -> None:
            self._fault = self._fault_type(
                self.nw1,
                c=self.config.volt_factor_max
            )

        def __call__(self, bus_id: str) -> Quantity:
            self._fault.set_faulted_node(bus_id)
            If_pu = self._fault.get_fault_current()
            bus = self.nw.graph.busses[bus_id]
            pu_sys = PerUnitSystem(self.config.S_base, bus.U_base)
            If = pu_sys.get_actual_current(If_pu)
            return abs(If)

        def print_seq_network(self) -> None:
            for branch in self.nw1.branches:
                print(branch)

        def print_impedance_matrix(self):
            self.nw1.show_impedance_matrix()


    class Min:

        def __init__(self, network: NetworkTopology, cfg: SCCalcConfig | None) -> None:
            self.nw = network
            if isinstance(cfg, SCCalcConfig):
                self.config.sc_case = "MIN"
            else:
                self.config = SCCalcConfig("MIN")

            self.nw1: SCNetwork | None = None
            self.nw2: SCNetwork | None = None
            self.nw0: SCNetwork | None = None

            # self._fault_type = LineToLineFault if self.nw.earthing_system.is_IT() else LineToGroundFault
            self._fault_type = LineToGroundFault
            self._fault = None

            self._create_sequence_networks()
            self._create_fault()

        def _create_sequence_networks(self) -> None:
            self.nw1 = build_sequence_network(
                self.nw.graph,
                sequence=1,
                config=self.config
            )
            self.nw2 = build_sequence_network(
                self.nw.graph,
                sequence=2,
                config=self.config
            )
            try:
                self.nw0 = build_sequence_network(
                    self.nw.graph,
                    sequence=0,
                    config=self.config
                )
            except MissingConnectionError:
                self.nw0 = None

        def _create_fault(self) -> None:
            self._fault = self._fault_type(
                [self.nw0, self.nw1, self.nw2],
                c=self.config.volt_factor_min
            )

        def __call__(self, bus_id: str) -> Quantity:
            self._fault.set_faulted_node(bus_id)
            try:
                If_pu = self._fault.get_fault_current_abc()
            except TypeError:
                # Happens when a bus is on a disconnected island in the Z0-network.
                return Q_(0.0, 'A')
            bus = self.nw.graph.busses[bus_id]
            pu_sys = PerUnitSystem(self.config.S_base, bus.U_base)
            if issubclass(self._fault_type, LineToGroundFault):
                If = pu_sys.get_actual_current(If_pu.flatten()[0])
            else:
                If = pu_sys.get_actual_current(If_pu.flatten()[1])
            return abs(If)

        def print_seq_network(self, seq: int) -> None:
            if seq == 1:
                nw_seq = self.nw1
            elif seq == 2:
                nw_seq = self.nw2
            elif seq == 0:
                nw_seq = self.nw0
            else:
                raise ValueError(f"seq must be 0, 1, or 2. Got {seq} instead.")

            for branch in nw_seq.branches:
                print(branch)

        def print_impedance_matrix(self, seq):
            if seq == 0:
                nw = self.nw0
            elif seq == 1:
                nw = self.nw1
            elif seq == 2:
                nw = self.nw2
            else:
                raise ValueError(f"seq must be 0, 1, or 2. Got {seq} instead.")

            nw.show_impedance_matrix()


    def __init__(
        self,
        network: NetworkTopology,
        cfg: SCCalcConfig | None = None
    ) -> None:
        self.max = ShortCircuitCalc.Max(network, cfg)
        self.min = ShortCircuitCalc.Min(network, cfg)
