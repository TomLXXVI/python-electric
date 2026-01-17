from typing import Any
from ... import Quantity
from ...short_circuit.network import Network as SCNetwork
from ...short_circuit.network import PerUnitSystem
from ...short_circuit.faults import ThreePhaseFault, LineToGroundFault
from ..network import Network
from .sequence_network_builder import CalcConfig, build_sequence_network

__all__ = ["ShortCircuitCalc"]


class ShortCircuitCalc:

    class Max:

        def __init__(self, network: Network, cfg: dict[str, Any] | None) -> None:
            self.nw = network
            if cfg:
                self.config = CalcConfig("MAX", **cfg)
            else:
                self.config = CalcConfig("MAX")
            self.nw1: SCNetwork | None = None
            self._fault_type = ThreePhaseFault
            self._fault = None
            self._create_sequence_networks()
            self._create_fault()

        def _create_sequence_networks(self) -> None:
            self.nw1 = build_sequence_network(
                self.nw,
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
            bus = self.nw.busses[bus_id]
            pu_sys = PerUnitSystem(self.config.S_base, bus.U_base)
            If = pu_sys.get_actual_current(If_pu)
            return If

    class Min:

        def __init__(self, network: Network, cfg: dict[str, Any] | None) -> None:
            self.nw = network
            if cfg:
                self.config = CalcConfig("MIN", **cfg)
            else:
                self.config = CalcConfig("MIN")
            self.nw1: SCNetwork | None = None
            self.nw2: SCNetwork | None = None
            self.nw0: SCNetwork | None = None
            self._fault_type = LineToGroundFault
            self._fault = None

            self._create_sequence_networks()
            self._create_fault()

        def _create_sequence_networks(self) -> None:
            self.nw1 = build_sequence_network(
                self.nw,
                sequence=1,
                config=self.config
            )
            self.nw2 = build_sequence_network(
                self.nw,
                sequence=2,
                config=self.config
            )
            self.nw0 = build_sequence_network(
                self.nw,
                sequence=0,
                config=self.config
            )

        def _create_fault(self) -> None:
            self._fault = self._fault_type(
                [self.nw0, self.nw1, self.nw2],
                c=self.config.volt_factor_min
            )

        def __call__(self, bus_id: str) -> Quantity:
            self._fault.set_faulted_node(bus_id)
            If_pu = self._fault.get_fault_current_abc()
            bus = self.nw.busses[bus_id]
            pu_sys = PerUnitSystem(self.config.S_base, bus.U_base)
            If = pu_sys.get_actual_current(If_pu)
            # noinspection PyUnresolvedReferences
            return If[0]

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


    def __init__(
        self,
        network: Network,
        max_cfg: dict[str, Any] | None = None,
        min_cfg: dict[str, Any] | None = None
    ) -> None:
        self.max = ShortCircuitCalc.Max(network, max_cfg)
        self.min = ShortCircuitCalc.Min(network, min_cfg)
