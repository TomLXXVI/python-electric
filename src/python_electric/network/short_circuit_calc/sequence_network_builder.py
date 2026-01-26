"""
short_circuit_builder_sequence.py

Build sequence short-circuit networks (Zbus-based) from the python-electric
`Network` graph (busses + connections + components).

This module generalizes the existing positive-sequence builder into a single
engine:

-   sequence=1 (positive): sources are included (grid always; machines optional).
-   sequence=2 (negative): sources are typically shorted/omitted by default.
-   sequence=0 (zero): requires topology rules (transformer winding connection
    and neutral grounding).

Key design choices
------------------
1)  Components return sequence impedances in Ohm (Quantity) via
    get_impedance(...). Components do NOT know any per-unit system.

2)  The builder determines the per-unit system and converts Ohm -> pu per
    component, potentially using different base voltages within a single
    Connection (e.g. when a transformer is present and impedances are referred
    to one side).

3)  A single BFS-based "network build engine" ensures branches are inserted in a
    valid order for short_circuit.network.Network.add_branch(...), which
    typically requires that at least one node is already known (except reference
    node None).

Zero-sequence (sequence=0) MVP
------------------------------
Zero-sequence modeling is not only about impedances; it is also about *which
branches exist*. In particular, transformer winding connection and neutral
grounding determine:

-   Whether Z0 current can transfer from one side to the other (delta blocks
    transfer).
-   Whether a shunt path from a bus to ground exists via a grounded neutral
    (YN + Zn).

This MVP implementation:
-   Supports at most one Transformer per Connection for Z0.
-   Adds all bus->ground shunts first (to seed the BFS).
-   Adds bus->bus series branches afterward where Z0 transfer exists.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Literal, cast

import logging

from ... import Quantity
from ...short_circuit.network.per_unit import PerUnitSystem
from ...short_circuit.network import Network as SCNetwork
from ..graph import NetworkGraph, Connection, Component
from ..config import SCCalcConfig
from ..components import Transformer

__all__ = [
    "build_sequence_network",
    "MissingConnectionError"
]


logger = logging.getLogger(__name__)


Sequence = Literal[0, 1, 2]

class MissingConnectionError(Exception):
    pass


# ------------------------------------------------------------------------------
# Data structure for emitted branches
# ------------------------------------------------------------------------------

@dataclass(frozen=True)
class BranchSpec:
    """
    Description of a short-circuit branch to be added via SCNetwork.add_branch(...).

    start/end are SC node IDs (where ground is represented by None).
    """
    Z_pu: complex
    start: str | None
    end: str | None
    has_source: bool = False


def is_ground_branch(bs: BranchSpec) -> bool:
    """Return True if the branch touches the reference node (ground)."""
    return (bs.start is None) ^ (bs.end is None)  # XOR: exactly one is None


def canonicalize_branch(bs: BranchSpec) -> BranchSpec:
    """
    Ensure a unique representation for ground-touching branches.

    Rule:
    -   If the branch touches ground, store it as (start=None, end=<bus>).
    -   If both ends are non-ground (series branch), leave as-is.
    -   If both ends are None (degenerate), keep as-is (will be ignored/handled
        elsewhere).
    """
    if not is_ground_branch(bs):
        return bs

    if bs.start is None:
        # Already canonical: None -> bus
        return bs

    # Swap to None -> bus
    return BranchSpec(
        Z_pu=bs.Z_pu,
        start=None,
        end=bs.start,
        has_source=bs.has_source,
    )

# ------------------------------------------------------------------------------
# Default hooks / policies
# ------------------------------------------------------------------------------

def _sc_volt_factor(cfg: SCCalcConfig) -> float:
    """Return the voltage factor for the configured short-circuit case."""
    case = cfg.sc_case.upper().strip()
    if case == "MAX":
        return cfg.volt_factor_max
    if case == "MIN":
        return cfg.volt_factor_min
    raise ValueError(f"sc_case must be 'MAX' or 'MIN', got {cfg.sc_case!r}")


def default_is_source_component(comp: Any, *, sequence: Sequence, cfg: SCCalcConfig) -> bool:
    """
    Decide whether a component is treated as a short-circuit source.

    - Positive sequence: Grid always, optionally machines.
    - Negative and zero sequence: sources are typically shorted/omitted -> False.
    """
    if sequence != 1:
        return False

    cls = comp.__class__.__name__

    if cls == "Grid":
        return True

    # Future extension points:
    if cls in {"SynchronousGenerator", "SynchronousMotor"}:
        return True

    if cls == "InductionMotor":
        if not cfg.include_induction_motor_sources:
            return False
        if hasattr(comp, "P_n"):
            try:
                return comp.P_n >= cfg.induction_motor_min_Pn
            except Exception:
                return True
        return True

    return False


def default_get_Z_ohm(comp: Component, *, sequence: Sequence, cfg: SCCalcConfig) -> Quantity:
    """
    Extract Z(sequence) in Ohm from a component.

    Expected methods in your codebase:
      - Grid.get_impedance(volt_factor=...) -> {0,1,2}
      - Transformer.get_impedance(volt_factor=...) -> {0,1,2}
      - BusBar.get_impedance(T) -> {0,1,2}
      - Cable.get_impedance(T) -> {0,1,2}
      - InductionMotor.get_impedance(...) -> {0,1,2}
    """
    cls = comp.__class__.__name__
    vf = _sc_volt_factor(cfg)

    case = cfg.sc_case.upper().strip()
    T_busbar = cfg.busbar_T_max if case == "MAX" else cfg.busbar_T_min
    T_cable = cfg.cable_T_max if case == "MAX" else cfg.cable_T_min

    if cls in {"Grid", "Transformer"}:
        Z = comp.get_impedance(vf)
        return Z[sequence]

    if cls == "BusBar":
        Z = comp.get_impedance(T_busbar)
        return Z[sequence]

    if cls == "Cable":
        Z = comp.get_impedance(T_cable)
        return Z[sequence]

    if cls == "InductionMotor":
        Z = comp.get_impedance()
        return Z[sequence]

    # Fallback: try get_impedance() without args
    if hasattr(comp, "get_impedance") and callable(getattr(comp, "get_impedance")):
        Z = comp.get_impedance()
        if isinstance(Z, dict) and sequence in Z:
            return Z[sequence]

    raise TypeError(
        f"Unsupported component type {cls!r} for Z{sequence} extraction. "
        "Extend default_get_Z_ohm(...) or pass a custom get_Z_ohm hook."
    )


def default_bus_u_base_from_object(bus: Any) -> Quantity | None:
    """If Bus has U_base attribute, return it."""
    if hasattr(bus, "U_base"):
        U = getattr(bus, "U_base")
        return U if U is not None else None
    return None


def u_base_for_component(comp: Any, *, default_u: Quantity | None) -> Quantity | None:
    """
    Return the most appropriate line-to-line base voltage for PU conversion.

    Conventions in this codebase:
    - Transformer.get_impedance(...) returns impedances referred to the secondary side,
      so we use Transformer.U_ls as its voltage base.
    - Grid/Cable/BusBar typically expose U_l.
    - Motors expose U_n.
    """
    cls = comp.__class__.__name__

    if cls == "Transformer" and hasattr(comp, "U_ls") and getattr(comp, "U_ls") is not None:
        return getattr(comp, "U_ls")

    if hasattr(comp, "U_l") and getattr(comp, "U_l") is not None:
        return getattr(comp, "U_l")

    if hasattr(comp, "U_n") and getattr(comp, "U_n") is not None:
        return getattr(comp, "U_n")

    for attr in ("U_nom", "U_nominal"):
        if hasattr(comp, attr) and getattr(comp, attr) is not None:
            return getattr(comp, attr)

    return default_u


def infer_bus_u_base(
    nw: NetworkGraph,
    *,
    bus_u_base_from_object: Callable[[Any], Quantity | None] = default_bus_u_base_from_object,
) -> dict[str, Quantity]:
    """
    Infer line-to-line base voltage (U_base) for each bus.

    Strategy (pragmatic MVP):
    - Use Bus.U_base if present.
    - Otherwise, BFS from GROUND and use voltage hints from components:
        - If a Transformer is present:
            * if upstream U matches U_lp, downstream U=U_ls
            * if upstream U matches U_ls, downstream U=U_lp
            * if upstream is GROUND and unknown, assume downstream is secondary (U_ls)
        - Otherwise, if any component has voltage attribute, use that.
    """
    ground = nw.GROUND_ID

    # Incidence map: bus_id -> list of conn_id
    inc: dict[str, list[str]] = {bus_id: [] for bus_id in nw.busses.keys()}
    for conn_id, conn in nw.connections.items():
        a = str(getattr(conn.start, "name", conn.start))
        b = str(getattr(conn.end, "name", conn.end))
        inc.setdefault(a, []).append(conn_id)
        inc.setdefault(b, []).append(conn_id)

    # Start with explicit U_base
    U_base: dict[str, Quantity] = {}
    for bus_id, bus in nw.busses.items():
        U = bus_u_base_from_object(bus)
        if U is not None:
            U_base[bus_id] = U

    visited: set[str] = {ground}
    queue: list[str] = [ground]

    def _voltage_hint(conn: Connection) -> Quantity | None:
        for comp in conn.components.values():
            U = u_base_for_component(comp, default_u=None)
            if U is not None:
                return U
        return None

    while queue:
        u = queue.pop(0)
        for conn_id in inc.get(u, []):
            conn = nw.connections[conn_id]
            v = other_bus_id(conn, u)
            if v in visited:
                continue

            inferred: Quantity | None = None

            # Check if the connection has a transformer
            transformer = None
            for comp in conn.components.values():
                if comp.__class__.__name__ == "Transformer":
                    transformer = comp
                    break

            if transformer is not None:
                transformer = cast(Transformer, transformer)
                U_lp = transformer.U_lp
                U_ls = transformer.U_ls
                if u in U_base:
                    Uu = U_base[u].to("V").m
                    if abs(Uu - U_lp.to("V").m) / max(U_lp.to("V").m, 1.0) < 1e-3:
                        inferred = U_ls
                    elif abs(Uu - U_ls.to("V").m) / max(U_ls.to("V").m, 1.0) < 1e-3:
                        inferred = U_lp
                else:
                    if u == ground:
                        inferred = U_ls  # convention for ground-connected source/trafo branch

            if inferred is None:
                inferred = _voltage_hint(conn)

            if inferred is not None:
                U_base[v] = inferred
                # Persist the inferred base voltage onto the bus object (useful downstream).
                nw.busses[v].U_base = inferred

            visited.add(v)
            queue.append(v)

    missing = [b for b in nw.busses.keys() if b != ground and b not in U_base]
    if missing:
        raise ValueError(
            "Could not infer U_base for all non-ground busses. "
            "Please set Bus.U_base for these busses or ensure voltage-hint components are connected. "
            f"Missing: {missing}"
        )

    return U_base


def ohm_to_pu(Z_ohm: Quantity, *, S_base: Quantity, U_base: Quantity) -> complex:
    """Convert an Ohmic impedance to per-unit using the provided PU base."""
    pu = PerUnitSystem(S_base=S_base, U_base=U_base)
    return complex(pu.get_per_unit_impedance(Z_ohm))


# ------------------------------------------------------------------------------
# Zero-sequence helpers (MVP)
# ------------------------------------------------------------------------------

def _norm_conn(value: Any) -> str | None:
    """Normalize winding connection values to 'Y', 'YN' or 'D' (upper-case)."""
    if value is None:
        return None
    s = str(value).strip().upper()
    return s


def _is_same_voltage(Ua: Quantity, Ub: Quantity, *, rel_tol: float = 1e-3) -> bool:
    """Return True if two line-to-line voltages are close within rel_tol."""
    a = Ua.to("V").m
    b = Ub.to("V").m
    return abs(a - b) / max(abs(a), abs(b), 1.0) <= rel_tol


def _find_transformers(conn: Connection) -> list[Transformer]:
    """Return all Transformer components on the connection."""
    out: list[Transformer] = []
    for comp in conn.components.values():
        if comp.__class__.__name__ == "Transformer":
            out.append(cast(Transformer, comp))
    return out


def _map_busses_to_transformer_sides(
    trafo: Transformer,
    bus_a: str,
    bus_b: str,
    *,
    U_base: dict[str, Quantity],
) -> tuple[str, str]:
    """
    Map endpoint busses to transformer primary/secondary sides based on base
    voltage.

    Returns (pri_bus, sec_bus).

    Raises if mapping is ambiguous.
    """
    Ua = U_base[bus_a]
    Ub = U_base[bus_b]

    U_lp = trafo.U_lp
    U_ls = trafo.U_ls

    # Determine which bus matches which transformer side by voltage.
    a_is_lp = _is_same_voltage(Ua, U_lp)
    a_is_ls = _is_same_voltage(Ua, U_ls)
    b_is_lp = _is_same_voltage(Ub, U_lp)
    b_is_ls = _is_same_voltage(Ub, U_ls)

    # Common case: one bus matches lp and the other matches ls.
    if a_is_lp and b_is_ls:
        return bus_a, bus_b
    if a_is_ls and b_is_lp:
        return bus_b, bus_a

    # If both busses are on the same voltage (e.g. modeling artifact), refuse
    # for MVP.
    raise ValueError(
        f"Cannot map busses {bus_a!r}-{bus_b!r} to transformer sides based on "
        f"voltage. Bus voltages: {Ua}, {Ub}; Transformer: U_lp={U_lp}, "
        f"U_ls={U_ls}."
    )


def _trafo_blocks_z0_transfer(trafo: Transformer) -> bool:
    """
    MVP rule for Z0 transfer:
    - Any delta winding blocks zero-sequence transfer between sides.
    """
    pri = _norm_conn(getattr(trafo, "pri_conn", None))
    sec = _norm_conn(getattr(trafo, "sec_conn", None))
    if pri is None or sec is None:
        raise ValueError(
            f"Transformer {getattr(trafo, 'name', '<unnamed>')!r}: "
            "pri_conn and sec_conn must be set for Z0 modelling."
        )
    return pri == "D" or sec == "D"


def _grounded_neutral_impedance(trafo: Transformer, side: str) -> Quantity | None:
    """
    Return the neutral-to-earth impedance for a grounded neutral, or None.

    side: 'pri' or 'sec'
    Condition: winding connection must be 'YN' and Zn must be provided.
    """
    if side not in {"pri", "sec"}:
        raise ValueError("side must be 'pri' or 'sec'")

    if side == "pri":
        conn = _norm_conn(getattr(trafo, "pri_conn", None))
        Zn = getattr(trafo, "Zn_pri", None)
    else:
        conn = _norm_conn(getattr(trafo, "sec_conn", None))
        Zn = getattr(trafo, "Zn_sec", None)

    if conn is None:
        raise ValueError(
            f"Transformer {getattr(trafo, 'name', '<unnamed>')!r}: "
            "pri_conn/sec_conn must be set for Z0 modelling."
        )

    if conn == "YN" and Zn is not None:
        return Zn.to("ohm")

    return None


# ------------------------------------------------------------------------------
# Builder implementation
# ------------------------------------------------------------------------------

class SequenceNetworkBuilder:
    """
    Builder object for sequence short-circuit networks.

    This class is a readability-focused refactor of the original
    build_sequence_network(...) function. Behavior is intended to be identical.
    """
    def __init__(
        self,
        nwgraph: NetworkGraph,
        *,
        sequence: Sequence,
        config: SCCalcConfig | None = None,
        # hooks:
        get_Z_ohm: Callable[[Any], Quantity] | None = None,
        is_source_component: Callable[[Any], bool] | None = None,
        bus_u_base: dict[str, Quantity] | None = None,
    ) -> None:
        self.nwgraph = nwgraph
        self.sequence = sequence
        self.cfg = config or SCCalcConfig()

        # Wrap hooks so they can see cfg/sequence without changing their
        # signature outside.
        self._get_Z_ohm = get_Z_ohm or (
            lambda comp: default_get_Z_ohm(
                comp,
                sequence=self.sequence,
                cfg=self.cfg
            )
        )
        self._is_source = is_source_component or (
            lambda comp: default_is_source_component(
                comp,
                sequence=self.sequence,
                cfg=self.cfg
            )
        )

        self.ground = self.nwgraph.GROUND_ID

        # incidence map: bus_id -> list of conn_id
        self.incidence: dict[str, list[str]] = self._build_incidence()

        if self.ground not in self.incidence:
            raise ValueError(
                f"GROUND_ID={self.ground!r} is not present "
                f"in network busses/incidence map."
            )

        # Determine U_base per bus
        self.U_base = bus_u_base or infer_bus_u_base(self.nwgraph)

        # The short-circuit network instance we are building
        self.nw_seq = SCNetwork()

        # Build-engine state
        self.processed_connections: set[str] = set()
        self.added_shunts: set[tuple[str, str]] = set()  # (conn_id, bus_id) to avoid duplicate shunts

        self.visited_buses: set[str] = {self.ground}
        self.queue: list[str] = []

    # --------------------------------------------------------------------------
    # Public API
    # --------------------------------------------------------------------------

    def build(self) -> SCNetwork:
        """Build and return the sequence Zbus short-circuit network."""
        if self.sequence in (1, 2):
            self._seed_seq12_from_ground()
        else:
            self._seed_seq0_from_ground_shunts()

        self._bfs_expand()

        self._assert_all_connections_processed()

        # Optional: add induction motor source branches ground->bus (sequence=1 only).
        if self.sequence == 1 and self.cfg.include_induction_motor_sources:
            self._add_induction_motor_sources()

        return self.nw_seq

    # --------------------------------------------------------------------------
    # Build engine
    # --------------------------------------------------------------------------

    def _seed_seq12_from_ground(self) -> None:
        """Seed BFS by adding all ground-touching branches for sequence 1/2."""
        root_conns: list[str] = []
        for conn_id, conn in self.nwgraph.connections.items():
            a = bus_id_of(conn.start)
            b = bus_id_of(conn.end)
            if a == self.ground or b == self.ground:
                root_conns.append(conn_id)

        # Add all ground-touching connections first
        for conn_id in root_conns:
            conn = self.nwgraph.connections[conn_id]
            a = bus_id_of(conn.start)
            b = bus_id_of(conn.end)

            if a == self.ground and b == self.ground:
                self.processed_connections.add(conn_id)
                continue

            non_ground_bus = b if a == self.ground else a

            Z_pu = self._z_pu_of_connection(conn, prefer_bus=non_ground_bus)
            self._add_branchspec(
                BranchSpec(
                    Z_pu=Z_pu,
                    start=None,
                    end=non_ground_bus,
                    has_source=self._conn_has_source(conn),
                )
            )
            self.processed_connections.add(conn_id)

            if non_ground_bus not in self.visited_buses:
                self.visited_buses.add(non_ground_bus)
                self.queue.append(non_ground_bus)

    def _seed_seq0_from_ground_shunts(self) -> None:
        """
        Sequence=0: add all grounding shunts (bus->ground) first, to seed BFS.

        If no grounding reference exists (no shunts found), Z0 Zbus is
        ill-defined.
        """
        for conn_id, conn in self.nwgraph.connections.items():
            u = bus_id_of(conn.start)
            v = bus_id_of(conn.end)

            for bs in self._emit_branches(conn, u, v):
                bs = canonicalize_branch(bs)

                if is_ground_branch(bs):
                    # Pre-seed only with ground branches.
                    if bs.end is None:
                        continue  # degenerate, should not happen

                    key = (conn_id, bs.end)
                    if key in self.added_shunts:
                        continue

                    self._add_branchspec(bs)
                    self.added_shunts.add(key)

                    if bs.end not in self.visited_buses:
                        self.visited_buses.add(bs.end)
                        self.queue.append(bs.end)

            # If this connection emits no series branch (delta-blocked etc.),
            # we will mark it as "done" later during BFS when encountered from a
            # visited bus. We do NOT mark it here because the series branch may
            # still exist.

        non_ground_busses = [
            b for b in self.nwgraph.busses.keys()
            if b != self.ground
        ]
        if non_ground_busses and not self.queue:
            raise ValueError(
                "Cannot build Z0 network: no grounding paths "
                "(bus->ground shunts) were found. Provide transformer "
                "neutral grounding (YN + Zn) or other grounding components."
            )

    def _bfs_expand(self) -> None:
        """BFS expansion for all sequences."""
        while self.queue:
            u = self.queue.pop(0)

            for conn_id in self.incidence.get(u, []):
                if conn_id in self.processed_connections:
                    continue

                conn = self.nwgraph.connections[conn_id]
                v = other_bus_id(conn, u)

                branches = self._emit_branches(conn, u, v)

                # 1) Add ground branches (shunts) first (avoid duplicates).
                for bs in branches:
                    bs = canonicalize_branch(bs)

                    if is_ground_branch(bs):
                        if bs.end is None:
                            continue

                        key = (conn_id, bs.end)
                        if key in self.added_shunts:
                            continue

                        self._add_branchspec(bs)
                        self.added_shunts.add(key)

                        if bs.end not in self.visited_buses:
                            self.visited_buses.add(bs.end)
                            self.queue.append(bs.end)

                # 2) Add the series branch (if any) to expand the graph.
                series_branches = [
                    canonicalize_branch(bs) for bs in branches
                    if not is_ground_branch(bs)
                ]
                if not series_branches:
                    # No series topology contribution (e.g. delta-blocked
                    # without transfer).
                    self.processed_connections.add(conn_id)
                    continue

                # For MVP, we expect at most one series branch per connection.
                if len(series_branches) > 1:
                    raise ValueError(
                        "Internal error: expected at most one series branch per "
                        f"connection. Got {len(series_branches)} for connection "
                        f"{conn_id!r}."
                    )

                bs = series_branches[0]
                self._add_branchspec(bs)
                self.processed_connections.add(conn_id)

                if v not in self.visited_buses:
                    self.visited_buses.add(v)
                    self.queue.append(v)

    def _assert_all_connections_processed(self) -> None:
        missing = set(self.nwgraph.connections.keys()) - self.processed_connections
        if missing:
            raise MissingConnectionError(
                "Not all connections could be added to the short-circuit network. "
                "Likely the network has a disconnected island not connected to the "
                "reference (GROUND) or no grounding exists for Z0. "
                f"Missing connections: {sorted(missing)}"
            )

    # --------------------------------------------------------------------------
    # SCNetwork write helpers
    # --------------------------------------------------------------------------

    def _add_branchspec(self, bs: BranchSpec) -> None:
        """Add a branch described by BranchSpec."""
        logger.debug(
            f"Add branch {bs.start}->{bs.end} with Z_pu = {bs.Z_pu:.3g} "
            f"to Z{self.sequence}-network."
        )

        self.nw_seq.add_branch(
            bs.Z_pu,
            start_node_ID=bs.start,
            end_node_ID=bs.end,
            has_source=bs.has_source
        )

    # --------------------------------------------------------------------------
    # Branch emission
    # --------------------------------------------------------------------------

    def _emit_branches(
        self,
        conn: Connection,
        u: str,
        v: str
    ) -> list[BranchSpec]:
        if self.sequence in (1, 2):
            return self._emit_branches_seq12(conn, u, v)
        return self._emit_branches_seq0(conn, u, v)

    def _emit_branches_seq12(
        self,
        conn: Connection,
        u: str,
        v: str
    ) -> list[BranchSpec]:
        """
        Emit branches for sequence 1/2.

        For sequence=1:
          - ground-touching connections may carry sources (has_source=True).
        For sequence=2:
          - sources are omitted by default policy (has_source=False).
        """
        Z_pu = self._z_pu_of_connection(conn, prefer_bus=u)
        return [
            BranchSpec(
                Z_pu=Z_pu,
                start=to_sc_node_id(u, self.ground),
                end=to_sc_node_id(v, self.ground),
                has_source=False,
            )
        ]

    def _emit_branches_seq0(
        self,
        conn: Connection,
        u: str,
        v: str
    ) -> list[BranchSpec]:
        """
        Emit branches for zero-sequence (sequence=0).

        Logic (MVP):
        -   If no transformer in the connection: emit a normal series branch
            u-v with Z0.
        -   If one transformer exists:
            *   If any delta winding exists -> no series transfer branch.
            *   If no delta -> emit series branch u-v with Z0.
            *   Additionally, if a side has grounded neutral (YN + Zn), emit
                shunt branch (ground -> bus).
        """
        transformers = _find_transformers(conn)
        if not transformers:
            Z_pu = self._z_pu_of_connection(conn, prefer_bus=u)
            return [BranchSpec(
                Z_pu=Z_pu,
                start=to_sc_node_id(u, self.ground),
                end=to_sc_node_id(v, self.ground))
            ]

        if len(transformers) > 1:
            raise ValueError(
                "Z0 builder MVP supports at most one Transformer per Connection. "
                f"Found {len(transformers)} transformers in connection between "
                f"{u!r} and {v!r}."
            )

        trafo = transformers[0]

        pri_bus, sec_bus = _map_busses_to_transformer_sides(trafo, u, v, U_base=self.U_base)
        blocks_transfer = _trafo_blocks_z0_transfer(trafo)

        branches: list[BranchSpec] = []

        if not blocks_transfer:
            Z_pu = self._z_pu_of_connection(conn, prefer_bus=u)
            branches.append(BranchSpec(
                Z_pu=Z_pu,
                start=to_sc_node_id(u, self.ground),
                end=to_sc_node_id(v, self.ground))
            )

        Zn_pri = _grounded_neutral_impedance(trafo, "pri")
        if Zn_pri is not None:
            Z_sh = self._z0_shunt_pu_for_side(conn, trafo, side_bus=pri_bus, Zn=Zn_pri)
            branches.append(BranchSpec(Z_pu=Z_sh, start=None, end=pri_bus))

        Zn_sec = _grounded_neutral_impedance(trafo, "sec")
        if Zn_sec is not None:
            Z_sh = self._z0_shunt_pu_for_side(conn, trafo, side_bus=sec_bus, Zn=Zn_sec)
            branches.append(BranchSpec(Z_pu=Z_sh, start=None, end=sec_bus))

        return branches

    # --------------------------------------------------------------------------
    # Electrical utilities
    # --------------------------------------------------------------------------

    def _build_incidence(self) -> dict[str, list[str]]:
        incidence: dict[str, list[str]] = {
            bus_id: [] for bus_id
            in self.nwgraph.busses.keys()
        }
        for conn_id, conn in self.nwgraph.connections.items():
            a = bus_id_of(conn.start)
            b = bus_id_of(conn.end)
            incidence.setdefault(a, []).append(conn_id)
            incidence.setdefault(b, []).append(conn_id)
        return incidence

    @staticmethod
    def _has_transformer(conn: Connection) -> bool:
        return any(
            comp.__class__.__name__ == "Transformer"
            for comp in conn.components.values()
        )

    def _z_pu_of_connection(
        self,
        conn: Connection,
        *,
        prefer_bus: str | None = None
    ) -> complex:
        """
        Compute total sequence impedance (pu) of a Connection by summing
        per-component impedances.

        When a transformer exists, conversion is done per component using
        u_base_for_component(). prefer_bus is used as a fallback to select a
        default U_base when endpoints have different U_base.
        """
        a = bus_id_of(conn.start)
        b = bus_id_of(conn.end)

        if a == self.ground and b != self.ground:
            U_conn = self.U_base[b]
        elif b == self.ground and a != self.ground:
            U_conn = self.U_base[a]
        else:
            Ua = self.U_base[a].to("V").m
            Ub = self.U_base[b].to("V").m
            if abs(Ua - Ub) / max(Ua, Ub, 1.0) > 1e-6 and not self._has_transformer(conn):
                raise ValueError(
                    "Connection links two different voltage levels but has no "
                    f"Transformer. Got {self.U_base[a]} vs {self.U_base[b]} for "
                    f"endpoints {a!r}-{b!r}."
                )
            U_conn = self.U_base[prefer_bus] if (prefer_bus is not None) else self.U_base[a]

        z_total: complex = 0 + 0j
        for comp in conn.components.values():
            if self.cfg.ignore_induction_motors_in_branch_Z and comp.__class__.__name__ == "InductionMotor":
                continue
            Z_ohm = self._get_Z_ohm(comp)
            U_comp = u_base_for_component(comp, default_u=U_conn) or U_conn
            z_total += ohm_to_pu(Z_ohm, S_base=self.cfg.S_base, U_base=U_comp)
        return z_total

    def _conn_has_source(self, conn: Connection) -> bool:
        """
        Source modeling policy:
        - Only for sequence=1 (positive)
        - Typically only on ground-touching branches
        """
        if self.sequence != 1:
            return False
        a = bus_id_of(conn.start)
        b = bus_id_of(conn.end)
        if not (a == self.ground or b == self.ground):
            return False
        return any(self._is_source(comp) for comp in conn.components.values())

    def _z0_shunt_pu_for_side(
        self,
        conn: Connection,
        trafo: Transformer,
        *,
        side_bus: str,
        Zn: Quantity
    ) -> complex:
        """
        Compute a Z0 shunt impedance (pu) from side_bus to ground.

        MVP approximation:
        -   Include Z0 of transformer.
        -   Include 3 * Zn (neutral grounding impedance contribution in zero-sequence).
        -   Include other non-transformer components in the same Connection that
            are on the same voltage level as side_bus (based on U_base matching).

        This yields a pragmatic "bus -> ground" branch impedance.
        """
        z_total: complex = 0 + 0j
        U_side = self.U_base[side_bus]

        # 1) Transformer Z0 contribution (always included)
        Z0_ohm = self._get_Z_ohm(trafo)
        U_trafo = u_base_for_component(trafo, default_u=U_side) or U_side
        z_total += ohm_to_pu(Z0_ohm, S_base=self.cfg.S_base, U_base=U_trafo)

        # 2) Grounding impedance contribution: 3 * Zn
        z_total += ohm_to_pu(3.0 * Zn, S_base=self.cfg.S_base, U_base=U_side)

        # 3) Add other components on the same voltage side (optional MVP behavior)
        for comp in conn.components.values():
            if comp is trafo:
                continue
            if self.cfg.ignore_induction_motors_in_branch_Z and comp.__class__.__name__ == "InductionMotor":
                continue

            U_comp = u_base_for_component(comp, default_u=None)
            if U_comp is None or _is_same_voltage(U_comp, U_side):
                Zc_ohm = self._get_Z_ohm(comp)
                Uc = U_comp or U_side
                z_total += ohm_to_pu(Zc_ohm, S_base=self.cfg.S_base, U_base=Uc)

        return z_total

    def _add_induction_motor_sources(self) -> None:
        for conn in self.nwgraph.connections.values():
            a = bus_id_of(conn.start)
            b = bus_id_of(conn.end)
            terminal_bus = a if a != self.ground else (b if b != self.ground else None)
            if terminal_bus is None:
                continue

            for comp in conn.components.values():
                if comp.__class__.__name__ != "InductionMotor":
                    continue
                if not self._is_source(comp):
                    continue

                Z_ohm = self._get_Z_ohm(comp)
                Z_pu = ohm_to_pu(Z_ohm, S_base=self.cfg.S_base, U_base=self.U_base[terminal_bus])
                self._add_branchspec(BranchSpec(Z_pu=Z_pu, start=None, end=terminal_bus, has_source=True))


def build_sequence_network(
    nwgraph: NetworkGraph,
    *,
    sequence: Sequence,
    config: SCCalcConfig | None = None,
    # hooks:
    get_Z_ohm: Callable[[Any], Quantity] | None = None,
    is_source_component: Callable[[Any], bool] | None = None,
    bus_u_base: dict[str, Quantity] | None = None
) -> SCNetwork:
    """
    Build the sequence Zbus network.

    Parameters
    ----------
    nwgraph:
        Network graph (busses + connections + components).
    sequence:
        1=positive, 2=negative, 0=zero.
    config:
        CalcConfig controlling base power, factors and temperatures.
    get_Z_ohm:
        Optional hook to extract component impedance for this sequence (Ohm).
        If None, default_get_Z_ohm(...) is used.
    is_source_component:
        Optional hook to detect sources. If None, default_is_source_component()
        is used.
    bus_u_base:
        Optional explicit mapping bus_id -> U_base (line-to-line).
        If None, inferred.

    Notes on sequence=0
    -------------------
    For Z0, topology depends on transformer winding connection and neutral
    grounding. This function therefore:
    (1)   adds all grounding (bus->ground) shunts first,
    (2)   then runs BFS to add remaining series branches where Z0 transfer
          exists.

    Returns
    -------
    SCNetwork
        A short-circuit network with branches added in a valid order.
    """
    return SequenceNetworkBuilder(
        nwgraph,
        sequence=sequence,
        config=config,
        get_Z_ohm=get_Z_ohm,
        is_source_component=is_source_component,
        bus_u_base=bus_u_base,
    ).build()

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

def bus_id_of(endpoint: Any) -> str:
    """Return the bus id from a connection endpoint."""
    return str(getattr(endpoint, "name", endpoint))

def other_bus_id(conn: Connection, bus_id: str) -> str:
    """Given one endpoint bus_id, return the other endpoint bus_id."""
    a = bus_id_of(conn.start)
    b = bus_id_of(conn.end)
    if a == bus_id:
        return b
    if b == bus_id:
        return a
    raise ValueError(f"Bus {bus_id!r} is not an endpoint of this connection.")

def to_sc_node_id(bus_id: str, ground_id: str) -> str | None:
    """Map python-electric bus ids to short_circuit node ids (GROUND -> None)."""
    return None if bus_id == ground_id else bus_id
