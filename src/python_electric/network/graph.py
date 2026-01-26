from __future__ import annotations

from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from typing import Type, TypeVar

from python_electric import Quantity

__all__ = ["Bus", "Connection", "Component", "NetworkGraph", "TComponent"]


@dataclass(slots=True)
class Bus:
    """A node in the network graph."""
    name: str

    U_base: Quantity = field(init=False, default=None)


@dataclass(slots=True)
class Connection:
    """
    A directed connection in the network graph (start -> end).
    """
    name: str
    start: Bus
    end: Bus
    components: dict[str, Component] = field(default_factory=dict)

    def add_component(self, comp: Component) -> None:
        if not comp.name:
            raise ValueError("Component.name must be set before adding it to a Connection.")

        if comp.name in self.components:
            raise KeyError(f"A component with id '{comp.name}' is already present on connection '{self.name}'.")

        self.components[comp.name] = comp
        comp.connection = self  # back-reference

    def remove_component(self, comp_id: str) -> Component:
        comp = self.components.pop(comp_id, None)
        if comp is not None:
            comp.connection = None
        return comp

    def get_component(self, comp_id: str) -> Component:
        comp = self.components.get(comp_id)
        if comp is None:
            raise KeyError(f"Component with id '{comp_id}' is not found on connection '{self.name}'.")
        return comp


class Component(ABC):
    """
    Base class for physical/electrical assets placed on a Connection.
    @DynamicAttrs
    """
    def __init__(self, name: str = "") -> None:
        self.name: str = name
        self.connection: Connection | None = None  # set when added to a Connection

    def __str__(self) -> str:
        return f"Component<{self.name}>"

    @abstractmethod
    def get_impedance(self, *args, **kwargs):
        ...


TComponent = TypeVar("TComponent", bound=Component)


class NetworkGraph:
    """
    Aggregate root. This class owns:
    - all busses
    - all connections

    Users should construct the model via NetworkGraph methods.
    """
    GROUND_ID = "ground"

    def __init__(self, name: str = "") -> None:
        self.name: str = name
        self._busses: dict[str, Bus] = {self.GROUND_ID: Bus(self.GROUND_ID)}
        self._connections: dict[str, Connection] = {}

    @classmethod
    def create(cls, name: str) -> NetworkGraph:
        return cls(name=name)

    @property
    def busses(self) -> dict[str, Bus]:
        return self._busses

    @property
    def connections(self) -> dict[str, Connection]:
        return self._connections

    def get_or_create_bus(self, bus_id: str) -> Bus:
        if not bus_id:
            raise ValueError("bus_id must be a non-empty string.")
        bus = self._busses.get(bus_id)
        if bus is None:
            bus = Bus(bus_id)
            self._busses[bus_id] = bus
        return bus

    def add_connection(self, conn_id: str, start_id: str, end_id: str) -> Connection:
        """
        Create and register a connection. Returns the created connection so
        users can chain operations if they want.
        """
        if not conn_id:
            raise ValueError("conn_id must be a non-empty string.")
        if conn_id in self._connections:
            raise ValueError(f"Network already contains a connection '{conn_id}'.")

        start = self.get_or_create_bus(start_id)
        end = self.get_or_create_bus(end_id)

        conn = Connection(name=conn_id, start=start, end=end)
        self._connections[conn_id] = conn
        return conn

    def get_connection(self, conn_id: str) -> Connection:
        conn = self._connections.get(conn_id)
        if conn is None:
            raise KeyError(f"Connection '{conn_id}' not found.")
        return conn

    def add_component(self, conn_id: str, comp: Component) -> None:
        """
        Centralized way to place components (preferred over calling
        conn.add_component directly).
        """
        self.get_connection(conn_id).add_component(comp)

    def get_component(self, conn_id: str, comp_id: str) -> Component:
        return self.get_connection(conn_id).get_component(comp_id)

    def find_components(
        self,
        comp_type: Type[TComponent],
        *,
        include_subclasses: bool = True,
    ) -> list[tuple[TComponent, str, str, str]]:
        """
        Find all components of a given type in the network.

        Returns
        -------
        list of tuples:
            (component, connection_id, start_id, end_id)
        """
        matches: list[tuple[TComponent, str, str, str]] = []

        for conn_id, conn in self._connections.items():
            start_id = conn.start.name
            end_id = conn.end.name

            for comp in conn.components.values():
                ok = isinstance(comp, comp_type) if include_subclasses else (type(comp) is comp_type)
                if ok:
                    matches.append((comp, conn_id, start_id, end_id))

        return matches

    def get_source_bus_ids(self) -> list[str]:
        """
        Return a list with the end bus IDs of all connections that start from
        ground.
        """
        source_bus_ids = [
            conn.end.name for conn in [
                conn for conn in self.connections.values()
                if conn.start.name == self.GROUND_ID
            ]
        ]
        return source_bus_ids

    def find_connection_path(
        self,
        start_bus_id: str,
        end_bus_id: str
    ) -> list[str]:
        """
        Find a directed path (sequence of connection IDs) from a start bus to an
        end bus.

        The method treats the network as a directed graph where every connection
        points from ``start`` to ``end``. The returned list contains the
        connection IDs in the order they must be traversed to go from
        ``start_bus_id`` to ``end_bus_id``.

        Parameters
        ----------
        start_bus_id:
            ID/name of the bus where the path starts.
        end_bus_id:
            ID/name of the bus where the path ends.

        Returns
        -------
        list[str]
            Ordered list of connection IDs forming the path. If ``start_bus_id``
            equals ``end_bus_id``, an empty list is returned.

        Raises
        ------
        KeyError
            If either bus does not exist in the network.
        ValueError
            If no directed path exists between the busses.
        """
        if start_bus_id not in self._busses:
            raise KeyError(f"Bus '{start_bus_id}' not found.")
        if end_bus_id not in self._busses:
            raise KeyError(f"Bus '{end_bus_id}' not found.")
        if start_bus_id == end_bus_id:
            return []

        # Build adjacency: start_bus -> [(conn_id, end_bus)]
        adjacency: dict[str, list[tuple[str, str]]] = {}
        for conn_id, conn in self._connections.items():
            adjacency.setdefault(conn.start.name, []).append((conn_id, conn.end.name))

        # Breadth-first search to find the shortest path in number of edges.
        from collections import deque

        queue: deque[str] = deque([start_bus_id])
        visited: set[str] = {start_bus_id}
        # For each visited bus (except start): store (previous_bus, via_connection_id)
        prev: dict[str, tuple[str, str]] = {}

        while queue:
            bus = queue.popleft()
            for conn_id, nxt in adjacency.get(bus, []):
                if nxt in visited:
                    continue
                visited.add(nxt)
                prev[nxt] = (bus, conn_id)
                if nxt == end_bus_id:
                    queue.clear()
                    break
                queue.append(nxt)

        if end_bus_id not in prev:
            raise ValueError(
                f"No directed path found from bus '{start_bus_id}' to bus '{end_bus_id}'."
            )

        # Reconstruct path (reverse)
        path_conn_ids: list[str] = []
        cur = end_bus_id
        while cur != start_bus_id:
            p_bus, p_conn = prev[cur]
            path_conn_ids.append(p_conn)
            cur = p_bus
        path_conn_ids.reverse()
        return path_conn_ids

    def find_reachable_final_busses(
        self,
        source_bus_id: str,
        *,
        include_ground: bool = False
    ) -> list[Bus]:
        """
        Find all *reachable* final busses starting from ``source_bus_id``.

        A "final bus" (leaf node) is defined as a bus with no outgoing
        connections (out-degree = 0) *within the directed network graph*.

        Parameters
        ----------
        source_bus_id:
            ID/name of the bus where the search starts.
        include_ground:
            If False (default), the special ground bus is excluded from the
            result.

        Returns
        -------
        list[Bus]
            List of reachable busses that have no outgoing connections.
            If the source bus itself has no outgoing connections, it will be
            included (unless it is the ground bus and include_ground=False).

        Raises
        ------
        KeyError
            If ``source_bus_id`` does not exist in the network.
        """
        if source_bus_id not in self._busses:
            raise KeyError(f"Bus '{source_bus_id}' not found.")

        # Build adjacency: bus -> next busses
        adjacency: dict[str, list[str]] = {}
        for conn in self._connections.values():
            adjacency.setdefault(conn.start.name, []).append(conn.end.name)

        from collections import deque
        q = deque([source_bus_id])
        visited = {source_bus_id}
        while q:
            b = q.popleft()
            for nxt in adjacency.get(b, []):
                if nxt not in visited:
                    visited.add(nxt)
                    q.append(nxt)

        start_bus_ids = {conn.start.name for conn in self._connections.values()}

        finals: list[Bus] = []
        for bus_id in visited:
            if not include_ground and bus_id == self.GROUND_ID:
                continue
            if bus_id not in start_bus_ids:
                finals.append(self._busses[bus_id])
        return finals

    @property
    def ground(self) -> Bus:
        return self._busses[self.GROUND_ID]

    def __str__(self) -> str:
        return f"Network<{self.name}>"
