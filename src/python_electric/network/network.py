from __future__ import annotations

from dataclasses import dataclass, field
from abc import ABC, abstractmethod

from python_electric import Quantity

__all__ = ["Bus", "Connection", "Component", "Network"]


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
    """
    def __init__(self, name: str = "") -> None:
        self.name: str = name
        self.connection: Connection | None = None  # set when added to a Connection

    def __str__(self) -> str:
        return f"Component<{self.name}>"

    @abstractmethod
    def get_impedance(self, *args, **kwargs):
        ...

class Network:
    """
    Aggregate root. This class owns:
    - all busses
    - all connections

    Users should construct the model via Network methods.
    """
    GROUND_ID = "ground"

    def __init__(self, name: str = "") -> None:
        self.name: str = name
        self._busses: dict[str, Bus] = {self.GROUND_ID: Bus(self.GROUND_ID)}
        self._connections: dict[str, Connection] = {}

    @classmethod
    def create(cls, name: str) -> Network:
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

    @property
    def ground(self) -> Bus:
        return self._busses[self.GROUND_ID]

    def __str__(self) -> str:
        return f"Network<{self.name}>"
