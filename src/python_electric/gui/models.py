from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Type

import pickle

from .. import Quantity, Q_
from ..network.topology import CableInput, TransformerInput
from ..network.components import Load
from ..protection.earthing_system import EarthingSystem


SCHEMA_VERSION = 1


@dataclass
class SourceModel:
    """User data for the source (grid) connection."""
    conn_id: str = "C0"
    end_id: str = "B1"
    U_l: Quantity = Q_(400, "V")
    S_sc: Quantity = Q_(500, "MVA")
    R_to_X: float = 0.1
    z0_r_factor: float = 0.0
    z0_x_factor: float = 0.0


ComponentInputType = Type[Any]  # CableInput, TransformerInput, ...


@dataclass
class ConnectionModel:
    """User data for a connection in the network."""
    conn_id: str
    start_id: str
    end_id: str
    load: Load | None = None
    # Map component input class -> instance (e.g. CableInput -> CableInput(...))
    components: dict[ComponentInputType, Any] = field(default_factory=dict)


@dataclass
class ProjectModel:
    """GUI project container that can be saved/loaded and built into NetworkTopology."""
    name: str = "network"
    U_n: Quantity = Q_(400, "V")
    earthing_system: EarthingSystem = EarthingSystem.TN
    neutral_distributed: bool = True

    source: SourceModel = field(default_factory=SourceModel)
    connections: dict[str, ConnectionModel] = field(default_factory=dict)

    # Defaults for component input dataclasses
    defaults: dict[ComponentInputType, Any] = field(default_factory=dict)

    # Metadata
    schema_version: int = SCHEMA_VERSION
    python_electric_version: str | None = None

    def __post_init__(self) -> None:
        # Defaults only created if not present (so loads from file keep their defaults)
        if not self.defaults:
            self.defaults = self._make_default_components()

    @staticmethod
    def _make_default_components() -> dict[ComponentInputType, Any]:
        # Keep MVP small: Cable + Transformer defaults.
        cable = CableInput(name="cable", L=Q_(1, "m"))
        # Provide transformer defaults that are common in examples; user can override in GUI.
        transformer = TransformerInput(
            name="transformer",
            S_n=Q_(630, "kVA"),
            U_lp=Q_(10, "kV"),
            U_ls=Q_(400, "V"),
            u_cc=Q_(6, "pct"),
            P_Cu=Q_(6, "kW"),
        )
        return {CableInput: cable, TransformerInput: transformer}

    # ------------ convenience ------------
    def add_connection(self, conn_id: str, start_id: str, end_id: str) -> None:
        if conn_id in self.connections:
            raise ValueError(f"Connection '{conn_id}' already exists.")
        self.connections[conn_id] = ConnectionModel(conn_id=conn_id, start_id=start_id, end_id=end_id)

    def delete_connection(self, conn_id: str) -> None:
        self.connections.pop(conn_id, None)

    def list_connection_ids(self) -> list[str]:
        # return sorted(self.connections.keys())
        return list(self.connections.keys())

    def set_component(self, conn_id: str, comp_input: Any) -> None:
        if conn_id not in self.connections:
            raise KeyError(f"Connection '{conn_id}' not found.")
        self.connections[conn_id].components[type(comp_input)] = comp_input

    def get_component(self, conn_id: str, comp_cls: ComponentInputType) -> Any | None:
        if conn_id not in self.connections:
            return None
        return self.connections[conn_id].components.get(comp_cls)

    def set_defaults(self, comp_input: Any) -> None:
        self.defaults[type(comp_input)] = comp_input

    def get_defaults(self, comp_cls: ComponentInputType) -> Any:
        return self.defaults[comp_cls]

    # ------------ persistence ------------
    def save(self, path: str) -> None:
        payload = {
            "schema_version": self.schema_version,
            "python_electric_version": self.python_electric_version,
            "project": self,
        }
        with open(path, "wb") as f:
            # noinspection PyTypeChecker
            pickle.dump(payload, f)

    @classmethod
    def load(cls, path: str) -> "ProjectModel":
        with open(path, "rb") as f:
            payload = pickle.load(f)
        project = payload.get("project")
        if not isinstance(project, ProjectModel):
            raise TypeError("Invalid project file.")
        # Schema check (soft)
        sv = payload.get("schema_version", None)
        if sv != SCHEMA_VERSION:
            # keep it simple: accept but store the number
            project.schema_version = sv
        return project
