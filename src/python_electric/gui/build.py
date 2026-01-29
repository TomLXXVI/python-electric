from __future__ import annotations

from .. import Q_
from ..network.topology import NetworkTopology
from ..network.components import Load
from .models import ProjectModel


def _ensure_load(project: ProjectModel, conn_id: str) -> Load:
    conn = project.connections[conn_id]
    if conn.load is not None:
        return conn.load
    # Minimal default load (so components like Cable can compute I_b, etc.)
    return Load(U_l=project.U_n, cos_phi=1.0, P_e=Q_(0, "W"))


def build_network_topology(project: ProjectModel) -> NetworkTopology:
    """Build a NetworkTopology instance from the GUI project model."""
    topo = NetworkTopology(
        name=project.name,
        U_n=project.U_n,
        earthing_system=project.earthing_system,
        neutral_distributed=project.neutral_distributed,
    )

    # Source (grid) connection
    s = project.source
    topo.add_source_grid_connection(
        conn_id=s.conn_id,
        end_id=s.end_id,
        U_l=s.U_l,
        S_sc=s.S_sc,
        R_to_X=s.R_to_X,
        z0_r_factor=s.z0_r_factor,
        z0_x_factor=s.z0_x_factor,
    )

    # Regular connections and their loads
    for conn_id, conn in project.connections.items():
        topo.add_connection(conn_id, conn.start_id, conn.end_id, load=_ensure_load(project, conn_id))

    # Components on each connection
    for conn_id, conn in project.connections.items():
        for comp in conn.components.values():
            topo.add_component(conn_id, comp)

    return topo
