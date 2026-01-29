"""
Tkinter GUI helpers for building python-electric projects.

This subpackage provides a minimal MVP GUI to:
- define network topology (source + connections)
- assign component input data to connections
- export/load the project model (pickle)
- build a NetworkTopology instance for use in notebooks
"""
from .models import ProjectModel, SourceModel, ConnectionModel
from .build import build_network_topology
