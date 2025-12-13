"""
python_electric

Core toolbox for designing electrical installations in buildings.
"""
from .pint_setup import UNITS, Quantity
from .general import VoltReference

from . import calc
from . import equipment
from . import materials
from . import sizing
from . import short_circuit
from . import general


__all__ = [
    "UNITS",
    "Quantity",
    "equipment",
    "materials",
    "sizing",
    "short_circuit",
    "calc",
    "general"
]


__version__ = "0.1.0"
