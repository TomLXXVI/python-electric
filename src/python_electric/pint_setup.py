from __future__ import annotations

import pint
from pint.facets.plain.quantity import PlainQuantity as Quantity

UNITS = pint.UnitRegistry()

Q_ = UNITS.Quantity

unit_definitions = [
    'fraction = [] = frac',
    'percent = 1e-2 frac = % = pct',
    'parts_per_million = 1e-6 fraction = ppm'
]
for ud in unit_definitions:
    UNITS.define(ud)

pint.set_application_registry(UNITS)

__all__ = ["UNITS", "Q_", "Quantity"]
