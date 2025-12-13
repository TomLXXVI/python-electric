from typing import TypeAlias, Sequence

from ... import equipment
from ... import sizing
from .per_unit import PerUnitSystem

__all__ = [
    "get_impedances",
]


NetworkEquipment: TypeAlias = (
    equipment.PowerGrid
    | equipment.Transformer
    | equipment.Cable
    | equipment.BusBar
    | equipment.Generator
    | equipment.SynchronousMotor
    | equipment.InductionMotor
)


NetworkComponent: TypeAlias = (
    sizing.Transformer
    | sizing.BusBars
    | sizing.SinglePhaseCable
    | sizing.ThreePhaseCable
)


Component: TypeAlias = NetworkEquipment | NetworkComponent


_which_mapping = {"MAX": "T_20", "MIN": "T_150"}


def get_impedances(
    components: Sequence[tuple[Component, PerUnitSystem]],
    which_short_circuit: str
) -> tuple[dict[str, tuple], dict[str, tuple]]:
    """
    Given a sequence of NetworkEquipment or NetworkComponent objects together
    with their associated per-unit system, returns the impedances and per-unit
    impedances of these components.

    Parameters
    ----------
    components: Sequence[tuple[Component, PerUnitSystem]]
        Each element in the sequence is actually a tuple composed of a network
        component and its associated per-unit system. Through the per-unit
        system, the per-unit impedance of each component in the sequence can be
        calculated.
    which_short_circuit: str, {"MAX", "MIN"}
        Indicates if the impedances will be used to calculate the maximum
        short-circuit current (which_short_circuit="MAX") or the minimum
        short-circuit current (which_short_circuit="MIN").

    Notes
    -----
    `Component` is type alias for all classes representing network components:
    1.  from package `equipment`: `PowerGrid`, `Transformer`, `Cable`, `BusBar`,
        `Generator`, `SynchronousMotor`, and `InductionMotor`.
    2.  from package `sizing`: `Transformer`, `BusBars`, `ThreePhaseCable`, and
        `SinglePhaseCable`.

    Returns
    -------
    Z_dict: dict[str, tuple]
        Dictionary of which the keys are the names of the network component
        objects and the values are the impedances of these objects.
    Z_pu_dict: dict[str, tuple]
        Dictionary of which the keys are the names of the network component
        objects and the values are the per-unit impedances of these objects.
    """
    try:
        which_short_circuit = _which_mapping[which_short_circuit.upper()]
    except KeyError:
        raise ValueError(
            f"Parameter `which` has unknown value: {which_short_circuit}."
            f"Should be either 'MAX' or 'MIN'."
        )
    Z_dict, Z_pu_dict = {}, {}
    for tup in components:
        component, per_unit_system = tup
        if isinstance(component, NetworkEquipment):
            Z0 = component.Z0
            Z1 = component.Z1
            Z2 = component.Z2
        elif isinstance(component, NetworkComponent):
            if isinstance(component, (sizing.BusBars, sizing.AbstractCable)):
                Z0 = component.Z[which_short_circuit][0]
                Z1 = component.Z[which_short_circuit][1]
                Z2 = component.Z[which_short_circuit][2]
            else:
                Z0 = component.Z[0]
                Z1 = component.Z[1]
                Z2 = component.Z[2]
        else:
            raise ValueError("Got an unknown component.")
        Z_dict[component.name] = (Z0, Z1, Z2)
        Z_pu_dict[component.name] = tuple([
            per_unit_system.get_per_unit_impedance(Z)
            for Z in (Z0, Z1, Z2)
        ])
    return Z_dict, Z_pu_dict
