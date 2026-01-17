from typing import TypeAlias, Sequence

from ... import Q_
from ... import equipment
from ... import network
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
    network.Transformer
    | network.BusBar
    | network.Cable
    | network.Grid
    | network.InductionMotor
)


Component: TypeAlias = NetworkEquipment | NetworkComponent


_config_dict = {
    "MAX": {"T": Q_(20, 'degC'), "c": 1.1},
    "MIN": {"T": Q_(150, 'degC'), "c": 0.95}
}


def get_impedances(
    components: Sequence[tuple[Component, PerUnitSystem]],
    sc_case: str
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
    sc_case: str, {"MAX", "MIN"}
        Indicates if the impedances will be used to calculate the maximum
        short-circuit current (sc_case="MAX") or the minimum
        short-circuit current (sc_case="MIN").

    Notes
    -----
    `Component` is type alias for all classes representing network components:
    1.  from package `equipment`: `PowerGrid`, `Transformer`, `Cable`, `BusBar`,
        `Generator`, `SynchronousMotor`, and `InductionMotor`.
    2.  from package `network.components`: `Transformer`, `BusBar`, `Cable`.

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
        cfg = _config_dict[sc_case.upper()]
    except KeyError:
        raise ValueError(
            f"Parameter `sc_case` has unknown value: {sc_case}."
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
            if isinstance(component, (network.BusBar, network.Cable)):
                Z_012 = component.get_impedance(T=cfg["T"])
                Z0 = Z_012[0]
                Z1 = Z_012[1]
                Z2 = Z_012[2]
            else:
                try:
                    Z_012 = component.get_impedance(volt_factor=cfg["c"])
                except TypeError:
                    Z_012 = component.get_impedance()
                Z0 = Z_012[0]
                Z1 = Z_012[1]
                Z2 = Z_012[2]
        else:
            raise ValueError("Got an unknown component.")
        Z_dict[component.name] = (Z0, Z1, Z2)
        Z_pu_dict[component.name] = tuple([
            per_unit_system.get_per_unit_impedance(Z)
            for Z in (Z0, Z1, Z2)
        ])
    return Z_dict, Z_pu_dict
