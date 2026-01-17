from dataclasses import dataclass, field

from python_electric import Quantity
from python_electric.network.network import Component

__all__ = ["Grid"]


@dataclass
class Grid(Component):
    """
    Represents the three-phase power grid to which the low-voltage network is
    connected.

    Parameters
    ----------
    name: str, optional
        Uniquely identifies the grid in the network.
    U_l: Quantity
        Rated line-to-line voltage of the grid.
    S_sc: Quantity
        Short-circuit power at the entry point of the low-voltage network.

    Attributes
    ----------
    R_to_X:
        Ratio of resistance to reactance of the power grid short-circuit
        impedance.
    Z_dict: dict[int, Quantity]
        Zero-sequence (key: 0), positive sequence (key: 1), and negative
        sequence (key: 2) impedances of the transformer referred to the
        secondary side.
    """
    name: str
    U_l: Quantity
    S_sc: Quantity

    R_to_X: float = 0.1

    Z_dict: dict[int, Quantity[float | complex]] = field(init=False, default_factory=dict)

    def __post_init__(self):
        super().__init__(self.name)
        self.Z_dict = self.get_impedance()

    def get_impedance(
        self,
        volt_factor: float = 1.1,
        z0_r_factor: float = 1.0,
        z0_x_factor: float = 3.0
    ) -> dict[int, Quantity]:
        from python_electric.equipment import PowerGrid
        grid = PowerGrid(
            line_voltage=self.U_l,
            short_circuit_power=self.S_sc,
            R_to_X_ratio=self.R_to_X,
            voltage_factor=volt_factor,
            z0_r_factor=z0_r_factor,
            z0_x_factor=z0_x_factor,
            name=self.name
        )
        return {1: grid.Z1, 2: grid.Z2, 0: grid.Z0}
