from dataclasses import dataclass, field

from python_electric import Quantity, Q_
from python_electric.network.graph import Component

__all__ = ["Generator", "SynchronousMotor"]


@dataclass
class Generator(Component):
    name: str
    U_n: Quantity
    S_n: Quantity
    x: float
    cos_phi: float = 0.8
    R_to_X: float = 0.15
    z0_factor: float = 1.0

    Z_dict: dict[int, Quantity[float | complex]] = field(init=False, default_factory=dict)

    def __post_init__(self):
        super().__init__(self.name)
        self.Z_dict = self.get_impedance()

    def get_impedance(self, volt_factor: float = 1.0) -> dict[int, Quantity]:
        from python_electric.equipment import Generator
        gen = Generator(
            nominal_voltage=self.U_n,
            nominal_power=self.S_n,
            per_unit_reactance=Q_(self.x, 'frac'),
            power_factor=self.cos_phi,
            R_to_X_ratio=self.R_to_X,
            voltage_factor=volt_factor,
            z0_factor=self.z0_factor
        )
        return {1: gen.Z1, 2: gen.Z2, 0: gen.Z0}



class SynchronousMotor(Generator):
    pass
