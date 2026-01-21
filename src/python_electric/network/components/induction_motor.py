from dataclasses import dataclass, field

from python_electric import Quantity
from python_electric.network.network import Component

__all__ = ["InductionMotor"]


@dataclass
class InductionMotor(Component):
    name: str
    U_n: Quantity
    I_n: Quantity
    P_m: Quantity
    eta: Quantity
    cos_phi: float
    k_start: float = 6.0
    R_to_X: float = 0.42
    z2_factor: float = 1.0
    z0_factor: float = 10.0

    Z_dict: dict[int, Quantity[float | complex]] = field(init=False, default_factory=dict)

    def __post_init__(self):
        super().__init__(self.name)
        self.Z_dict = self.get_impedance()

    def get_impedance(self) -> dict[int, Quantity]:
        from python_electric.equipment import InductionMotor
        motor = InductionMotor(
            nominal_voltage=self.U_n,
            nominal_current=self.I_n,
            locked_rotor_current=self.k_start * self.I_n,
            P_m=self.P_m,
            efficiency=self.eta,
            power_factor=self.cos_phi,
            R_to_X_ratio=self.R_to_X,
            z2_factor=self.z2_factor,
            z0_factor=self.z0_factor,
            name=self.name
        )
        return {1: motor.Z1, 2: motor.Z2, 0: motor.Z0}
