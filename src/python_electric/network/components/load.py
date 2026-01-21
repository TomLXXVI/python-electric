from dataclasses import dataclass

from ... import Quantity, Q_
from ... import calc

__all__ = ["Load"]


@dataclass
class Load:
    U_l: Quantity
    cos_phi: float
    P_e: Quantity | None = None
    P_m: Quantity | None = None
    eta: Quantity = Q_(100, 'pct')
    I_b: Quantity | None = None

    def __post_init__(self):
        self.I_b = self._get_load_current()

    def _get_load_current(self) -> Quantity:
        if isinstance(self.I_b, Quantity):
            return self.I_b
        if isinstance(self.P_e, Quantity):
            self.I_b = calc.load_current(self.U_l, self.cos_phi, self.P_e)
            return self.I_b
        if isinstance(self.P_m, Quantity):
            self.I_b = calc.load_current(self.U_l, self.cos_phi, P_m=self.P_m, eta=self.eta)
            return self.I_b
        raise ValueError("Load is undefined.")