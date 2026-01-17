from dataclasses import dataclass

from ... import Quantity, Q_
from ... import calc

__all__ = ["Load"]


@dataclass
class Load:
    U_lv: Quantity
    cos_phi: float
    P_e: Quantity | None = None
    P_m: Quantity | None = None
    eta: Quantity = Q_(100, 'pct')
    I_load: Quantity | None = None

    def get_load_current(self) -> Quantity:
        if isinstance(self.I_load, Quantity):
            return self.I_load
        if isinstance(self.P_e, Quantity):
            self.I_load = calc.load_current(self.U_lv, self.cos_phi, self.P_e)
            return self.I_load
        if isinstance(self.P_m, Quantity):
            self.I_load = calc.load_current(self.U_lv, self.cos_phi, P_m=self.P_m, eta=self.eta)
            return self.I_load
        raise ValueError("Load is undefined.")