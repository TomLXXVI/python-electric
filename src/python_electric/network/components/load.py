from __future__ import annotations
from dataclasses import dataclass, field

from ... import Quantity, Q_
from ... import calc

__all__ = ["Load"]


@dataclass
class Load:
    U_l: Quantity | None = None
    cos_phi: float = 0.8
    P_e: Quantity | None = None
    P_m: Quantity | None = None
    eta: Quantity = Q_(100, 'pct')
    k_u: float = 1.0  # utilization factor, ratio of actual used power to nominal power

    _I_b: Quantity = field(init=False, default=None)

    def _get_load_current(self) -> Quantity:
        if isinstance(self._I_b, Quantity):
            return self._I_b
        if isinstance(self.P_e, Quantity):
            self._I_b = calc.load_current(self.U_l, self.cos_phi, self.k_u * self.P_e)
            return self._I_b
        if isinstance(self.P_m, Quantity):
            self._I_b = calc.load_current(self.U_l, self.cos_phi, P_m=self.k_u * self.P_m, eta=self.eta)
            return self._I_b
        raise ValueError("Load is undefined.")

    @property
    def I_b(self) -> Quantity:
        if self._I_b is None:
            self._I_b = self._get_load_current()
        return self._I_b

    @I_b.setter
    def I_b(self, v: Quantity) -> None:
        self._I_b = v

    @classmethod
    def add(cls, *loads: Load) -> Load:
        I_b = sum([load.I_b for load in loads])
        load = cls()
        load.I_b = I_b
        return load
