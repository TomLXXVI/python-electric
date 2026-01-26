import math
from dataclasses import dataclass, field

from ... import Quantity, Q_, VoltReference
from ...materials import ConductorMaterial
from ...protection import CircuitBreaker
from ...calc import voltage_drop

from python_electric.network.graph import Component

__all__ = ["BusBar"]


@dataclass
class BusBar(Component):
    name: str
    S: Quantity
    I_z: Quantity
    L: Quantity
    U_l: Quantity
    I_b: Quantity
    cos_phi: float
    conductor_material: ConductorMaterial = ConductorMaterial.COPPER
    T_max_cont: Quantity = Q_(100, 'degC')

    I_n: Quantity = field(init=False, default=None)
    I_sc_max: Quantity = field(init=False, default=None)
    I_sc_min: Quantity = field(init=False, default=None)

    Z_dict: dict[str, dict[int, Quantity]] = field(init=False, default_factory=dict)
    z0_r_factor: float = 1.0
    z0_x_factor: float = 1.0

    dU: Quantity = field(init=False, default=None)
    dU_rel: Quantity = field(init=False, default=None)

    circuit_breaker: CircuitBreaker = field(init=False, default=None)

    def __post_init__(self):
        super().__init__(self.name)
        self.I_n = _get_nominal_current(self.I_z)
        self.Z_dict = self._calc_impedance_dict()
        self.dU, self.dU_rel = self.get_voltage_drop(self.U_l, self.I_b, self.cos_phi)

    def get_impedance(
        self,
        T: Quantity
    ) -> dict[int, Quantity]:
        from python_electric.equipment import BusBar
        bus_bar = BusBar(
            length=self.L,
            cross_section_area=self.S,
            conductor_material=self.conductor_material,
            temperature=T,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor
        )
        return {1: bus_bar.Z1, 2: bus_bar.Z2, 0: bus_bar.Z0}

    def _calc_impedance_dict(self) -> dict[str, dict[int, Quantity]]:
        return {
            "T20": self.get_impedance(Q_(20, 'degC')),
            "T_n": self.get_impedance(self.T_max_cont),
            "T150": self.get_impedance(Q_(150, 'degC'))
        }

    def get_voltage_drop(
        self,
        U_l: Quantity,
        I_b: Quantity,
        cos_phi: float,
        volt_ref: VoltReference = VoltReference.PH3_GROUND_TO_LINE
    ) -> tuple[Quantity, Quantity]:
        if volt_ref == VoltReference.PH3_GROUND_TO_LINE:
            U = U_l.to('V') / math.sqrt(3)
        elif volt_ref == VoltReference.PH3_LINE_TO_LINE:
            U = U_l.to('V')
        else:
            raise ValueError("volt_ref must refer to a 3-phase network.")

        R = Q_(self.Z_dict["T_n"][1].m.real, 'ohm')
        X = Q_(self.Z_dict["T_n"][1].m.imag, 'ohm')

        dU = voltage_drop(R, X, I_b, volt_ref, cos_phi)
        dU_rel = dU / U
        return dU, dU_rel.to('pct')

    def connect_circuit_breaker(
        self,
        standard: CircuitBreaker.Standard,
        category: CircuitBreaker.Category,
        I_cu: Quantity,
        I_sc_max: Quantity | None = None,
        I_sc_min: Quantity | None = None,
        k_m: float | None = None,
        E_t: Quantity | None = None,
        t_m: Quantity | None = None
    ) -> None:
        I_sc_max = I_sc_max if I_sc_max is not None else self.I_sc_max
        I_sc_min = I_sc_min if I_sc_min is not None else self.I_sc_min

        cb = CircuitBreaker(
            standard=standard,
            category=category,
            I_b=self.I_b,
            I_n=self.I_n,
            I_z=self.I_z,
            I2t=None,
            I_cu=I_cu,
            E_t=E_t,
            k_m=k_m,
            t_m=t_m
        )

        c1 = cb.check_overload_protection()
        c2 = cb.check_shortcircuit_protection(I_sc_max, I_sc_min)
        if not (c1 and c2):
            raise ValueError("Circuit breaker")

        self.circuit_breaker = cb
        self.I_sc_max = I_sc_max
        self.I_sc_min = I_sc_min
        return None

    @property
    def I_b_tot(self) -> Quantity:
        return self.I_b

    @I_b_tot.setter
    def I_b_tot(self, v: Quantity) -> None:
        self.I_b = v

    @property
    def I_n_tot(self) -> Quantity:
        return self.I_n

    @I_n_tot.setter
    def I_n_tot(self, v: Quantity) -> None:
        self.I_n = v

    @property
    def I_z_tot(self) -> Quantity:
        return self.I_z

    @I_z_tot.setter
    def I_z_tot(self, v: Quantity) -> None:
        self.I_z = v


def _get_nominal_current(I_z: Quantity) -> Quantity:
    """
    Returns the standardized nominal current which is just smaller than the
    given current-carrying capacity `I_z` of the bus bar.
    """
    nom_values = [1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0, 5.0, 6.3, 8.0]
    I_z = I_z.to('A').m

    if 0 < I_z <= 10:
        I_z = round(I_z, 2)
        c = 1
    elif 10 < I_z <= 100:
        I_z = round(I_z / 10, 2)
        c = 10
    elif 100 < I_z <= 1000:
        I_z = round(I_z / 100, 2)
        c = 100
    elif 1000 < I_z <= 8000:
        I_z = round(I_z / 1000, 2)
        c = 1000
    else:
        raise ValueError("Ampacity cannot exceed 8 kA.")

    nom_values.append(I_z)
    nom_values.sort()
    i = nom_values.index(I_z)
    if i == 0:
        if c > 1:
            I_n = nom_values[-1] * c * 0.1
        else:
            raise ValueError("Nominal current out of range.")
    else:
        I_n = nom_values[i - 1] * c
    return Q_(I_n, 'A')
