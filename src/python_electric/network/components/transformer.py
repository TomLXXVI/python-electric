import math
from dataclasses import dataclass, field
from enum import StrEnum

from python_electric import Quantity, Q_, VoltReference
from python_electric import sizing
from python_electric.calc import voltage_drop

from python_electric.network.graph import Component

__all__ = ["Transformer"]


@dataclass
class Transformer(Component):
    """
    Represents a three-phase transformer in a low-voltage network.

    Parameters
    ----------
    name: str, optional
        Uniquely identifies the transformer in the network.
    S_n: Quantity
        Rated apparent power of the transformer.
    U_lp: Quantity
        Rated primary line-to-line voltage of the transformer.
    U_ls: Quantity
        Rated secondary line-to-line voltage of the transformer.
    u_cc: Quantity
        Relative short-circuit voltage of the transformer.
    P_Cu: Quantity
        Rated copper loss of the transformer.
    I_b: Quantity
        Load current.
    cos_phi: Quantity
        Power factor of the connected load.
    pri_conn: WindingConn
        Primary windings connection. See enum WindingConn.
    sec_conn: WindingConn
        Secondary windings connection.
    Zn_pri: Quantity, optional
        Neutral-to-earth impedance of primary windings, if present.
    Zn_sec: Quantity, optional
        Neutral-to-earth impedance of secondary windings, if present.
    z0_r_factor: float, default 3.0
        Scaling factor applied to the positive-sequence resistance of the
        transformer to determine its zero-sequence resistance.
    z0_x_factor: float, default 3.0
        Scaling factor applied to the negative-sequence reactance of the
        transformer to determine its zero-sequence reactance.

    Attributes
    ----------
    I_p: Quantity
        Calculated rated primary current of the transformer.
    I_s: Quantity
        Calculated rated secondary current of the transformer.
    I_sn: Quantity
        Standardized nominal current of the transformer which is just greater
        than the calculated rated secondary current. This current is used to
        size the downstream cable or bus bars.
    I_b: Quantity
        Rated secondary load current delivered by the transformer.
    Z_dict: dict[int, Quantity]
        Zero-sequence (key: 0), positive sequence (key: 1), and negative
        sequence (key: 2) impedances of the transformer referred to the
        secondary side.
    dU: Quantity
        Absolute voltage drop across the transformer. By default, the voltage
        drop is measured with respect to ground.
    dU_rel: Quantity
        Relative voltage drop.

    Methods
    -------
    get_voltage_drop:
        Returns the voltage drop across the transformer.
    get_impedance:
        Returns the impedance of the transformer referred to its secondary side.
    """
    class WindingConn(StrEnum):
        Y = "Y"
        YN = "YN"
        D = "D"

    name: str
    S_n: Quantity
    U_lp: Quantity
    U_ls: Quantity
    u_cc: Quantity
    P_Cu: Quantity
    I_b: Quantity
    cos_phi: float

    pri_conn: WindingConn | None = None
    sec_conn: WindingConn | None = None
    Zn_pri: Quantity = Q_(0, 'ohm')
    Zn_sec: Quantity = Q_(0, 'ohm')

    I_p: Quantity = field(init=False)
    I_s: Quantity = field(init=False)
    I_sn: Quantity = field(init=False)

    dU: Quantity = field(init=False)
    dU_rel: Quantity = field(init=False)

    Z_dict: dict[int, Quantity] = field(init=False, default_factory=dict)
    z0_r_factor: float = 1.0
    z0_x_factor: float = 1.0

    def __post_init__(self):
        super().__init__(self.name)
        self.I_p = self._get_rated_pri_current()
        self.I_s = self._get_rated_sec_current()
        self.I_sn = self._get_stand_sec_current()
        self.Z_dict = self.get_impedance()
        self.dU, self.dU_rel = self.get_voltage_drop(self.U_ls, self.I_b, self.cos_phi)

    def _get_rated_pri_current(self) -> Quantity:
        I_np = self.S_n / (math.sqrt(3) * self.U_lp)
        return I_np.to('A')

    def _get_rated_sec_current(self) -> Quantity:
        I_ns = self.S_n / (math.sqrt(3) * self.U_ls)
        return I_ns.to('A')

    def _get_stand_sec_current(self) -> Quantity:
        I_sn = sizing.get_nominal_current(self.I_s.to('A').m)
        return Q_(I_sn, 'A')

    def get_voltage_drop(
        self,
        U_l: Quantity,
        I_b: Quantity,
        cos_phi: float,
        volt_ref: VoltReference = VoltReference.PH3_GROUND_TO_LINE
    ) -> tuple[Quantity, Quantity]:
        """
        Returns the absolute and relative voltage drop across the transformer.

        Parameters
        ----------
        U_l: Quantity, optional
            Line-to-line voltage.
        I_b: Quantity, optional
            Load current.
        cos_phi: float, optional
            Power factor.
        volt_ref: VoltReference, {VoltRef.PH3_LINE_TO_LINE, VoltRef.PH3_GROUND_TO_LINE (default)}
            Reference of the voltage measurement: either between two lines
            (line-to-line), or between ground (neutral) and a line
            (ground-to-line).

        Returns
        -------
        tuple[Quantity, Quantity]
            First element is the absolute voltage drop, second element the
            relative voltage drop.
        """
        if volt_ref == VoltReference.PH3_GROUND_TO_LINE:
            U = U_l.to('V') / math.sqrt(3)
        elif volt_ref == VoltReference.PH3_LINE_TO_LINE:
            U = U_l.to('V')
        else:
            raise ValueError("volt_ref must refer to a 3-phase network.")

        R = Q_(self.Z_dict[1].m.real, 'ohm')
        X = Q_(self.Z_dict[1].m.imag, 'ohm')

        dU = voltage_drop(R, X, I_b, volt_ref, cos_phi)
        dU_rel = dU / U
        return dU, dU_rel.to('pct')

    def get_impedance(self, volt_factor: float = 1.1) -> dict[int, Quantity]:
        """
        Returns the positive, negative, and zero-sequence impedance of the
        transformer referred to the secondary side of the transformer.

        Parameters
        ----------
        volt_factor: float, default 1.1
            Voltage factor used to calculate the impedance of the transformer.

        Returns
        -------
        dict[int, Quantity]
            A dictionary of which the keys 1, 2, 0 map to the positive,
            negative, and zero-sequence impedance respectively.
        """
        from python_electric.equipment import Transformer
        t = Transformer(
            nominal_voltage=self.U_ls,
            nominal_power=self.S_n,
            percent_short_circuit_voltage=self.u_cc,
            copper_loss=self.P_Cu,
            voltage_factor=volt_factor,
            z0_r_factor=self.z0_r_factor,
            z0_x_factor=self.z0_x_factor,
            name=self.name
        )
        return {1: t.Z1, 2: t.Z2, 0: t.Z0}
