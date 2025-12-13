import math
from dataclasses import dataclass, field

from ... import Quantity
from ... general import VoltReference
from ...calc import get_3ph_line_current, voltage_drop
from ..cable import get_nominal_current

__all__ = [
    "Transformer",
]

Q_ = Quantity


@dataclass
class Transformer:
    """
    Represents a transformer in a low-voltage three-phase system.

    Calculates the impedance of the transformer on instantiation of the class.
    To calculate the voltage drop across the transformer, use method
    `voltage_drop()`.

    Parameters
    ----------
    S_nom: Quantity
        Rated apparent power of the transformer.
    U_line_pri: Quantity
        Rated primary line-to-line voltage of the transformer.
    U_line: Quantity
        Rated secondary line-to-line voltage of the transformer.
    u_cc: Quantity
        Relative short-circuit voltage of the transformer.
    P_Cu: Quantity
        Rated copper loss of the transformer.
    P_load: Quantity
        Rated three-phase active power of the connected load.
    cos_phi: Quantity
        Rated power factor of the connected load.
    name: str, optional
        To identify the unique object.

    Attributes
    ----------
    I_nom_pri: Quantity
        Calculated rated primary current of the transformer.
    I_nom_sec: Quantity
        Calculated rated secondary current of the transformer.
    I_nom: Quantity
        Standardized nominal current of the transformer which is just greater
        than the calculated rated secondary current. This current is used to
        size the downstream cable or bus bars.
    I_load: Quantity
        Rated secondary load current delivered by the transformer.
    Z: dict[int, Quantity]
        Zero-sequence (key: 0), positive sequence (key: 1), and negative
        sequence (key: 2) impedances of the transformer referred to the
        secondary side.
    U_drop: Quantity
        Absolute drop in voltage across the transformer in volts. Available
        after method voltage_drop() has been called.
    U_drop_pct: Quantity
        Relative drop in voltage in percent. Available after method
        voltage_drop() has been called.
    """
    S_nom: Quantity
    U_line_pri: Quantity
    U_line: Quantity
    u_cc: Quantity
    P_Cu: Quantity
    P_load: Quantity
    cos_phi: Quantity
    U_drop: Quantity = field(init=False)
    U_drop_pct: Quantity = field(init=False)
    name: str = field(default="transformer")

    def __post_init__(self):
        self.Z: dict[int, Quantity] = self._calculate_sec_impedance()
        self.I_nom_pri = self._calculate_pri_nom_current()
        self.I_nom_sec, self.I_nom = self._calculate_sec_nom_current()
        self.I_load = self._calculate_sec_load_current()

    def _calculate_sec_impedance(self) -> dict[int, Quantity]:
        from python_electric.equipment import Transformer
        tr = Transformer(
            nominal_voltage=self.U_line,
            nominal_power=self.S_nom,
            percent_short_circuit_voltage=self.u_cc,
            copper_loss=self.P_Cu,
            voltage_factor=1.1
        )
        return {1: tr.Z1, 2: tr.Z2, 0: tr.Z0}

    def _calculate_pri_nom_current(self) -> Quantity:
        return self.S_nom / (math.sqrt(3) * self.U_line_pri)

    def _calculate_sec_nom_current(self) -> tuple[Quantity, Quantity]:
        I_nom_sec = self.S_nom / (math.sqrt(3) * self.U_line)
        I_nom = get_nominal_current(I_nom_sec.to('A').m)
        return I_nom_sec, Q_(I_nom, 'A')

    def _calculate_sec_load_current(self) -> Quantity:
        return Q_(get_3ph_line_current(self.P_load, self.U_line, self.cos_phi), 'A')

    def voltage_drop(
        self,
        U_line: Quantity | None = None,
        I_load: Quantity | None = None,
        cos_phi: float | None = None,
        volt_ref: VoltReference = VoltReference.PH3_LINE_TO_LINE
    ) -> None:
        """
        Calculates the absolute and relative voltage drop across the
        transformer. Results are stored in attributes `U_drop` and `U_drop_pct`.

        Parameters
        ----------
        U_line: optional
            Three-phase line-to-line voltage at the entry.
        I_load: optional
            Load current.
        cos_phi: optional
            Power factor.
        volt_ref: {VoltRef.PH3_LINE_TO_LINE, VoltRef.PH3_GROUND_TO_LINE}
            Indicates the type of voltage drop in a three-phase system, and to
            which voltage the relative voltage drop is referred.
            If VoltRef.PH3_LINE_TO_LINE, the voltage drop is the drop in voltage
            measured between two line conductors, and the relative voltage is
            referred to the line-to-line voltage of the three-phase system.
            If VoltRef.PH3_GROUND_TO_LINE, the voltage drop is the drop in
            voltage measured across a single line conductor, and the relative
            voltage is referred to the ground-to-line voltage of the three-phase
            system.

        Returns
        -------
        None
        """
        if U_line is None:
            U_line = self.U_line
        if I_load is None:
            I_load = self.I_load
        if cos_phi is None:
            cos_phi = self.cos_phi
        R = self.Z[1].real
        X = self.Z[1].imag
        U_drop = voltage_drop(R, X, I_load, volt_ref, cos_phi)
        if volt_ref == VoltReference.PH3_GROUND_TO_LINE:
            U = U_line / math.sqrt(3)
        elif volt_ref == VoltReference.PH3_LINE_TO_LINE:
            U = U_line
        else:
            raise ValueError(
                "Parameter `volt_ref` must refer to a three-phase system."
                "Either PH3_GROUND_TO_LINE or PH3_LINE_TO_LINE."
            )
        U_drop_pct = 100 * U_drop / U.to('V').m

        self.U_drop = Q_(U_drop, "V")
        self.U_drop_pct = Q_(U_drop_pct, "pct")
        return None
