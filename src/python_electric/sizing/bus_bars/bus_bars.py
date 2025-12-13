import math
from dataclasses import dataclass, field

from ... import Quantity
from ...materials import ConductorMaterials
from ...protection import CircuitBreaker
from ...general import VoltReference
from ...calc import get_3ph_line_current, voltage_drop

__all__ = [
    "BusBars",
]

Q_ = Quantity


def _get_nominal_current(I_z: Quantity) -> Quantity:
    """
    Returns the standardized nominal current which is just smaller than the
    given current-carrying capacity.
    """
    I_z = I_z.to('A').m
    nom_sizes = [1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0, 5.0, 6.3, 8.0]
    if 0 < I_z <= 10:
        I_z = round(I_z, 2)
        c = 1
    elif 10 < I_z <= 100:
        I_z = round(I_z / 10, 2)
        c = 10
    elif 100 < I_z <= 1000:
        I_z = round(I_z / 100, 2)
        c = 100
    elif 1000 < I_z <= 8_000:
        I_z = round(I_z / 1000, 2)
        c = 1000
    else:
        raise ValueError("Ampacity cannot exceed 8 kA.")
    nom_sizes.append(I_z)
    nom_sizes.sort()
    i = nom_sizes.index(I_z)
    if i == 0:
        if c > 1:
            I_nom = nom_sizes[-1] * c * 0.1
        else:
            raise ValueError("Nominal current out of range.")
    else:
        I_nom = nom_sizes[i - 1] * c
    return Q_(I_nom, 'A')


@dataclass
class BusBars:
    """
    Represents bus bars in a low-voltage three-phase system.

    Calculates the impedance of the bus bars on instantiation of the class.

    Parameters
    ----------
    S: Quantity
        Cross-sectional area of a single bus bar selected based on the maximum
        continuous secondary current of the transformer (not the possible
        smaller rated current of the connected load).
    I_z: Quantity
        Current-carrying capacity of a single bus bar that goes with the
        selected cross-sectional area.
    L: Quantity
        Length of the bus bars.
    U_line: Quantity
        Three-phase line-to-line voltage.
    I_load: Quantity | None = None
        Rated load current through the bus bars. If left to None, parameter P
        must be set.
    P_load: Quantity | None = None
        Three-phase active power that flows through the bus bars. If left to
        None, parameter I_b must be set.
    cos_phi: float = 1.0
        Power factor of the load, i.e. the cosine of the phase shift between
        voltage and current.
    name: str, optional
        To identify the unique object.

    Attributes
    ----------
    I_nom: Quantity
        Standardized nominal current which is just smaller than the
        current-carrying ampacity of the bus bar.
    Z: dict[str, dict[int, Quantity]]
        Dictionary that contains the zero-sequence (key 0), positive sequence
        (key 1), and negative sequence (key 2) impedances of the cable at three
        different temperatures:
        -   at 20°C (key "T_20") used for calculating the maximum short-circuit
            current.
        -   at the maximum continuous temperature (key "T_nom").
        -   at 150°C (key "T_150") used for calculating the minimum
            short-circuit current.
        E.g. to get Z1 at 150 °C, use Z["T_150"][1].
    U_drop: Quantity
        Absolute drop in voltage across the cable in volts. Available after
        method `voltage_drop()` has been called.
    U_drop_pct: Quantity
        Relative drop in voltage in percent. Available after method
        `voltage_drop()` has been called.
    circuit_breaker: CircuitBreaker
        Circuit breaker connected to the cable. Available after method
        `connect_circuit_breaker()` has been called.

    Methods
    -------
    voltage_drop:
        Calculates the voltage drop across the bus bars.
    connect_circuit_breaker:
        Connects a circuit breaker to the bus bars.
    """
    S: Quantity
    I_z: Quantity
    L: Quantity
    U_line: Quantity
    I_load: Quantity | None = None
    P_load: Quantity | None = None
    cos_phi: float = 1.0
    T_max_cont: Quantity = Q_(100, 'degC')
    name: str = "bus_bars"

    U_drop: Quantity = field(init=False)
    U_drop_pct: Quantity = field(init=False)
    circuit_breaker: CircuitBreaker = field(init=False, default=None)

    def __post_init__(self):
        if self.I_load is None:
            if self.P_load is not None:
                self.I_load = Q_(get_3ph_line_current(self.P_load, self.U_line, self.cos_phi), 'A')
            else:
                raise ValueError(
                    "Parameters I_load and P_load cannot both be None."
                )
        
        self.I_nom = _get_nominal_current(self.I_z)
        
        self.Z: dict[str, dict[int, Quantity]] = {
            "T_20": self._calc_impedance(Q_(20, 'degC')),
            "T_nom": self._calc_impedance(self.T_max_cont),
            "T_150": self._calc_impedance(Q_(150, 'degC'))
        }

    def _calc_impedance(self, T: Quantity) -> dict[int, Quantity]:
        from python_electric.equipment import BusBar
        busbar = BusBar(
            length=self.L,
            cross_section_area=self.S,
            conductor_material=ConductorMaterials.COPPER,
            temperature=T
        )
        return {1: busbar.Z1, 2: busbar.Z2, 0: busbar.Z0}

    def voltage_drop(
        self,
        U_line: Quantity | None = None,
        I_load: Quantity | None = None,
        cos_phi: float | None = None,
        volt_ref: VoltReference = VoltReference.PH3_LINE_TO_LINE
    ) -> None:
        """
        Calculates the absolute and relative voltage drop across the
        bus bars. Results are stored in attributes `U_drop` and `U_drop_pct`.

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
        R = self.Z["T_nom"][1].real
        X = self.Z["T_nom"][1].imag
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
        self.U_drop=Q_(U_drop, "V")
        self.U_drop_pct=Q_(U_drop_pct, "pct")
        return None

    def connect_circuit_breaker(
        self,
        standard: CircuitBreaker.Standard,
        category: CircuitBreaker.Category,
        ultim_break_capacity: Quantity,
        I_sc_max: Quantity,
        I_sc_min: Quantity,
        k_magn_trip: float | None = None,
        let_through_energy: Quantity | None = None,
        t_m_lim: Quantity | None = None
    ) -> None:
        """
        Connects a circuit breaker with the bus bars.

        Parameters
        ----------
        standard: CircuitBreaker.Standard
            Either the residential standard (IEC 60898-1) or the industrial
            standard (IEC 60947-2), which specifies the minimal performance
            requirements of the circuit breaker. See enum
            CircuitBreaker.Standard.
        category: CircuitBreaker.Category
            Specifies the category of the circuit breaker depending on its
            magnetic trip threshold. Either category B, C, D, or AJUSTABLE. See
            enum CircuitBreaker.Category. Note that category ADJUSTABLE also
            demands that the circuit breaker is of the industrial type.
        ultim_break_capacity: Quantity
            Ultimate breaking capacity of the circuit breaker.
        I_sc_max: Quantity
            Calculated maximum short-circuit current at the location of the
            circuit breaker (most often caused by a three-phase fault).
        I_sc_min: Quantity
            Calculated minimum short-circuit current at the end of the cable
            (usually caused by a single line-to-ground fault).
        k_magn_trip: float, optional
            Multiplication factor that determines the rated magnetic trip
            current as a multiple of the thermal current setting if the circuit
            breaker is of the industrial type and adjustable.
        let_through_energy: Quantity, optional
            Thermal energy in A².s let through by the circuit breaker at the
            calculated maximum short-circuit current.
        t_m_lim: Quantity, optional
            Upper limit of instantaneous magnetic tripping time with regard to
            short-circuits. By default, this time limit is set to 100 ms
            according to the residential standard IEC 60898-1.
        """
        cb = CircuitBreaker(
            standard=standard,
            category=category,
            load_current=self.I_load,
            nom_current=self.I_nom,
            ampacity=self.I_z,
            joule_integral=None,
            ultim_break_capacity=ultim_break_capacity,
            let_through_energy=let_through_energy,
            k_magn_trip=k_magn_trip,
            t_m_lim=t_m_lim
        )
        overload_ok = cb.check_overload_protection()
        short_circuit_ok = cb.check_shortcircuit_protection(I_sc_max, I_sc_min)
        if not (overload_ok and short_circuit_ok):
            raise ValueError("circuit-breaker protection fails.")
        self.circuit_breaker = cb
        return None
