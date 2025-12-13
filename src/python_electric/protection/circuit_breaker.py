from __future__ import annotations

from enum import StrEnum
import warnings

from .. import Quantity


__all__ = [
    "CircuitBreaker",
    "check_current_based_selectivity",
    "ProtectionWarning",
    "SelectivityWarning"
]

Q_ = Quantity


class ProtectionWarning(Warning):
    pass


class SelectivityWarning(Warning):
    pass


class CircuitBreaker:
    """
    Represents a circuit-breaker, according to either the residential standard
    (IEC 60898-1) or the industrial standard (IEC 60947-2).

    Attributes
    ----------
    I_load: Quantity
        Rated load current through the protected cable.
    I_nom: Quantity
        Standardized nominal current which is just greater than or equal to
        the rated load current.
    I_r: Quantity
        Thermal current setting of the circuit breaker. In case of an adjustable
        circuit breaker, the thermal current setting of the circuit breaker is
        set equal to the rated load current. In all other cases, it is set equal
        to the nominal current.
    I_z: Quantity
        Current-carrying capacity of the protected cable.
    I_cu: Quantity
        Ultimate breaking capacity of the circuit breaker.
    I_nf: Quantity
        Conventional non-tripping current with regard to overload.
    I_f: Quantity
        Conventional tripping current with regard to overload.
    t_conv: Quantity
        Conventional tripping time with regard to overload.
    I_m_min: Quantity
        Lower limit of magnetic tripping current with regard to short-circuits.
    I_m_max: Quantity
        Upper limit of magnetic tripping current with regard to short-circuits.
    E_through: Quantity
        Thermal energy (I²t) let through by the circuit breaker at the
        calculated maximum short-circuit current during the time that is
        needed to interupt the current.
    """

    class Standard(StrEnum):
        RESIDENTIAL = "IEC 60898-1"
        INDUSTRIAL = "IEC 60947-2"

    class Category(StrEnum):
        B = "B"
        C = "C"
        D = "D"
        ADJUSTABLE = "adjustable"

    def __init__(
        self,
        standard: CircuitBreaker.Standard,
        category: CircuitBreaker.Category,
        load_current: Quantity,
        nom_current: Quantity,
        ampacity: Quantity,
        joule_integral: Quantity | None,
        ultim_break_capacity: Quantity,
        let_through_energy: Quantity | None = None,
        k_magn_trip: float | None = None,
        t_m_lim: Quantity | None = None
    ) -> None:
        """
        Creates a `CircuitBreaker` object.

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
        load_current: Quantity
            Rated load current through the protected cable.
        nom_current: Quantity
            Standardized nominal current which is just greater than or equal to
            the rated load current.
        ampacity: Quantity
            Current-carrying capacity of the protected cable.
        joule_integral: Quantity | None
            Joule-integral of the cable, i.e. the maximum thermal energy,
            usually expressed in square amperes x seconds, the conductors are
            capable to store without thermal degradation of the insulation
            assuming adiabatic heating (i.e. without heat transfer to the
            ambient). Set to None, if the joule integral is not known.
        ultim_break_capacity: Quantity
            Ultimate breaking capacity of the circuit breaker.
        let_through_energy: Quantity, optional
            Thermal energy (I²t) let through by the circuit breaker at the
            calculated maximum short-circuit current during the time that is
            needed to interupt the current.
        k_magn_trip: float, optional
            Multiplication factor that determines the rated magnetic trip 
            current as a multiple of the thermal current setting if the circuit 
            breaker is of the industrial type and adjustable.
        t_m_lim: Quantity, optional
            Upper limit of instantaneous magnetic tripping time with regard to
            short-circuits. By default, this time limit is set to 100 ms
            according to the residential standard IEC 60898-1.
        """
        self.standard = standard
        self.category = category
        self.I_load = load_current.to('A')
        self.I_nom = nom_current.to('A')

        if (self.standard == CircuitBreaker.Standard.INDUSTRIAL
            and self.category == CircuitBreaker.Category.ADJUSTABLE):
            self.I_r = self.I_load
        else:
            self.I_r = self.I_nom

        self.I_z = ampacity.to('A')

        self.joule_integral = None
        if joule_integral is not None:
            self.joule_integral = joule_integral.to('A ** 2 * s')

        self.I_cu = ultim_break_capacity.to('A')

        self.E_through = None
        if let_through_energy is not None:
            self.E_through = let_through_energy.to('A ** 2 * s')

        self.k_magn_trip = k_magn_trip

        self.t_m_lim = t_m_lim
        if self.t_m_lim is None:
            self.t_m_lim = self._get_magn_trip_time_limit()

        # Characteristics of the circuit breaker according to standards:
        self.I_nf = self._get_conv_non_trip_current()
        self.I_f = self._get_conv_trip_current()
        self.t_conv = self._get_conv_trip_time()
        self.I_m_min = self._get_min_magn_trip_threshold()
        self.I_m_max = self._get_max_magn_trip_threshold()

    def _get_conv_non_trip_current(self) -> Quantity:
        if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
            return 1.13 * self.I_r
        else:
            return 1.05 * self.I_r

    def _get_conv_trip_current(self) -> Quantity:
        if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
            return 1.45 * self.I_r
        else:
            return 1.30 * self.I_r

    def _get_conv_trip_time(self) -> Quantity:
        if self.I_r.to('A').m <= 63.0:
            return Q_(1.0, 'h')
        else:
            return Q_(2.0, 'h')

    def _get_min_magn_trip_threshold(self) -> Quantity:
        match self.category:
            case CircuitBreaker.Category.B:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 3.0 * self.I_r
                else:
                    return 4.0 * self.I_r * (1.0 - 0.2)
            case CircuitBreaker.Category.C:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 5.0 * self.I_r
                else:
                    return 8.0 * self.I_r * (1.0 - 0.2)
            case CircuitBreaker.Category.D:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 10 * self.I_r
                else:
                    return 12 * self.I_r * (1.0 - 0.2)
            case CircuitBreaker.Category.ADJUSTABLE:
                if self.standard == CircuitBreaker.Standard.INDUSTRIAL:
                    if self.k_magn_trip is not None:
                        return self.k_magn_trip * self.I_r * (1.0 - 0.2)
                    else:
                        raise AttributeError(
                            "Parameter `K_magn_trip` cannot be None."
                        )
                else:
                    raise AttributeError(
                        "Adjustable circuit breaker cannot be residential."
                    )

    def _get_max_magn_trip_threshold(self) -> Quantity:
        match self.category:
            case CircuitBreaker.Category.B:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 5.0 * self.I_r
                else:
                    return 4.0 * self.I_r * (1.0 + 0.2)
            case CircuitBreaker.Category.C:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 10.0 * self.I_r
                else:
                    return 8.0 * self.I_r * (1.0 + 0.2)
            case CircuitBreaker.Category.D:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 20.0 * self.I_r
                else:
                    return 12.0 * self.I_r * (1.0 + 0.2)
            case CircuitBreaker.Category.ADJUSTABLE:
                if self.standard == CircuitBreaker.Standard.INDUSTRIAL:
                    if self.k_magn_trip is not None:
                        return self.k_magn_trip * self.I_r * (1.0 + 0.2)
                    else:
                        raise AttributeError(
                            "Parameter `K_magn_trip` cannot be None."
                        )
                else:
                    raise AttributeError(
                        "An adjustable circuit breaker cannot be residential."
                    )

    @staticmethod
    def _get_magn_trip_time_limit() -> Quantity:
        return Q_(0.1, 's')

    @property
    def has_adjustable_delay(self) -> bool:
        """
        Return True if this breaker can have an intentional short-time delay.
        """
        return (
                self.standard == CircuitBreaker.Standard.INDUSTRIAL
                and self.category == CircuitBreaker.Category.ADJUSTABLE
        )

    @property
    def is_industrial(self) -> bool:
        return self.standard == CircuitBreaker.Standard.INDUSTRIAL

    def check_overload_protection(self) -> bool:
        """
        Checks if all overload protection requirements for the cable are
        fulfilled by the selected circuit breaker.

        Warnings
        --------
        ProtectionWarning

        Returns
        -------
        bool
        """
        I_b = self.I_load
        I_n = self.I_r
        I_z = self.I_z
        I_nf = self.I_nf
        I_f = self.I_f
        if I_b <= I_n <= I_z:
            pass
        else:
            warnings.warn(
                f"I_nom {I_n.to('A'):~P.1f} is not in the range "
                f"[{I_b.to('A'):~P.1f}, {I_z.to('A'):~P.1f}].",
                category=ProtectionWarning
            )
            return False
        if I_nf <= 1.15 * I_z:
            pass
        else:
            warnings.warn(
                f"I_nf {I_nf.to('A'):~P.2f} > {(1.15 * I_z).to('A'):~P.2f}.",
                category=ProtectionWarning
            )
            return False
        if I_f <= 1.45 * I_z:
            pass
        else:
            warnings.warn(
                f"I_f {I_f.to('A'):~P.2f} > {(1.45 * I_z).to('A'):~P.2f}.",
                category=ProtectionWarning
            )
            return False
        return True

    def check_shortcircuit_protection(
        self,
        I_sc_max: Quantity,
        I_sc_min: Quantity
    ) -> bool:
        """
        Checks if the short-circuit protection requirements for the cable are
        fulfilled by the selected circuit breaker.

        Parameters
        ----------
        I_sc_max: Quantity
            Calculated maximum short-circuit current at the location of the
            circuit breaker (most often caused by a three-phase fault).
        I_sc_min: Quantity
            Calculated minimum short-circuit current at the location of the next
            downstream circuit breaker or electrical consumer (usually caused by
            a single line-to-ground fault).

        Warnings
        --------
        ProtectionWarning

        Returns
        -------
        bool
        """
        I_sc_max = I_sc_max.to('A')
        I_sc_min = I_sc_min.to('A')
        self.I_cu.ito('A')
        self.I_m_min.ito('A')

        if I_sc_max <= self.I_cu:
            pass
        else:
            warnings.warn(
                "I_sc_max > ultimate breaking capacity",
                category=ProtectionWarning
            )
            return False

        if I_sc_min >= self.I_m_min:
            pass
        else:
            warnings.warn(
                f"I_sc_min < minimum limit of magnetic trip current.",
                category=ProtectionWarning
            )
            return False

        if self.joule_integral is not None:
            self.joule_integral.ito('A ** 2 * s')

            if self.E_through is not None:
                E_through = self.E_through.to('A ** 2 * s')

                if self.joule_integral >= E_through:
                    pass
                else:
                    warnings.warn(
                        "Let-through energy exceeds Joule-integral of the "
                        "cable.",
                        category=ProtectionWarning
                    )
                    return False

            t_allow = self.joule_integral / (I_sc_min ** 2)
            if self.t_m_lim.to('s') <= t_allow.to('s'):
                pass
            else:
                warnings.warn(
                    "Magnetic interruption time exceeds the allowable fault "
                    "duration.",
                    category=ProtectionWarning
                )
                return False
            return True
        return True
    
    def __str__(self):
        s =  f"I_n     = {self.I_nom.to('A'):~P.1f}\n"
        s += f"I_r     = {self.I_r.to('A'):~P.1f}\n"
        s += f"I_nf    = {self.I_nf.to('A'):~P.1f}\n"
        s += f"I_f     = {self.I_f.to('A'):~P.1f}\n"
        s += f"t_conv  = {self.t_conv.to('h'):~P.1f}\n"
        s += f"I_m_min = {self.I_m_min.to('A'):~P.1f}\n"
        s += f"I_m_max = {self.I_m_max.to('A'):~P.1f}\n"
        s += f"t_m_lim = {self.t_m_lim.to('s'):~P.1f}"
        return s


def check_current_based_selectivity(
    cb_down: CircuitBreaker,
    cb_up: CircuitBreaker,
) -> bool:
    """
    Checks if current-based selectivity between a downstream circuit breaker and
    its nearest upstream circuit breaker is total.

    Parameters
    ----------
    cb_down:
        Downstream circuit breaker.
    cb_up:
        Nearest upstream circuit breaker.

    Notes
    -----
    Current-based selectivity is considered to be total if the following two
    conditions are both met:
    1.  The conventional tripping current `I_f` of the downstream circuit
        breaker is smaller than the conventional non-tripping current `I_nf` of
        the upstream circuit breaker.
    2.  The maximum limit of the magnetic tripping current `I_m_max` of the
        downstream circuit breaker is smaller than the minimum limit `I_m_min`
        of the magnetic tripping current of the upstream circuit breaker.

    Warnings
    --------
    SelectivityWarning

    Returns
    -------
    bool
    """
    I_f_down = cb_down.I_f.to('A')
    I_nf_up = cb_up.I_nf.to('A')
    I_m_down = cb_down.I_m_max.to('A')
    I_m_up = cb_up.I_m_min.to('A')
    if I_f_down < I_nf_up:
        pass
    else:
        warnings.warn(
            "Downstream I_f > upstream I_nf.",
            category=SelectivityWarning
        )
        return False
    if I_m_down < I_m_up:
        pass
    else:
        warnings.warn(
            "Downstream I_m_max > upstream I_m_min.",
            category=SelectivityWarning
        )
        return False
    return True
