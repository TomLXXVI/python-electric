from __future__ import annotations

from enum import StrEnum
import warnings

from .. import Quantity, Q_


__all__ = [
    "CircuitBreaker",
    "check_current_based_selectivity",
    "ProtectionWarning",
    "SelectivityWarning"
]


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
    I_b: Quantity
        Rated load current through the protected cable.
    I_n: Quantity
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
    E_t: Quantity
        Thermal energy (I²t) let through by the circuit breaker at the
        calculated maximum short-circuit current during the time needed to
        interupt the current.
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
        I_b: Quantity,
        I_n: Quantity,
        I_z: Quantity,
        I2t: Quantity | None,
        I_cu: Quantity,
        E_t: Quantity | None = None,
        k_m: float | None = None,
        t_m: Quantity | None = None
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
        I_b: Quantity
            Rated load current through the protected cable.
        I_n: Quantity
            Standardized nominal current which is just greater than or equal to
            the rated load current.
        I_z: Quantity
            Current-carrying capacity of the protected cable.
        I2t: Quantity | None
            Joule-integral of the cable, i.e. the maximum thermal energy,
            usually expressed in square amperes x seconds, the conductors are
            capable to store without thermal degradation of the insulation
            assuming adiabatic heating (i.e. without heat transfer to the
            ambient). Set to None, if the joule integral is not known.
        I_cu: Quantity
            Ultimate breaking capacity of the circuit breaker.
        E_t: Quantity, optional
            Thermal energy (I²t) let through by the circuit breaker at the
            calculated maximum short-circuit current during interruption time.
        k_m: float, optional
            Multiplication factor that determines the rated short-circuit trip
            current I_m as a multiple of the thermal current setting I_r if the
            circuit breaker is of the industrial type and adjustable.
        t_m: Quantity, optional
            Upper limit for the instantaneous magnetic tripping time with regard
            to short-circuits. By default, this time limit is set to 100 ms
            according to the residential standard IEC 60898-1.
        """
        self.standard = standard
        self.category = category
        self.I_b = I_b.to('A')
        self.I_n = I_n.to('A')

        if (self.standard == CircuitBreaker.Standard.INDUSTRIAL
            and self.category == CircuitBreaker.Category.ADJUSTABLE):
            self.I_r = self.I_b
        else:
            self.I_r = self.I_n

        self.I_z = I_z.to('A')

        self.I2t = None
        if I2t is not None:
            self.I2t = I2t.to('A ** 2 * s')

        self.I_cu = I_cu.to('A')

        self.E_t = None
        if E_t is not None:
            self.E_t = E_t.to('A ** 2 * s')

        self.k_m = k_m

        self.t_m = t_m
        if self.t_m is None:
            self.t_m = self._get_magn_trip_time_limit()

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
                    return 4.0 * self.I_r * 0.8
            case CircuitBreaker.Category.C:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 5.0 * self.I_r
                else:
                    return 8.0 * self.I_r * 0.8
            case CircuitBreaker.Category.D:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 10 * self.I_r
                else:
                    return 12 * self.I_r * 0.8
            case CircuitBreaker.Category.ADJUSTABLE:
                if self.standard == CircuitBreaker.Standard.INDUSTRIAL:
                    if self.k_m is not None:
                        return self.k_m * self.I_r * 0.8
                    else:
                        raise AttributeError(
                            "Parameter `k_m` cannot be None."
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
                    return 4.0 * self.I_r * 1.2
            case CircuitBreaker.Category.C:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 10.0 * self.I_r
                else:
                    return 8.0 * self.I_r * 1.2
            case CircuitBreaker.Category.D:
                if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
                    return 20.0 * self.I_r
                else:
                    return 12.0 * self.I_r * 1.2
            case CircuitBreaker.Category.ADJUSTABLE:
                if self.standard == CircuitBreaker.Standard.INDUSTRIAL:
                    if self.k_m is not None:
                        return self.k_m * self.I_r * 1.2
                    else:
                        raise AttributeError(
                            "Parameter `k_m` cannot be None."
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
        I_b = self.I_b
        I_n = self.I_r
        I_z = self.I_z
        I_nf = self.I_nf
        I_f = self.I_f

        # check 1
        if I_b <= I_n <= I_z:
            pass
        else:
            warnings.warn(
                f"I_n {I_n.to('A'):~P.1f} is not in the interval "
                f"[{I_b.to('A'):~P.1f}, {I_z.to('A'):~P.1f}].",
                category=ProtectionWarning
            )
            return False

        # check 2
        if I_nf <= 1.15 * I_z:
            pass
        else:
            warnings.warn(
                f"I_nf {I_nf.to('A'):~P.2f} > {(1.15 * I_z).to('A'):~P.2f}.",
                category=ProtectionWarning
            )
            return False

        # check 3
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

        # Check the ultimate breaking capacity of circuit breaker.
        if I_sc_max <= self.I_cu:
            pass
        else:
            warnings.warn(
                f"I_sc_max ({I_sc_max:~P.0f}) > I_cu ({self.I_cu:~P.0f})",
                category=ProtectionWarning
            )
            return False

        # Check whether the circuit breaker still recognizes I_sc_min as being a
        # short-circuit (ensure magnetic tripping).
        if I_sc_min >= self.I_m_min:
            pass
        else:
            warnings.warn(
                f"I_sc_min ({I_sc_min:~P.0f} < I_m_min ({self.I_m_min:~P.0f})",
                category=ProtectionWarning
            )
            return False

        if self.I2t is not None:
            self.I2t.ito('A ** 2 * s')

            # Check whether the circuit breaker prevents excessive heating of
            # the cable in case a short-circuit with the maximum current should
            # occur.
            if self.E_t is not None:
                E_through = self.E_t.to('A ** 2 * s')
            else:
                E_through = I_sc_max ** 2 * self.t_m.to('s')

            if self.I2t >= E_through:
                pass
            else:
                warnings.warn(
                    f"Let-through energy ({E_through:~P.4e}) > "
                    f"Joule-integral (I2t) of cable ({self.I2t:~P.4e})",
                    category=ProtectionWarning
                )
                return False

            # Check whether the circuit breaker is also able to magnetically
            # interrupt I_sc_min in time to prevent too excessive heating of the
            # cable.
            t_allow = self.I2t / (I_sc_min ** 2)
            if self.t_m.to('s') <= t_allow.to('s'):
                pass
            else:
                warnings.warn(
                    f"Magnetic tripping-time limit t_m "
                    f"({self.t_m.to('s'):~P.0f}) > allowable fault duration "
                    f"({t_allow.to('s'):~P.0f}).",
                    category=ProtectionWarning
                )
                return False
            return True
        return True
    
    def __str__(self):
        s =  f"I_n     = {self.I_n.to('A'):~P.1f}\n"
        s += f"I_r     = {self.I_r.to('A'):~P.1f}\n"
        s += f"I_nf    = {self.I_nf.to('A'):~P.1f}\n"
        s += f"I_f     = {self.I_f.to('A'):~P.1f}\n"
        s += f"t_conv  = {self.t_conv.to('h'):~P.1f}\n"
        s += f"I_m_min = {self.I_m_min.to('kA'):~P.1f}\n"
        s += f"I_m_max = {self.I_m_max.to('kA'):~P.1f}\n"
        s += f"t_m     = {self.t_m.to('ms'):~P.1f}"
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
