from python_electric import Quantity, Q_
from python_electric.materials import (
    ConductorMaterial,
    CONDUCTOR_PROPS
)
from python_electric.protection import SafetyCurve
from .earthing_system import (
    EarthingSystem,
    get_max_allow_disconnect_time,
    IndirectContactProtResult
)

__all__ = [
    "LineToExposedConductivePartFault",
    "LineToExtraneousConductivePartFault",
    "check_indirect_contact",
    "check_earthing_resistance"
]


class LineToExposedConductivePartFault:
    """
    Calculation of the fault current and fault voltage due to a line-to-ground
    fault via an exposed conductive part. Uses the "simplified" calculation
    method.
    """
    def __init__(
        self,
        U_phase: Quantity,
        L: Quantity,
        S_phase: Quantity,
        S_pe: Quantity,
        conductor_material: ConductorMaterial
    ) -> None:
        """
        Creates a `LineToGroundFault` calculation object.

        Parameters
        ----------
        U_phase: Quantity
            Line-to-ground (phase) system voltage.
        L: Quantity.
            Cable length measured from the electrical distribution board to the
            electrical appliance with the exposed conductive part.
        S_phase: Quantity
            Cross-sectional area of the phase conductors.
        S_pe: Quantity
            Cross-sectional area of PE-conductor.
        conductor_material: ConductorMaterial
            Material the cable conductors are made of. See enum
            ConductorMaterials in /materials.
        """
        U_phase = U_phase.to('V').m
        L = L.to('m').m
        S_phase = S_phase.to('mm ** 2').m
        S_pe = S_pe.to('mm ** 2').m

        self.U_phase = U_phase
        self.L = L
        self.S_phase = S_phase
        self.S_pe = S_pe
        conductor_material = CONDUCTOR_PROPS[conductor_material]
        self.rho = 1.25 * conductor_material.resistivity20
        self.R_phase = self._get_R(L, S_phase, self.rho)
        self.R_pe = self._get_R(L, S_pe, self.rho)

    @staticmethod
    def _get_R(L: float, S: float, rho: float) -> float:
        R = rho * L / S
        if 120.0 < S <= 150.0:
            return 1.15 * R
        elif 150.0 < S <= 185.0:
            return 1.20 * R
        elif S > 185.0:
            return 1.25 * R
        else:
            return R

    @property
    def I_fault(self) -> Quantity:
        """
        Returns the short-circuit current due to a line-to-ground fault at the
        end of the cable (see attribute `self.L`).

        Returns
        -------
        Quantity
        """
        I_f = 0.8 * self.U_phase / (self.R_phase + self.R_pe)
        return Q_(I_f, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the fault voltage (touch potential rise) between ground (0 V)
        and the exposed conductive part.

        Returns
        -------
        Quantity
        """
        I_f = self.I_fault.to('A').m
        U_f = self.R_pe * I_f
        return Q_(U_f, 'V')

    def L_max(self, I_m_cb: Quantity) -> Quantity:
        """
        Returns the maximum allowable cable length, so that the circuit breaker
        with magnetic tripping current `I_m_cb` is just able to clear the fault
        in "short-circuit" time (the magnetic tripping-time limit `t_m`).
        """
        I_m_cb = I_m_cb.to('A').m
        m = self.S_phase / self.S_pe
        L_max = 0.8 * self.U_phase * self.S_phase / (self.rho * I_m_cb * (1 + m))
        return Q_(L_max, 'm')


class LineToExtraneousConductivePartFault:
    """
    Calculation of the fault current and fault voltage due to a line-to-ground
    fault between a phase conductor and an extraneous conductive part. The fault
    current now flows through earth instead of the PE-conductor(s).
    """
    def __init__(
        self,
        U_phase: Quantity,
        R_e: Quantity,
        R_e_extr: Quantity = Q_(5, 'ohm'),
        skin_condition: str = "BB2"
    ) -> None:
        """
        Creates a `LineToGroundFault` calculation object.

        Parameters
        ----------
        U_phase: Quantity
            Line-to-ground (phase) system voltage.
        R_e: Quantity
            Earth-spreading resistance of the low-voltage network.
        R_e_extr: Quantity
            Contact resistance of any extraneous conductive part to earth.
            By default, set to 5 ohm (according to AREI art. 80.03).
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the human skin condition: either dry (BB1) or
            wet (BB2).
        """
        U_phase = U_phase.to('V').m
        R_e = R_e.to('ohm').m
        R_e_extr = R_e_extr.to('ohm').m

        self.U_phase = U_phase
        self.R_e = R_e
        self.R_e_extr = R_e_extr
        self.skin_condition = skin_condition

    @property
    def I_fault(self) -> Quantity:
        """
        Returns the fault current that flows through the network
        earthing-resistance due to an insulation fault with an extraneous
        conductive part.

        Returns
        -------
        Quantity
        """
        I_f = self.U_phase / (self.R_e + self.R_e_extr)
        return Q_(I_f, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the fault voltage (touch potential rise) between earth (0 V) and
        the exposed conductive parts in the network due to a fault current
        flowing through the earthing-resistance of the network.

        Returns
        -------
        Quantity
        """
        I_f = self.I_fault.to('A').m
        U_f = self.R_e * I_f
        return Q_(U_f, 'V')

    @property
    def R_e_max(self) -> Quantity:
        """
        Returns the maximum allowable earth-spreading resistance of the network
        so that should an insulation fault occur through an extraneous
        conductive part, the conventional absolute limit-voltage cannot be
        exceeded.

        Returns
        -------
        Quantity
        """
        safety_curve = SafetyCurve(
            voltage_type="AC",
            skin_condition=self.skin_condition
        )
        U_l = safety_curve.conv_abs_limit_voltage.to('V').m
        R_e_max = self.R_e_extr * U_l / (self.U_phase - U_l)
        return Q_(R_e_max, 'ohm')

    @property
    def t_contact_max(self) -> Quantity:
        """
        Returns the allowable maximum contact duration with the fault voltage
        according to the applicable safety curve.

        Returns
        -------
        Quantity
        """
        safety_curve = SafetyCurve(
            voltage_type="AC",
            skin_condition=self.skin_condition
        )
        t_contact_max = safety_curve.max_contact_duration(self.U_fault)
        return t_contact_max


def check_indirect_contact(
    U_phase: Quantity,
    L_cable: Quantity,
    S_phase: Quantity,
    conductor_material: ConductorMaterial,
    I_m_cb: Quantity | None = None,
    S_pe: Quantity | None = None,
    skin_condition: str = "BB2",
    final_circuit: bool = True
) -> IndirectContactProtResult:
    """
    Determines the requirements so that a circuit breaker would also protect
    against indirect contact in a low-voltage distribution network with
    TN-earthing system when an insulation fault with an exposed conductive
    part should occur.

    Parameters
    ----------
    U_phase: Quantity
        Line-to-ground (phase) system voltage.
    L_cable: Quantity.
        Cable length between the electrical distribution board and the exposed
        conductive part.
    S_phase: Quantity
        Cross-sectional area of the loaded phase conductors.
    conductor_material: ConductorMaterial
        See enum ConductorMaterials in /materials.
    I_m_cb: Quantity, optional
        Maximum short-circuit tripping current of the circuit breaker. If
        specified (not None), the maximum allowable cable length L_max is
        calculated such that the short-circuit current due to a fault occuring
        at the cable end becomes equal to the circuit breaker's maximum
        short-circuit tripping current.
    S_pe: Quantity, optional
        Cross-sectional area of the PE-conductor. If None, `S_pe` is set equal
        to `S_phase`.
    skin_condition: str, {"BB1", "BB2" (default)}
        Code that identifies the human skin condition: either dry (BB1) or wet
        (BB2).
    final_circuit: bool, default True
        Indicates that the cable directly feeds a consumer (of which the
        nominal current does not exceed 32 A).

    Returns
    -------
    IndirectContactProtResult
    """
    fault = LineToExposedConductivePartFault(
        U_phase=U_phase,
        L=L_cable,
        S_phase=S_phase,
        S_pe=S_phase if S_pe is None else S_pe,
        conductor_material=conductor_material
    )

    if I_m_cb is not None:
        L_max = fault.L_max(I_m_cb)
        R_pe_max = (fault.U_fault / I_m_cb).to('ohm')
    else:
        L_max = None
        R_pe_max = None

    safety_curve = SafetyCurve(
        voltage_type="AC",
        skin_condition=skin_condition
    )
    t_c_max = safety_curve.max_contact_duration(fault.U_fault)

    if final_circuit:
        t_c_max_iec = get_max_allow_disconnect_time(U_phase, EarthingSystem.TN)
        t_c_max = min(t_c_max.to('ms'), t_c_max_iec.to('ms'))

    return IndirectContactProtResult(
        I_f=fault.I_fault,
        U_f=fault.U_fault,
        L_max=L_max,
        t_contact_max=t_c_max,
        R_pe_max=R_pe_max
    )


def check_earthing_resistance(
    U_phase: Quantity,
    R_e: Quantity,
    R_e_extr: Quantity = Q_(5, 'ohm'),
    skin_condition: str = "BB2",
    final_circuit: bool = True
) -> IndirectContactProtResult:
    """
    Determines the requirements regarding the earthing resistance of the
    TN-earthed network by considering an insulation fault via an extraneous
    conductive part.

    Parameters
    ----------
    U_phase: Quantity
        Phase voltage (line-to-ground voltage) of the network.
    R_e: Quantity
        Earth-spreading resistance of the network.
    R_e_extr: Quantity, optional
        Earth-contact resistance between ground and any extraneous conductive
        part. By default, set to 5 ohm (according to AREI art. 80.03).
    skin_condition: str, {"BB1", "BB2" (default)}
        Code that identifies the human skin condition: either dry (BB1) or wet
        (BB2).
    final_circuit: bool, default True
        Indicates that the cable feeds a consumer of which the nominal current
        does not exceed 32 A.

    Returns
    -------
    IndirectContactProtResult
    """
    fault = LineToExtraneousConductivePartFault(U_phase, R_e, R_e_extr, skin_condition)
    if final_circuit:
        t_c_max_iec = get_max_allow_disconnect_time(U_phase, EarthingSystem.TN)
        t_c_max = min(fault.t_contact_max.to('ms'), t_c_max_iec.to('ms'))
    else:
        t_c_max = fault.t_contact_max
    return IndirectContactProtResult(
        I_f=fault.I_fault,
        U_f=fault.U_fault,
        R_e_max=fault.R_e_max,
        t_contact_max=t_c_max
    )
