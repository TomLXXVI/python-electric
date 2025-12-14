from python_electric import Quantity
from python_electric.materials import (
    ConductorMaterials,
    CONDUCTOR_MATERIALS
)
from python_electric.protection import SafetyCurve
from .earthing_system import (
    EarthingSystem,
    get_max_allow_disconnect_time,
    IndirectContactProtectionResult
)

__all__ = [
    "LineToExposedConductivePartFault",
    "LineToExtraneousConductivePartFault",
    "check_indirect_contact",
    "check_earthing_resistance"
]

Q_ = Quantity


class LineToExposedConductivePartFault:
    """
    Implements the calculation of the fault current and fault voltage due to
    a line-to-ground fault in an electrical appliance with an exposed conductive
    part using the simplified method. The fault current flows through the
    earth protective conductor.
    """
    def __init__(
        self,
        U_phase: Quantity,
        L: Quantity,
        S_phase: Quantity,
        S_pe: Quantity,
        conductor_material: ConductorMaterials
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
        conductor_material: ConductorMaterials
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
        conductor_material = CONDUCTOR_MATERIALS[conductor_material]
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
        Returns the current due to a line-to-ground fault at the end of the
        given cable length (see attribute `self.L`).

        Returns
        -------
        Quantity
        """
        I_f = 0.8 * self.U_phase / (self.R_phase + self.R_pe)
        return Q_(I_f, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the voltage (touch potential rise) between ground (0 V) and the
        exposed conductive part of the electrical appliance.

        Returns
        -------
        Quantity
        """
        I_f = self.I_fault.to('A').m
        U_f = self.R_pe * I_f
        return Q_(U_f, 'V')

    def L_max(self, I_m_cb: Quantity) -> Quantity:
        """
        Returns the maximum length of the cable calculated so that the circuit
        breaker with magnetic tripping current `I_m_cb` is still able to clear
        the line-to-ground fault current in "short-circuit" time (i.e. the
        magnetic tripping time upper limit, `t_m_lim`).
        """
        I_m_cb = I_m_cb.to('A').m
        m = self.S_phase / self.S_pe
        L_max = 0.8 * self.U_phase * self.S_phase / (self.rho * I_m_cb * (1 + m))
        return Q_(L_max, 'm')


class LineToExtraneousConductivePartFault:
    """
    Implements the calculation of the fault current and fault voltage due to
    a line-to-ground fault between a phase conductor and an extraneous
    conductive part. The fault current flows through earth, instead of the
    earth-protective conductor of the low-voltage distribution network.
    """
    def __init__(
        self,
        U_phase: Quantity,
        R_b: Quantity,
        R_e: Quantity = Q_(5, 'ohm'),
        skin_condition: str = "BB2"
    ) -> None:
        """
        Creates a `LineToGroundFault` calculation object.

        Parameters
        ----------
        U_phase: Quantity
            Line-to-ground (phase) system voltage.
        R_b: Quantity
            Measured earth spreading resistance of the low-voltage distribution
            network.
        R_e: Quantity
            Contact resistance of any extraneous conductive part to earth.
            By default, it is set to 5 ohm (according to AREI art. 80.03).
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the condition of the human skin: either dry
            (BB1) or wet (BB2).
        """
        U_phase = U_phase.to('V').m
        R_b = R_b.to('ohm').m
        R_e = R_e.to('ohm').m

        self.U_phase = U_phase
        self.R_b = R_b
        self.R_e = R_e
        self.skin_condition = skin_condition

    @property
    def I_fault(self) -> Quantity:
        """
        Returns the fault current that flows through the earthing resistance of
        the low-voltage distribution network due to an insulation fault with an
        extraneous conductive part.

        Returns
        -------
        Quantity
        """
        I_f = self.U_phase / (self.R_b + self.R_e)
        return Q_(I_f, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the voltage (touch potential rise) between earth (0 V) and the
        exposed conductive parts of the low-voltage distribution network due to
        a fault current flowing through the earthing resistance of the
        low-voltage network.

        Returns
        -------
        Quantity
        """
        I_f = self.I_fault.to('A').m
        U_f = self.R_b * I_f
        return Q_(U_f, 'V')

    @property
    def R_b_max(self) -> Quantity:
        """
        Returns the maximum allowable earth spreading resistance of the
        low-voltage distribution network so that when an insulation fault to 
        earth happens along an extraneous conductive part, the conventional 
        absolute limit voltage cannot be exceeded.

        Returns
        -------
        Quantity
        """
        safety_curve = SafetyCurve(
            voltage_type="AC",
            skin_condition=self.skin_condition
        )
        U_l = safety_curve.conv_abs_limit_voltage.to('V').m
        R_b_max = self.R_e * U_l / (self.U_phase - U_l)
        return Q_(R_b_max, 'ohm')

    @property
    def t_c_max(self) -> Quantity:
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
        t_c_max = safety_curve.max_contact_duration(self.U_fault)
        return t_c_max


def check_indirect_contact(
    U_phase: Quantity,
    L_cable: Quantity,
    S_phase: Quantity,
    conductor_material: ConductorMaterials,
    I_m_cb: Quantity | None = None,
    S_pe: Quantity | None = None,
    skin_condition: str = "BB2",
    final_circuit: bool = True
) -> IndirectContactProtectionResult:
    """
    Determines the requirements so that a circuit breaker would also protect
    against indirect contact in a low-voltage distribution network with
    TN-earthing system when an insulation fault to an exposed conductive
    part should happen.

    Parameters
    ----------
    U_phase: Quantity
        Line-to-ground (phase) system voltage.
    L_cable: Quantity.
        Cable length measured from the electrical distribution board to the
        electrical appliance with the exposed conductive part.
    S_phase: Quantity
        Cross-sectional area of the loaded phase conductors.
    conductor_material: ConductorMaterials
        Material the cable conductors are made of. See enum
        ConductorMaterials in /materials.
    I_m_cb: Quantity, optional
        Maximum limit of the magnetic tripping current of the circuit
        breaker. If not None, the maximum allowable length of the cable is
        calculated so that the minimum short-circuit current at the end of
        the cable equals the maximum limit of the magnetic trip current of
        the circuit breaker.
    S_pe: Quantity, optional
        Cross-sectional area of the earth protective conductor. If None, `S_pe`
        is set equal to `S_phase`.
    skin_condition: str, {"BB1", "BB2" (default)}
        Code that identifies the condition of the human skin: either dry
        (BB1) or wet (BB2).
    final_circuit: bool, default True
        Indicates that the cable directly feeds a consumer (of which the
        nominal current does not exceed 32 A).

    Returns
    -------
    IndirectContactProtectionResult
        I_f: Quantity | None
            Fault current.
        U_f: Quantity | None
            Voltage (touch potential rise) between earth (0 V) and an
            exposed conductive part of the low-voltage distribution network.
        L_max: Quantity | None
            Maximum allowable length of the cable between the distribution
            board and the appliance/sub-distribution board, so that the
            minimum short-circuit current calculated at the end of the cable
            would equal the maximum limit of the magnetic trip current of
            the circuit breaker.
        t_c_max: Quantity | None
            Maximum allowable contact duration with touch voltage according
            to the applicable safety curve.
        R_pe_max: Quantity | None
            Maximum allowable resistance of the PE-conductor between the
            main earthing terminal and the most distant exposed conductive
            part in the distribution network, so that the fault current
            would equal the maximum limit of the magnetic trip current of
            the circuit breaker just upstream of the most distant exposed
            conductive part.
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
    return IndirectContactProtectionResult(
        I_f=fault.I_fault,
        U_f=fault.U_fault,
        L_max=L_max,
        t_c_max=t_c_max,
        R_pe_max=R_pe_max
    )


def check_earthing_resistance(
    U_phase: Quantity,
    R_b: Quantity,
    R_e: Quantity = Q_(5, 'ohm'),
    skin_condition: str = "BB2",
    final_circuit: bool = True
) -> Quantity:
    """
    Determines the requirements regarding the earthing resistance of the
    low-voltage distribution network by considering an insulation fault via
    an extraneous conductive part.

    Parameters
    ----------
    U_phase: Quantity
        Phase voltage (line-to-ground voltage) of the low-voltage
        distribution network.
    R_b: Quantity
        Measured earth spreading resistance of the low-voltage distribution
        network earthing system.
    R_e: Quantity, optional
        The earth contact resistance between ground and any extraneous
        conductive part. By default, this is set to 5 ohm (according to AREI
        art. 80.03).
    skin_condition: str, {"BB1", "BB2" (default)}
        Code that identifies the condition of the human skin: either dry
        (BB1) or wet (BB2).
    final_circuit: bool, default True
        Indicates that the cable is feeding an electrical consumer (of which
        the nominal current does not exceed 32 A).

    Returns
    -------
    IndirectContactProtectionResult.
        I_f: Quantity
            Fault current.
        U_f: Quantity | None
            Voltage (touch potential rise) between earth (0 V) and an
            exposed conductive part of the low-voltage distribution network.
        t_c_max: Quantity | None
            Maximum allowable contact duration with touch voltage according
            to the applicable safety curve.
        R_b_max: Quantity | None
            Maximum allowable earth spreading resistance of the low-voltage
            distribution network.
    """
    fault = LineToExtraneousConductivePartFault(U_phase, R_b, R_e, skin_condition)
    if final_circuit:
        t_c_max_iec = get_max_allow_disconnect_time(U_phase, EarthingSystem.TN)
        t_c_max = min(fault.t_c_max.to('ms'), t_c_max_iec.to('ms'))
    else:
        t_c_max = fault.t_c_max
    return IndirectContactProtectionResult(
        I_f=fault.I_fault,
        U_f=fault.U_fault,
        R_b_max=fault.R_b_max,
        t_c_max=t_c_max
    )
