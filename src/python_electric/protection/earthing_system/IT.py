import math

from ... import Quantity
from ... import calc
from ...materials import ConductorMaterials
from ..safety_curve import SafetyCurve
from .earthing_system import (
    IndirectContactProtectionResult,
    get_max_allow_disconnect_time,
    EarthingSystem
)
from . import TN

__all__ = [
    "ITSingleFault",
    "ITDoubleFault",
    "ITExtraneousConductivePartFault",
    "check_indirect_contact",
    "check_earthing_resistance"
]

Q_ = Quantity
PI = math.pi


class ITSingleFault:
    """
    Implements the calculation of the fault current and fault voltage in a
    low-voltage, IT-earthed distribution network due to a single line-to-ground
    fault with an exposed conductive part. The fault loop is closed through
    parasitic capacitances of the network cabling.
    """
    def __init__(
        self,
        U_phase: Quantity,
        Z_n: Quantity,
        R_e: Quantity,
        L_tot: Quantity,
        c_cable: Quantity = Q_(25, 'µF / km'),
        neutral_distributed: bool = True,
        C_filter: Quantity = Q_(0.0, 'F'),
        R_f: Quantity = Q_(0.0, 'ohm'),
        f: Quantity = Q_(50, 'Hz')
    ) -> None:
        """
        Creates a `ITSingleFault` calculation object.

        Parameters
        ----------
        U_phase: Quantity
            Line-to-ground (phase) system voltage.
        Z_n: Quantity
            Earthing impedance of the low-voltage, IT-earthed distribution
            network.
        R_e: Quantity
            Separate earthing resistance of the exposed conductive parts in
            the installation (including the total resistance of the
            PE-conductors between the earthing main terminal and the exposed
            conductive part that is part of the fault loop).
        L_tot:
            Approximate total cable length of the low-voltage network in the
            building.
        c_cable: Quanity
            Specific capacitance (i.e. per unit length of cable) of the
            conductors with respect to ground.
        neutral_distributed: bool, default True
            Indicates whether the neutral conductor is also distributed in the
            network or not.
        C_filter: Quantity, default 0 F
            Total capacitance of all capacitive filters in electronic devices
            that conduct current harmonics to earth.
        R_f: Quantity, default 0 ohm
            Resistance of the insulation fault between a phase of the network
            and an exposed conductive part connected to the PE-conductor. 
        f: Quantity, default 50 Hz
            Network voltage frequency.
        
        Returns
        -------
        None
        """
        self.U_phase = U_phase.to('V').m
        self.Z_n = Z_n.to('ohm').m
        self.R_e = R_e.to('ohm').m
        self.L_tot = L_tot.to('km').m
        self.c_cable = c_cable.to('F / km').m
        self.neutral_distributed = neutral_distributed
        self.C_filter = C_filter.to('F').m
        self.R_f = R_f.to('ohm').m
        self.f = f.to('Hz').m
        self.omega = 2 * PI * self.f

        self._E = calc.phasor(-self.U_phase, 0.0)
        self._C = self._calc_C()

    def _calc_C(self) -> float:
        C_cable = self.c_cable * self.L_tot
        C = 3 * C_cable
        if self.neutral_distributed:
            C += C_cable
        C += self.C_filter
        return C

    @property
    def I_fault(self) -> Quantity:
        """
        Returns the fault current due to a single line-to-ground fault through
        an exposed conductive part somewhere in the installation.

        Returns
        -------
        Quantity
        """
        a = complex(0, self.omega * self._C)
        b = a * self.Z_n
        I_f = -self._E * b / (self.Z_n + self.R_f * (1 + b))
        I_f_mag = abs(I_f)
        return Q_(I_f_mag, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the voltage (touch potential rise) between ground (0 V) and the
        exposed conductive part, due to the earthing resistance `R_e` of the
        installation.

        Returns
        -------
        Quantity
        """
        I_f = self.I_fault.to('A').m
        U_f = self.R_e * I_f
        return Q_(U_f, 'V')


class ITDoubleFault:
    """
    Implements the calculation of the fault current and voltage in a low-voltage
    IT-earthed distribution network in case of a fault loop wherein two
    line-to-ground faults exist.

    Notes
    -----
    To analyze a double line-to-ground fault, it is assumed that for any circuit
    in the network an identical circuit is present.
    """
    def __init__(
        self,
        U_phase: Quantity,
        L: Quantity,
        S_phase: Quantity,
        S_pe: Quantity,
        conductor_material: ConductorMaterials,
        neutral_distributed: bool = True,
        R_e: Quantity | None = None
    ) -> None:
        """
        Creates a `ITDoubleFault` calculation object.

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
        neutral_distributed: bool, default True
            Indicates whether the neutral conductor is also distributed in the
            network or not.
        R_e: Quantity | None
            Spreading resistance of the earth electrodes if the two exposed
            conductive parts that are part of the fault loop are both connected
            to different, separate earth electrodes. It is assumend that all
            earth electrodes in the installation have the same spreading
            resistance.
        """
        if not neutral_distributed:
            U = math.sqrt(3) * U_phase
        else:
            U = U_phase
        self._TN_fault = TN.LineToExposedConductivePartFault(
            U_phase=U,
            L=2 * L,
            S_phase=S_phase,
            S_pe=S_pe,
            conductor_material=conductor_material
        )
        if R_e is None:
            self._R_e = None
        else:
            self._R_e = R_e.to('ohm').m
            R_phase = self._TN_fault.R_phase
            R_pe = self._TN_fault.R_pe
            U = U.to('V').m
            self._I_fault = 0.8 * U / (R_phase + R_pe + 2 * self._R_e)
            self._U_fault = (R_pe + self._R_e) * self._I_fault

    @property
    def I_fault(self) -> Quantity:
        """
        Returns the fault current in the circuit due to a double line-to-ground
        fault.

        Returns
        -------
        Quantity
        """
        if self._R_e is None:
            return self._TN_fault.I_fault
        else:
            return Q_(self._I_fault, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the voltage (touch potential rise) between ground (0 V) and the
        faulted exposed conductive part in the circuit.

        Returns
        -------
        Quantity
        """
        if self._R_e is None:
            return self._TN_fault.U_fault
        else:
            return Q_(self._U_fault, 'V')

    def L_max(self, I_m_cb: Quantity) -> Quantity:
        """
        Returns the maximum length of the cable calculated so that the circuit
        breaker with magnetic tripping current `I_m_cb` is still able to clear
        the line-to-ground fault current in "short-circuit" time (i.e. the
        magnetic tripping time upper limit, `t_m_lim`).
        """
        I_m_cb = I_m_cb.to('A').m
        rho = self._TN_fault.rho
        R_e = 0.0 if self._R_e is None else self._R_e
        U = self._TN_fault.U_phase
        S_ph = self._TN_fault.S_phase
        S_pe = self._TN_fault.S_pe
        n = 0.8 * U / I_m_cb - 2 * R_e
        d = rho * (1 / S_ph + 1 / S_pe)
        L_max = n / d
        return Q_(L_max, 'm')


class ITExtraneousConductivePartFault:
    """
    Implements the calculation of the fault current and fault voltage due to
    a double line-to-ground fault in a low-voltage, IT-earthed distribution
    network where one of the two faults exists between a phase conductor and an
    extraneous conductive part. The fault current flows through earth towards
    the earthing resistance of the installation.
    """
    def __init__(
        self,
        U_phase: Quantity,
        R_b: Quantity,
        R_e: Quantity = Q_(5, 'ohm'),
        skin_condition: str = "BB2",
        neutral_distributed: bool = True
    ) -> None:
        """
        Creates a `ITExtraneousConductivePartFault` calculation object.

        Parameters
        ----------
        U_phase: Quantity
            Line-to-ground (phase) system voltage.
        R_b: Quantity
            Measured earth spreading resistance of the installation.
        R_e: Quantity
            Contact resistance of any extraneous conductive part to earth.
            By default, it is set to 5 ohm (according to AREI art. 80.03).
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the condition of the human skin: either dry
            (BB1) or wet (BB2).
        neutral_distributed: bool, default True
            Indicates whether the neutral conductor is also distributed in the
            network or not.
        """
        if not neutral_distributed:
            U = math.sqrt(3) * U_phase
        else:
            U = U_phase
        self._TN_fault = TN.LineToExtraneousConductivePartFault(
            U_phase=U,
            R_b=R_b,
            R_e=R_e,
            skin_condition=skin_condition
        )

    @property
    def I_fault(self) -> Quantity:
        return self._TN_fault.I_fault

    @property
    def U_fault(self) -> Quantity:
        return self._TN_fault.U_fault

    @property
    def R_b_max(self) -> Quantity:
        return self._TN_fault.R_b_max

    @property
    def t_c_max(self) -> Quantity:
        return self._TN_fault.t_c_max


def check_indirect_contact(
    U_phase: Quantity,
    L_cable: Quantity,
    S_phase: Quantity,
    conductor_material: ConductorMaterials,
    I_m_cb: Quantity | None = None,
    S_pe: Quantity | None = None,
    skin_condition: str = "BB2",
    final_circuit: bool = True,
    neutral_distributed: bool = True,
    R_e: Quantity | None = None
) -> IndirectContactProtectionResult:
    """
    Determines the requirements so that the circuit breaker of a circuit would
    also protect against indirect contact in a low-voltage distribution network
    with IT-earthing system where, by assumption, a second insulation fault
    through an exposed conductive part occurs in circuit which is identical to
    the circuit under investigation.

    Parameters
    ----------
    U_phase: Quantity
        Line-to-ground (phase) system voltage.
    L_cable: Quantity.
        Cable length measured from the electrical distribution board to the
        electrical appliance with the exposed conductive part in the circuit
        under investigation.
    S_phase: Quantity
        Cross-sectional area of the loaded phase conductors.
    conductor_material: ConductorMaterials
        Material the cable conductors are made of. See enum
        ConductorMaterials in /materials.
    I_m_cb: Quantity, optional
        Maximum limit of the magnetic tripping current of the circuit
        breaker in the circuit under investigation. If not None, the maximum
        allowable length of the cable is calculated so that the minimum
        short-circuit current at the end of the cable equals the maximum limit
        of the magnetic trip current of the circuit breaker.
    S_pe: Quantity, optional
        Cross-sectional area of the earth protective conductor. If None, `S_pe`
        is set equal to `S_phase`.
    skin_condition: str, {"BB1", "BB2" (default)}
        Code that identifies the condition of the human skin: either dry
        (BB1) or wet (BB2).
    final_circuit: bool, default True
        Indicates that the cable directly feeds a consumer (of which the
        nominal current does not exceed 32 A).
    neutral_distributed: bool, default True
        Indicates whether the neutral conductor is also distributed in the
        network or not.
    R_e: Quantity | None
        Spreading resistance of the earth electrodes if the two exposed
        conductive parts that are part of the fault loop are both connected
        to different, separate earth electrodes. It is assumend that all
        earth electrodes in the installation have the same spreading
        resistance.

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
    fault = ITDoubleFault(
        U_phase=U_phase,
        L=L_cable,
        S_phase=S_phase,
        S_pe=S_phase if S_pe is None else S_pe,
        conductor_material=conductor_material,
        neutral_distributed=neutral_distributed,
        R_e=R_e
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
    final_circuit: bool = True,
    neutral_distributed: bool = True
) -> Quantity:
    """
    Determines the requirements regarding the earthing resistance of the
    consumer installation by considering a second insulation fault via an
    extraneous conductive part.

    Parameters
    ----------
    U_phase: Quantity
        Phase voltage (line-to-ground voltage) of the low-voltage
        distribution network.
    R_b: Quantity
        Measured earth spreading resistance of the consumer installation.
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
    neutral_distributed: bool, default True
        Indicates whether the neutral conductor is also distributed in the
        network or not.

    Returns
    -------
    IndirectContactProtectionResult.
        I_f: Quantity
            Fault current.
        U_f: Quantity | None
            Voltage (touch potential rise) between earth (0 V) and an
            exposed conductive part of the consumer installation.
        t_c_max: Quantity | None
            Maximum allowable contact duration with touch voltage according
            to the applicable safety curve.
        R_b_max: Quantity | None
            Maximum allowable earth spreading resistance of the consumer
            installation.
    """
    fault = ITExtraneousConductivePartFault(U_phase, R_b, R_e, skin_condition, neutral_distributed)
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
