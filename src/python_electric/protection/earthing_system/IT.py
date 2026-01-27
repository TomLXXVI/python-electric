import math

from ... import Quantity, Q_
from ... import calc
from ...materials import ConductorMaterial
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

PI = math.pi


class ITSingleFault:
    """
    Calculation of the fault current and fault voltage in an IT-earthed network
    due to a single line-to-ground fault via an exposed conductive part. The
    fault loop is closed by parasitic capacitances of the network cabling.
    """
    def __init__(
        self,
        U_phase: Quantity,
        Z_n: Quantity,
        R_e: Quantity,
        L_tot: Quantity,
        c_cable: Quantity = Q_(0.25, 'µF / km'),
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
            Earth-spreading resistance of the consumer installation, including
            the resistance of the PE-conductors between the earthing main
            terminal and the exposed conductive part in the fault loop.
        L_tot:
            Approximate total cable length of the low-voltage network in the
            building.
        c_cable: Quanity, default 0.25 µF
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
        any exposed conductive part in the network.

        Returns
        -------
        Quantity
        """
        a = complex(1, self.omega * self._C * self.Z_n)
        I_f = a * -self._E / (self.Z_n + a * self.R_f)
        I_f_mag = abs(I_f)
        return Q_(I_f_mag, 'A')

    @property
    def U_fault(self) -> Quantity:
        """
        Returns the fault voltage (i.e. the touch potential rise) between ground
        (0 V) and the exposed conductive part due to the earthing resistance
        `R_e` of the consumer installation.

        Returns
        -------
        Quantity
        """
        I_f = self.I_fault.to('A').m
        U_f = self.R_e * I_f
        return Q_(U_f, 'V')


class ITDoubleFault:
    """
    Calculation of the fault current and fault voltage in an IT-earthed network
    where a fault loop arises due to two separate line-to-ground faults in the
    network.

    To analyze the double line-to-ground fault in an IT network, it is assumed
    that for each circuit in the network, an identical circuit is present in the
    same network.
    """
    def __init__(
        self,
        U_phase: Quantity,
        L: Quantity,
        S_phase: Quantity,
        S_pe: Quantity,
        conductor_material: ConductorMaterial,
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
            Cable length between the electrical distribution board and the
            exposed conductive part.
        S_phase: Quantity
            Cross-sectional area of the phase conductors.
        S_pe: Quantity
            Cross-sectional area of PE-conductor.
        conductor_material: ConductorMaterial
            See enum ConductorMaterials in materials.
        neutral_distributed: bool, default True
            Indicates whether the neutral conductor is distributed in the
            network.
        R_e: Quantity | None
            Earth-spreading resistance of an individual earth electrode in the
            consumer installation. Must be set when the two exposed conductive
            parts in the fault loop are connected to different, separate earth
            electrodes. It is assumed that all earth electrodes in the consumer
            installation have the same earth-spreading resistance.
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
        Returns the fault voltage (i.e. the touch potential rise) between ground
        (0 V) and the faulted exposed conductive part.

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
        Returns the maximum allowable length of the cable, so that the circuit
        breaker with short-circuit tripping current `I_m_cb` is just able to
        clear fault current in "short-circuit" time.
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
    Calculation of the fault current and fault voltage due to two separate
    line-to-ground faults in an IT-earthed network. One of the two
    faults occurs between a phase conductor and an extraneous conductive part.
    The fault current passes through earth and through the earthing resistance
    of the consumer installation.
    """
    def __init__(
        self,
        U_phase: Quantity,
        R_e: Quantity,
        R_e_extr: Quantity = Q_(5, 'ohm'),
        skin_condition: str = "BB2",
        neutral_distributed: bool = True
    ) -> None:
        """
        Creates a `ITExtraneousConductivePartFault` calculation object.

        Parameters
        ----------
        U_phase: Quantity
            Line-to-ground (phase) system voltage.
        R_e: Quantity
            Earth-spreading resistance of the consumer installation.
        R_e_extr: Quantity
            Contact resistance of any extraneous conductive part to earth.
            By default, set to 5 ohm (according to AREI art. 80.03).
        skin_condition: str, {"BB1", "BB2" (default)}
            Code that identifies the human skin condition: either dry (BB1) or
            wet (BB2).
        neutral_distributed: bool, default True
            Indicates whether the neutral conductor is distributed in the
            network.
        """
        if not neutral_distributed:
            U = math.sqrt(3) * U_phase
        else:
            U = U_phase
        self._TN_fault = TN.LineToExtraneousConductivePartFault(
            U_phase=U,
            R_e=R_e,
            R_e_extr=R_e_extr,
            skin_condition=skin_condition
        )

    @property
    def I_fault(self) -> Quantity:
        return self._TN_fault.I_fault

    @property
    def U_fault(self) -> Quantity:
        return self._TN_fault.U_fault

    @property
    def R_e_max(self) -> Quantity:
        return self._TN_fault.R_e_max

    @property
    def t_contact_max(self) -> Quantity:
        return self._TN_fault.t_contact_max


def check_indirect_contact(
    U_phase: Quantity,
    L_cable: Quantity,
    S_phase: Quantity,
    conductor_material: ConductorMaterial,
    I_m_cb: Quantity | None = None,
    S_pe: Quantity | None = None,
    skin_condition: str = "BB2",
    final_circuit: bool = True,
    neutral_distributed: bool = True,
    R_e: Quantity | None = None
) -> IndirectContactProtectionResult:
    """
    Determines the requirements so that a circuit breaker would also protect
    against indirect contact in a low-voltage network with IT-earthing system in
    case two insulation faults occur in the network. It is assumed that the
    second insulation fault via an exposed conductive part occurs in a circuit
    being identical to the circuit under investigation.

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
    conductor_material: ConductorMaterial
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
        Earth-spreading resistance of the earth electrodes in the consumer
        installation. Must be set when two exposed conductive parts in the
        fault loop may be connected to different, separate earthelectrodes. It
        is assumed that all earth electrodes in the consumer installation have
        the same earth-spreading resistance.

    Returns
    -------
    IndirectContactProtectionResult
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
        t_contact_max=t_c_max,
        R_pe_max=R_pe_max
    )


def check_earthing_resistance(
    U_phase: Quantity,
    R_e: Quantity,
    R_e_extr: Quantity = Q_(5, 'ohm'),
    skin_condition: str = "BB2",
    final_circuit: bool = True,
    neutral_distributed: bool = True
) -> IndirectContactProtectionResult:
    """
    Determines the requirements regarding the earthing resistance of the
    consumer installation by considering a second insulation fault via an
    extraneous conductive part.

    Parameters
    ----------
    U_phase: Quantity
        Phase voltage (line-to-ground voltage) of the low-voltage
        distribution network.
    R_e: Quantity
        Earth spreading resistance of the consumer installation.
    R_e_extr: Quantity, optional
        Earth-contact resistance of any extraneous conductive part. By default,
        this is set to 5 ohm (according to AREI art. 80.03).
    skin_condition: str, {"BB1", "BB2" (default)}
        Code that identifies the human skin condition: either dry (BB1) or wet
        (BB2).
    final_circuit: bool, default True
        Indicates the cable is feeding an electrical consumer (of which the
        nominal current does not exceed 32 A).
    neutral_distributed: bool, default True
        Indicates whether the neutral conductor is distributed in the network.

    Returns
    -------
    IndirectContactProtectionResult
    """
    fault = ITExtraneousConductivePartFault(
        U_phase,
        R_e, R_e_extr,
        skin_condition,
        neutral_distributed
    )
    if final_circuit:
        t_c_max_iec = get_max_allow_disconnect_time(U_phase, EarthingSystem.TN)
        t_c_max = min(fault.t_contact_max.to('ms'), t_c_max_iec.to('ms'))
    else:
        t_c_max = fault.t_contact_max
    return IndirectContactProtectionResult(
        I_f=fault.I_fault,
        U_f=fault.U_fault,
        R_e_max=fault.R_e_max,
        t_contact_max=t_c_max
    )
