from typing import Optional
import math
import cmath

from .. import Quantity
from ..materials import (
    ConductorMaterials,
    CONDUCTOR_MATERIALS
)
from ..sizing.cable.sizing import (
    CABLE_UNIT_REACTANCES,
    CableArrangement
)

__all__ = [
    "PowerGrid",
    "Transformer",
    "Generator",
    "SynchronousMotor",
    "Cable",
    "BusBar",
    "InductionMotor"
]


class Equipment:

    def __init__(self, name: str) -> None:
        self.name = name


class PowerGrid(Equipment):

    def __init__(
        self,
        line_voltage: Quantity,
        short_circuit_power: Quantity,
        R_to_X_ratio: float = 0.1,
        voltage_factor: float = 1.0,
        z0_r_factor: float = 1.0,
        z0_x_factor: float = 3.0,
        name = "grid"
    ) -> None:
        super().__init__(name)
        self.U_line = line_voltage
        self.S_sc = short_circuit_power
        self.c = voltage_factor
        self.R_to_X_ratio = R_to_X_ratio
        # Zero-sequence scaling factors
        self.z0_r_factor = z0_r_factor
        self.z0_x_factor = z0_x_factor

    @property
    def Z1(self) -> Quantity:
        U_line = self.U_line.to('V').magnitude
        S_sc = self.S_sc.to('VA').magnitude
        Z_mag = self.c * U_line ** 2 / S_sc
        X = Z_mag / math.sqrt(1 + self.R_to_X_ratio ** 2)
        R = self.R_to_X_ratio * X
        Z_ph = math.atan2(X, R)
        Z = cmath.rect(Z_mag, Z_ph)
        return Quantity(Z, 'ohm')

    @property
    def Z2(self) -> Quantity:
        """
        Negative-sequence impedance; approximated equal to Z1 for this
        component.
        """
        return self.Z1

    @property
    def Z0(self) -> Quantity:
        """
        Approximate zero-sequence impedance of the utility grid infeed
        from its positive-sequence impedance.

        The model scales the positive-sequence resistance and reactance:

            R0 ≈ z0_r_factor * R1
            X0 ≈ z0_x_factor * X1

        Default values:
            z0_r_factor = 1.0
            z0_x_factor = 3.0

        These defaults are typical for many LV/MV distribution networks
        where zero-sequence currents return through earth and neutral
        conductors, resulting in a zero-sequence reactance that is several
        times larger than the positive-sequence reactance.

        Limitations:
            - Intended as a design-level approximation only.
            - Actual Z0/Z1 ratio depends on earthing system, network
              topology, and soil resistivity.
            - For accurate short-circuit studies, Z0 should be obtained
              from the grid operator or from a detailed network model.
        """
        z1 = self.Z1.to('ohm').magnitude
        r1, x1 = z1.real, z1.imag
        r0 = self.z0_r_factor * r1
        x0 = self.z0_x_factor * x1
        return Quantity(complex(r0, x0), 'ohm')


class Transformer(Equipment):

    def __init__(
        self,
        nominal_voltage: Quantity,
        nominal_power: Quantity,
        percent_short_circuit_voltage: Quantity,
        copper_loss: Optional[Quantity] = None,
        voltage_factor: float = 1.0,
        z0_factor: float = 1.0,
        name: str = "transformer"
    ) -> None:
        super().__init__(name)
        self.U_nom = nominal_voltage
        self.S_nom = nominal_power
        self.u_sc = percent_short_circuit_voltage
        self.P_Cu = copper_loss
        self.c = voltage_factor
        # Zero-sequence scaling factor
        self.z0_factor = z0_factor

    @property
    def Z1(self) -> Quantity:
        U_nom = self.U_nom.to('V').magnitude
        S_nom = self.S_nom.to('VA').magnitude
        u_sc = self.u_sc.to('frac').magnitude
        if self.P_Cu is not None:
            P_Cu = self.P_Cu.to('W').magnitude
            R = P_Cu * (U_nom / S_nom) ** 2
        else:
            R = 0.0
        Z_mag = u_sc * U_nom ** 2 / S_nom
        X = math.sqrt(Z_mag ** 2 - R ** 2)
        sin_phi = 0.6
        x = X / (U_nom ** 2 / S_nom)
        Z_mag *= 0.95 * self.c / (1.0 + x * sin_phi)
        Z_ph = math.atan2(X, R)
        Z = cmath.rect(Z_mag, Z_ph)
        return Quantity(Z, 'ohm')

    @property
    def Z2(self) -> Quantity:
        """
        Negative-sequence impedance; approximated equal to Z1 for this
        component.
        """
        return self.Z1

    @property
    def Z0(self) -> Quantity:
        """
        Approximate zero-sequence impedance of the transformer from its
        positive-sequence impedance.

            Z0 ≈ z0_factor * Z1

        Default:
            z0_factor = 1.0

        This assumes that, on the side of interest, the transformer has
        an effectively grounded star point and that the zero-sequence
        leakage path is of the same order as the positive-sequence one.

        Important limitations:
            - This model does NOT implement vector-group logic. A delta
              winding blocks zero-sequence currents on that side, which
              corresponds to Z0 → very large (not represented here).
            - Neutral earthing impedance (neutral resistor/reactor) is
              not included explicitly. In reality:
                  Z0 = Z0,transformer + 3 * Zn
            - Use this as a rough approximation when no manufacturer
              data for Z0 are available and the neutral is effectively
              grounded on the considered side.
            - For detailed studies, a dedicated transformer model with
              vector group and neutral earthing should be used.
        """
        z1 = self.Z1.to('ohm').magnitude
        return Quantity(self.z0_factor * z1, 'ohm')


class Generator(Equipment):

    def __init__(
        self,
        nominal_voltage: Quantity,
        nominal_power: Quantity,
        per_unit_reactance: Quantity,
        power_factor: float = 0.8,
        R_to_X_ratio: float = 0.15,
        voltage_factor: float = 1.0,
        z0_factor: float = 1.0,
        name: str = "generator"
    ) -> None:
        super().__init__(name)
        self.U_nom = nominal_voltage
        self.S_nom = nominal_power
        self.x = per_unit_reactance
        self.cos_phi = power_factor
        self.R_to_X_ratio = R_to_X_ratio
        self.c = voltage_factor
        # Zero-sequence scaling factor
        self.z0_factor = z0_factor

    @property
    def Z1(self) -> Quantity:
        U_nom = self.U_nom.to('V').magnitude
        S_nom = self.S_nom.to('VA').magnitude
        x = self.x.to('frac').magnitude
        Z_nom = U_nom ** 2 / S_nom
        X = x * Z_nom
        R = self.R_to_X_ratio * X
        Z_mag = math.sqrt(R ** 2 + X ** 2)
        sin_phi = math.sqrt(1 - self.cos_phi ** 2)
        Z_mag *= self.c / (1 + x * sin_phi)
        Z_ph = math.atan2(X, R)
        Z = cmath.rect(Z_mag, Z_ph)
        return Quantity(Z, 'ohm')

    @property
    def Z2(self) -> Quantity:
        """
        Negative-sequence impedance; approximated equal to Z1 for this
        component.
        """
        return self.Z1

    @property
    def Z0(self) -> Quantity:
        """
        Approximate zero-sequence impedance of a synchronous generator
        from its positive-sequence impedance:

            Z0 ≈ z0_factor * Z1

        Default:
            z0_factor = 1.0

        This reflects the common assumption that, for an effectively
        grounded neutral, the zero-sequence impedance is of the same
        order of magnitude as the positive-sequence impedance.

        Limitations:
            - Valid only if the generator neutral is brought out and
              effectively grounded (solidly or via a neutral impedance).
            - If the neutral is isolated or not accessible, zero-sequence
              currents cannot flow and the effective Z0 is very large.
            - Actual Z0 depends on detailed machine design and neutral
              earthing impedance; manufacturer data are preferable when
              available.
        """
        z1 = self.Z1.to('ohm').magnitude
        return Quantity(self.z0_factor * z1, 'ohm')


class SynchronousMotor(Generator):
    """
    Synchronous motor model reusing the Generator impedance formulation.

    The Z0 property inherited from Generator uses the same approximation:

        Z0 ≈ z0_factor * Z1

    Limitations are similar:
        - In many LV applications the stator neutral of synchronous motors
          is not brought out; in such cases the contribution to zero-sequence
          fault currents is negligible (Z0 effectively very large).
        - If a grounded neutral does exist, z0_factor can be adjusted to
          match manufacturer data or study assumptions.
    """
    pass


class Cable(Equipment):

    def __init__(
        self,
        length: Quantity,
        cross_section_area: Quantity,
        conductor_material: ConductorMaterials,
        cable_arrangement: CableArrangement,
        temperature: Quantity,
        z0_r_factor: float = 3.0,
        z0_x_factor: float = 3.0,
        name: str = "cable"
    ) -> None:
        super().__init__(name)
        self.L = length
        self.A = cross_section_area
        conductor_material = CONDUCTOR_MATERIALS[conductor_material]
        self.r = conductor_material.resistivity(temperature.to('degC').m)
        self.x = CABLE_UNIT_REACTANCES[cable_arrangement]
        # Zero-sequence scaling factors
        self.z0_r_factor = z0_r_factor
        self.z0_x_factor = z0_x_factor

    @property
    def Z1(self) -> Quantity:
        L = self.L.to('m').magnitude
        A = self.A.to('mm ** 2').magnitude
        R = self.r * L / A
        X = self.x * L
        Z_mag = math.sqrt(R ** 2 + X ** 2)
        Z_ph = math.atan2(X, R)
        Z = cmath.rect(Z_mag, Z_ph)
        return Quantity(Z, 'ohm')

    @property
    def Z2(self) -> Quantity:
        """
        Negative-sequence impedance; approximated equal to Z1 for this
        component.
        """
        return self.Z1

    @property
    def Z0(self) -> Quantity:
        """
        Approximate zero-sequence impedance of an LV cable from its
        positive-sequence impedance by scaling R1 and X1:

            R0 ≈ z0_r_factor * R1
            X0 ≈ z0_x_factor * X1

        Defaults:
            z0_r_factor = 3.0
            z0_x_factor = 3.0

        These values are a common rule of thumb for typical underground
        3-core or 4-core LV cables where the return path of zero-sequence
        currents includes soil and metallic sheaths or PEN conductors.

        Limitations:
            - Intended for typical LV distribution cables in ground.
            - Not valid for overhead lines, special sheath bonding
              arrangements, or unusual installation methods.
            - Actual Z0 depends strongly on soil resistivity, cable
              geometry, and sheath/PEN configuration.
            - When available, manufacturer or utility data for Z0 should
              be used instead of this approximation.
        """
        z1 = self.Z1.to('ohm').magnitude
        r1, x1 = z1.real, z1.imag
        r0 = self.z0_r_factor * r1
        x0 = self.z0_x_factor * x1
        return Quantity(complex(r0, x0), 'ohm')


class BusBar(Cable):

    def __init__(
        self,
        length: Quantity,
        cross_section_area: Quantity,
        conductor_material: ConductorMaterials,
        temperature: Quantity,
        z0_r_factor: float = 1.0,
        z0_x_factor: float = 1.0
    ) -> None:
        super().__init__(
            length, cross_section_area, conductor_material,
            CableArrangement.BUS_BAR, temperature,
            z0_r_factor=z0_r_factor,
            z0_x_factor=z0_x_factor
        )


class InductionMotor(Equipment):

    def __init__(
        self,
        nominal_voltage: Quantity,
        nominal_current: Quantity,
        locked_rotor_current: Quantity,
        P_m: Quantity,
        efficiency: Quantity,
        power_factor: float,
        R_to_X_ratio: float = 0.42,
        z2_factor: float = 1.0,
        z0_factor: float = 10.0,
        name: str = "induction_motor"
    ) -> None:
        super().__init__(name)
        self.U_nom = nominal_voltage
        self.I_nom = nominal_current
        self.I_start = locked_rotor_current
        self.P_m = P_m
        self.e = efficiency
        self.cos_phi = power_factor
        self.R_to_X_ratio = R_to_X_ratio
        # Scaling factor for the negative-sequence impedance
        self.z2_factor = z2_factor
        # Zero-sequence scaling factor (large ⇒ negligible contribution)
        self.z0_factor = z0_factor

    @property
    def Z1(self) -> Quantity:
        P_m = self.P_m.to('W').magnitude
        e = self.e.to('frac').magnitude
        I_nom = self.I_nom.to('A').magnitude
        I_start = self.I_start.to('A').magnitude
        U_nom = self.U_nom.to('V').magnitude
        S_nom = (P_m / e) / self.cos_phi
        Z_mag = (I_nom / I_start) * U_nom ** 2 / S_nom
        X = Z_mag / math.sqrt(1 + self.R_to_X_ratio ** 2)
        R = self.R_to_X_ratio * X
        Z_ph = math.atan2(X, R)
        Z = cmath.rect(Z_mag, Z_ph)
        return Quantity(Z, 'ohm')

    @property
    def Z2(self) -> Quantity:
        """
        Approximate negative-sequence impedance of the induction motor.

        In a detailed equivalent circuit of a cage induction motor, the
        negative-sequence impedance differs from the positive-sequence
        impedance because the effective slip for negative-sequence currents
        is (2 - s) instead of s. This affects the rotor branch impedance and
        therefore Z2.

        However, this simplified motor model does not include the separate
        rotor and stator parameters required to compute Z2 explicitly. For
        that reason, Z2 is here obtained by scaling the positive-sequence
        impedance:

            Z2 ≈ z2_factor * Z1

        Default:
            z2_factor = 1.0  →  Z2 ≈ Z1

        This is generally acceptable for design-level short-circuit studies,
        where the negative-sequence impedance of cage induction motors is
        often taken to be of the same order as the locked-rotor
        positive-sequence impedance.

        Limitations:
            - This is an approximation based on a lumped-impedance model.
            - It does not replace a full sequence model with explicit rotor
              and stator parameters.
            - If detailed manufacturer data (or a full equivalent circuit)
              are available, Z2 should preferably be computed from that
              data instead of using this scaling approach.
        """
        z1 = self.Z1.to('ohm').magnitude
        return Quantity(self.z2_factor * z1, 'ohm')

    @property
    def Z0(self) -> Quantity:
        """
        Very approximate 'zero-sequence impedance' of an induction motor.

        In most LV installations, three-phase induction motors do not have
        a neutral connection and therefore do not provide a significant path
        for zero-sequence currents. In network studies, their Z0 contribution
        is usually neglected, i.e. Z0 is treated as very large.

        This property returns:

            Z0 ≈ z0_factor * Z1

        with a large default z0_factor:

            z0_factor = 10.0

        so that the motor's contribution to zero-sequence fault currents is
        numerically small but still representable with a finite complex value.

        Limitations:
            - This is a modelling convenience rather than a physical formula.
            - If the motor is actually connected in a way that allows
              zero-sequence currents (e.g. special winding or grounded
              neutral), a dedicated sequence model should be used and
              z0_factor adjusted or replaced accordingly.
        """
        z1 = self.Z1.to('ohm').magnitude
        return Quantity(self.z0_factor * z1, 'ohm')
