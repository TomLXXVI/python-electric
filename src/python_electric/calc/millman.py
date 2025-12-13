__all__ = ["MillmanTheorem"]


class MillmanTheorem:

    def __init__(
        self,
        phase_voltages: complex | tuple[complex, ...],
        source_impedances: complex | tuple[complex, ...] | None,
        load_impedances: complex | tuple[complex, ...],
        neutral_impedance: complex
    ) -> None:
        """
        Initializes a `MillmanTheorem` calculation object.

        Parameters
        ----------
        phase_voltages: tuple[complex, ...]
            Phase (line-to-ground) voltages (E1, E2, E3) of the three-phase 
            system in volts.
        source_impedances: tuple[complex, ...]
            Internal phase impedances (Z_s1, Z_s2, Z_s3) of the three-phase 
            voltage source in ohms.
        load_impedances: tuple[complex, ...]
            Phase impedances of the load (Z_l1, Z_l2, Z_l3) in ohms.
        """
        if not isinstance(phase_voltages, tuple):
            phase_voltages = tuple([phase_voltages for _ in range(3)])
        if source_impedances is None:
            source_impedances = tuple([0.0 for _ in range(3)])
        elif not isinstance(source_impedances, tuple):
            source_impedances = tuple([source_impedances for _ in range(3)])
        if not isinstance(load_impedances, tuple):
            load_impedances = tuple([load_impedances for _ in range(3)])

        E_1 = phase_voltages[0]
        E_2 = phase_voltages[1]
        E_3 = phase_voltages[2]
        Z_s1 = source_impedances[0]
        Z_s2 = source_impedances[1]
        Z_s3 = source_impedances[2]
        Z_l1 = load_impedances[0]
        Z_l2 = load_impedances[1]
        Z_l3 = load_impedances[2]
        Z_0 = neutral_impedance
        Y_1 = 1 / (Z_s1 + Z_l1)
        Y_2 = 1 / (Z_s2 + Z_l2)
        Y_3 = 1 / (Z_s3 + Z_l3)
        Y_0 = 1 / Z_0

        n = Y_1 * E_1 + Y_2 * E_2 + Y_3 * E_3
        d = Y_1 + Y_2 + Y_3 + Y_0
        self._U_nn = n / d
        self._I = [
            Y_i * (E_i - self._U_nn)
            for Y_i, E_i in zip((Y_1, Y_2, Y_3), (E_1, E_2, E_3))
        ]
        self._I_n = sum(self._I)

    @property
    def deltaU_neutral(self) -> complex:
        """
        Returns the voltage difference of the load-side neutral with respect to
        the source-side neutral in volts.

        Returns
        -------
        complex
        """
        return self._U_nn

    @property
    def I_line(self) -> tuple[complex, ...]:
        """
        Returns the line currents in amperes.

        Returns
        -------
        tuple[complex, ...]
        """
        return tuple(self._I)

    @property
    def I_n(self) -> complex:
        """
        Returns the neutral current in amperes.

        Returns
        -------
        complex
        """
        return self._I_n
