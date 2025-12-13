from scipy.interpolate import interp1d

from .. import Quantity
from ..utils.charts import LineChart

__all__ = [
    "SafetyCurve",
]


Q_ = Quantity


class SafetyCurve:
    """
    Implements the safety curves that indicate which combinations of fault
    voltage and fault duration are still acceptable to ensure that a person, in
    case of indirect contact with an exposed conductive part that has become
    energized due to an insulation fault, does not receive a dangerous current
    through the body.
    """
    t: list[float] = [5.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01]

    class AC:
        BB1: list[float] = [50.0, 72.0, 87.0, 207.0, 340.0, 465.0, 520.0, 543.0, 565.0]
        BB2: list[float] = [25.0, 43.0, 50.0, 109.0, 170.0, 227.0, 253.0, 263.0, 275.0]

    @classmethod
    def _validate_input(cls, voltage_type: str, skin_code: str) -> None:
        match voltage_type:
            case "AC":
                pass
            case _:
                raise ValueError(f"Unknown voltage type: {voltage_type}")
        match skin_code:
            case "BB1":
                pass
            case "BB2":
                pass
            case _:
                raise ValueError(f"Unknown skin code: {skin_code}")

    @classmethod
    def _select(cls, voltage_type: str, skin_code: str) -> tuple[list[float], list[float]]:
        if voltage_type == "AC":
            if skin_code == "BB1":
                return cls.t, cls.AC.BB1
            if skin_code == "BB2":
                return cls.t, cls.AC.BB2
        raise ValueError("Invalid input")

    @classmethod
    def plot_safety_curve(cls, voltage_type: str, skin_code: str) -> LineChart:
        cls._validate_input(voltage_type, skin_code)
        x = None
        if voltage_type == "AC":
            if skin_code == "BB1":
                x = cls.AC.BB1
            if skin_code == "BB2":
                x = cls.AC.BB2
        curve = LineChart()
        curve.add_xy_data(
            label="",
            x1_values=x,
            y1_values=cls.t
        )
        curve.x1.add_title("voltage, V")
        curve.y1.add_title("maximum contact time, s")
        return curve

    def __init__(
        self,
        voltage_type: str,
        skin_condition: str
    ) -> None:
        """
        Creates a `SafetyCurve` object that contains the safety curve for the
        given conditions.

        Parameters
        ----------
        voltage_type: str, {"AC"}
            Voltage type: either alternating current (AC) or direct current (DC).
            However, note that only AC has been implemented yet.
        skin_condition: str, {"BB1", "BB2"}
            Code that identifies the condition of the human skin: either dry
            (BB1) or wet (BB2).
        """
        self._validate_input(voltage_type, skin_condition)
        self.voltage_type = voltage_type
        self.skin_condition = skin_condition
        t_lst, UL_lst = self._select(voltage_type, skin_condition)
        self._fn_t = interp1d(
            UL_lst, t_lst,
            bounds_error=False,
            fill_value=(UL_lst[0], UL_lst[-1])
        )
        self._fn_UL = interp1d(
            t_lst, UL_lst,
            bounds_error=False,
            fill_value=(t_lst[0], t_lst[-1])
        )

    def max_contact_duration(self, U_f: Quantity) -> Quantity:
        """
        Returns the maximum allowable contact duration in seconds for the given
        fault voltage `U_f`.

        Returns
        -------
        Quantity
        """
        U_f = U_f.to('V').m
        t_c_max = self._fn_t(U_f)
        return Q_(t_c_max, 's')

    def max_touch_voltage(self, t_c: Quantity) -> Quantity:
        """
        Returns the conventional touch voltage limit in volts that goes with the
        given contact duration `t_c`.

        Returns
        -------
        Quantity
        """
        t_c = t_c.to('s').m
        U_touch_max = self._fn_UL(t_c)
        return Q_(U_touch_max, 'V')

    @property
    def conv_abs_limit_voltage(self) -> Quantity:
        return self.max_touch_voltage(t_c=Q_(5, 's'))
