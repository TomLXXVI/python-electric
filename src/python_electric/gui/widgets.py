from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from dataclasses import dataclass
from typing import Callable, Sequence, Any

from .. import Quantity


@dataclass(frozen=True)
class QuantitySpec:
    """Describes a quantity widget configuration."""
    units: Sequence[str]
    default_unit: str


class LabeledEntry(ttk.Frame):
    """Label + Entry convenience widget."""

    def __init__(self, master, *, text: str, width: int = 18, label_width_px: int = 140, **kwargs) -> None:
        super().__init__(master, **kwargs)
        self.columnconfigure(0, minsize=label_width_px)
        self.columnconfigure(1, weight=1)
        self.var = tk.StringVar(value="")
        ttk.Label(self, text=text).grid(row=0, column=0, sticky="w", padx=(0, 6))
        self.entry = ttk.Entry(self, textvariable=self.var, width=width)
        self.entry.grid(row=0, column=1, sticky="we")

    def get(self) -> str:
        return self.var.get().strip()

    def set(self, value: str) -> None:
        self.var.set("" if value is None else str(value))


class LabeledSpinbox(ttk.Frame):
    """Label + Spinbox convenience widget."""

    def __init__(self, master, *, text: str, from_: int, to: int, width: int = 8, label_width_px: int = 140, **kwargs):
        super().__init__(master, **kwargs)
        self.columnconfigure(0, minsize=label_width_px)
        self.columnconfigure(1, weight=1)
        self.var = tk.StringVar(value=str(from_))
        ttk.Label(self, text=text).grid(row=0, column=0, sticky="w", padx=(0, 6))
        self.spin = ttk.Spinbox(self, textvariable=self.var, from_=from_, to=to, width=width)
        self.spin.grid(row=0, column=1, sticky="w")

    def get_int(self) -> int:
        raw = self.var.get().strip()
        return int(raw)


class LabeledCombobox(ttk.Frame):
    """Label + readonly Combobox convenience widget."""

    def __init__(self, master, *, text: str, values: Sequence[str], width: int = 16, label_width_px: int = 140, **kwargs):
        super().__init__(master, **kwargs)
        self.columnconfigure(0, minsize=label_width_px)
        self.columnconfigure(1, weight=1)
        self.var = tk.StringVar(value=values[0] if values else "")
        ttk.Label(self, text=text).grid(row=0, column=0, sticky="w", padx=(0, 6))
        self.combo = ttk.Combobox(self, textvariable=self.var, values=list(values), state="readonly", width=width)
        self.combo.grid(row=0, column=1, sticky="we")

    def get(self) -> str:
        return self.var.get()

    def set(self, value: str) -> None:
        self.var.set(value)


class QuantityInputField(ttk.Frame):
    """Numeric entry + unit combobox that returns a pint Quantity via injected Q_."""

    def __init__(
        self,
        master,
        *,
        label: str,
        units: Sequence[str],
        default_unit: str | None = None,
        width_value: int = 10,
        width_unit: int = 7,
        label_width_px: int = 140,
        allow_negative: bool = False,
        Q_: Callable[[float, str], Any] | None = None,
        **kwargs
    ):
        super().__init__(master, **kwargs)
        self.columnconfigure(0, minsize=label_width_px)
        self.columnconfigure(1, weight=1)
        if not units:
            raise ValueError("units must be non-empty.")
        self._units = list(units)
        self._Q_ = Q_
        self._allow_negative = allow_negative

        self.value_var = tk.StringVar(value="")
        self.unit_var = tk.StringVar(value=default_unit or self._units[0])

        ttk.Label(self, text=label).grid(row=0, column=0, sticky="w", padx=(0, 6))
        self.entry = ttk.Entry(self, textvariable=self.value_var, width=width_value)
        self.entry.grid(row=0, column=1, sticky="we")
        self.combo = ttk.Combobox(
            self, textvariable=self.unit_var, values=self._units, state="readonly", width=width_unit
        )
        self.combo.grid(row=0, column=2, sticky="w", padx=(6, 0))

        self.entry.bind("<FocusIn>", lambda e: self.entry.selection_range(0, tk.END))

    def set_quantity(self, q: Any) -> None:
        """Populate widget from an existing Quantity (tries to convert to one of the supported units)."""
        if q is None:
            self.value_var.set("")
            self.unit_var.set(self._units[0])
            return
        elif isinstance(q, Quantity):
            self.value_var.set(f"{q.magnitude:g}")
            self.unit_var.set(f"{q.units:~P}")
            # for u in self._units:
            #     try:
            #         q_u = q.to(u)
            #     except Exception:
            #         continue
            #     self.value_var.set(f"{q_u.magnitude:g}")
            #     self.unit_var.set(u)
            return
        else:
            # fallback
            self.value_var.set(f"{getattr(q, 'magnitude', q)}")
            self.unit_var.set(self._units[0])

    def get_quantity(self) -> Any:
        """Return Quantity using Q_(value, unit). Raises ValueError on invalid input."""
        if self._Q_ is None:
            raise RuntimeError("Q_ was not injected into QuantityInputField.")
        raw = self.value_var.get().strip().replace(",", ".")
        if raw == "":
            raise ValueError("Value is empty.")
        try:
            value = float(raw)
        except ValueError as e:
            raise ValueError(f"Invalid number: {raw}") from e
        if (not self._allow_negative) and value < 0:
            raise ValueError("Negative values are not allowed.")
        unit = self.unit_var.get()
        if unit not in self._units:
            raise ValueError(f"Unsupported unit: {unit}")
        return self._Q_(value, unit)


def show_error(master: tk.Misc, title: str, message: str) -> None:
    from tkinter import messagebox
    parent = master.winfo_toplevel()
    parent.update()
    parent.lift()
    parent.focus_force()
    messagebox.showerror(title, message, parent=parent)


def show_info(master: tk.Misc, title: str, message: str) -> None:
    from tkinter import messagebox
    parent = master.winfo_toplevel()
    parent.update()
    parent.lift()
    parent.focus_force()
    messagebox.showinfo(title, message, parent=parent)
