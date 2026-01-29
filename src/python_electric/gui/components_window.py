from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from dataclasses import replace
from enum import StrEnum, IntEnum
from typing import Type

from .. import Q_
from ..protection.earthing_system import EarthingSystem
from ..network.topology import CableInput, TransformerInput
from ..network.components import Load, Transformer
from ..network.components.cable import (
    ConductorMaterial,
    InsulationMaterial,
    InstallMethod,
    CableMounting,
    CableArrangement,
    PhaseSystem,
)

from .models import ProjectModel
from .widgets import (
    LabeledEntry,
    LabeledCombobox,
    QuantityInputField,
    show_error,
    show_info,
)


class EnumMapper:

    def __init__(self, e: Type[StrEnum] | Type[IntEnum]):
        self._labels = {m: m.name.replace("_", " ").title() for m in e}
        self._label_to_enum = {label: m for m, label in self._labels.items()}
        self._enum_to_label = {m: label for m, label in self._labels.items()}

    def get_enum(self, label: str) -> StrEnum | IntEnum:
        return self._label_to_enum[label]

    def get_label(self, enum: StrEnum | IntEnum) -> str:
        return self._enum_to_label[enum]

    @property
    def labels(self) -> list[str]:
        return list(self._labels.values())


CableMountingMapper = EnumMapper(CableMounting)
CableArrangementMapper = EnumMapper(CableArrangement)
PhaseSystemMapper = EnumMapper(PhaseSystem)


# noinspection PyTypeChecker
class ComponentsWindow(tk.Toplevel):
    """
    Window 2: assign component input data to existing connections.
    """
    def __init__(self, master: tk.Misc, project: ProjectModel):
        super().__init__(master)
        self.title("python-electric – Components Builder (MVP)")
        self.geometry("980x560")
        self.project = project

        self.selected_conn_id: str | None = None
        self._suppress_tree_select: bool = False

        self._build_layout()
        self._populate_connection_list()

    def _build_layout(self) -> None:
        root = ttk.Frame(self, padding=10)
        root.pack(fill="both", expand=True)

        paned = ttk.Panedwindow(root, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left = ttk.Frame(paned, padding=8)
        right = ttk.Frame(paned, padding=8)
        paned.add(left, weight=1)
        paned.add(right, weight=3)

        # Left: connection tree (connections + assigned components)
        ttk.Label(left, text="Connections").pack(anchor="w")

        self.conn_tree = ttk.Treeview(
            left,
            columns=("route", "summary"),
            show="tree headings",
            height=20,
        )
        self.conn_tree.heading("#0", text="Connection")
        self.conn_tree.heading("route", text="Route")
        self.conn_tree.heading("summary", text="Components")
        self.conn_tree.column("#0", width=120, anchor="w")
        self.conn_tree.column("route", width=170, anchor="w")
        self.conn_tree.column("summary", width=180, anchor="w")
        self.conn_tree.pack(fill="both", expand=True, pady=(6, 0))
        self.conn_tree.bind("<<TreeviewSelect>>", self._on_select_connection_tree)

        # Right: notebook with component tabs
        header = ttk.Frame(right)
        header.pack(fill="x")

        self.conn_label = ttk.Label(header, text="Selected: —", font=("Segoe UI", 10, "bold"))
        self.conn_label.pack(side="left")

        self.nb = ttk.Notebook(right)
        self.nb.pack(fill="both", expand=True, pady=(10, 0))

        self.tab_cable = ttk.Frame(self.nb, padding=10)
        self.tab_transformer = ttk.Frame(self.nb, padding=10)
        self.tab_load = ttk.Frame(self.nb, padding=10)

        self.nb.add(self.tab_cable, text="Cable")
        self.nb.add(self.tab_transformer, text="Transformer")
        self.nb.add(self.tab_load, text="Load")

        self._build_cable_tab(self.tab_cable)
        self._build_transformer_tab(self.tab_transformer)
        self._build_load_tab(self.tab_load)

    # -------- connection selection --------
    def _populate_connection_list(self) -> None:
        self._refresh_connection_tree()

        # select first connection by default
        ids = self.project.list_connection_ids()
        if ids:
            self._select_conn_by_id(ids[0])

    def _refresh_connection_tree(self) -> None:
        """Refresh the tree view (connections + assigned components)."""
        self.conn_tree.delete(*self.conn_tree.get_children())

        for conn_id in self.project.list_connection_ids():
            c = self.project.connections[conn_id]
            route = f"{c.start_id} - {c.end_id}"

            # Build a small summary
            labels: list[str] = []
            if CableInput in c.components:
                labels.append("Cable")
            if TransformerInput in c.components:
                labels.append("Transformer")
            if c.load is not None:
                labels.append("Load")
            summary = ", ".join(labels) if labels else "—"

            # Parent node: the connection
            self.conn_tree.insert(
                "",
                tk.END,
                iid=conn_id,  # stable id for easy selection
                text=conn_id,
                values=(route, summary),
                open=True,
            )

            # Child nodes: assigned components
            # for comp_cls in c.components.keys():
            #     friendly = comp_cls.__name__.removesuffix("Input")
            #     self.conn_tree.insert(conn_id, tk.END, text=friendly, values=("", ""))
            #
            # if c.load is not None:
            #     self.conn_tree.insert(conn_id, tk.END, text="Load", values=("", ""))

    def _select_conn_by_id(self, conn_id: str, *, update_tree: bool = True) -> None:
        if conn_id not in self.project.connections:
            return

        # If already selected, do nothing (prevents pointless UI churn)
        if conn_id == getattr(self, "selected_conn_id", None):
            return

        self.selected_conn_id = conn_id
        c = self.project.connections[conn_id]
        self.conn_label.config(text=f"Selected: {c.conn_id}   ({c.start_id} - {c.end_id})")
        self._load_forms_for_connection(conn_id)

        if not update_tree:
            return

        # reflect selection in tree (programmatic) WITHOUT triggering recursion
        self._suppress_tree_select = True
        try:
            self.conn_tree.selection_set(conn_id)
            self.conn_tree.focus(conn_id)
            self.conn_tree.see(conn_id)
        except tk.TclError:
            pass
        finally:
            # reset after this event-loop tick
            self.after_idle(lambda: setattr(self, "_suppress_tree_select", False))

    def _on_select_connection_tree(self, _evt) -> None:
        if getattr(self, "_suppress_tree_select", False):
            return

        sel = self.conn_tree.selection()
        if not sel:
            return

        item = sel[0]
        parent = self.conn_tree.parent(item)
        conn_id = parent if parent else item

        # Here: user already selected it in the tree, so update_tree=False
        self._select_conn_by_id(conn_id, update_tree=False)

    # -------- forms: Cable --------
    def _build_cable_tab(self, tab: ttk.Frame) -> None:
        # minimal but useful set of fields
        self.cab_name = LabeledEntry(tab, text="name", label_width_px=170)
        self.cab_name.grid(row=0, column=0, sticky="we", pady=4)

        self.cab_L = QuantityInputField(tab, label="L", units=["m", "km"], default_unit="m", Q_=Q_, label_width_px=170)
        self.cab_L.grid(row=1, column=0, sticky="we", pady=4)

        self.cab_mat = LabeledCombobox(tab, text="conductor", values=[m.value for m in ConductorMaterial], label_width_px=170)
        self.cab_mat.grid(row=2, column=0, sticky="we", pady=4)

        self.cab_ins = LabeledCombobox(tab, text="insulation", values=[m.value for m in InsulationMaterial], label_width_px=170)
        self.cab_ins.grid(row=3, column=0, sticky="we", pady=4)

        self.cab_im = LabeledCombobox(tab, text="install method", values=[m.value for m in InstallMethod], width=12, label_width_px=170)
        self.cab_im.grid(row=4, column=0, sticky="we", pady=4)

        self.cab_mount = LabeledCombobox(tab, text="mounting", values=list(CableMountingMapper.labels), label_width_px=170, width=24)
        self.cab_mount.grid(row=5, column=0, sticky="we", pady=4)

        self.cab_arr = LabeledCombobox(tab, text="arrangement", values=list(CableArrangementMapper.labels), width=14, label_width_px=170)
        self.cab_arr.grid(row=6, column=0, sticky="we", pady=4)

        self.cab_phase = LabeledCombobox(tab, text="phase system", values=list(PhaseSystemMapper.labels), width=10, label_width_px=170)
        self.cab_phase.grid(row=7, column=0, sticky="we", pady=4)

        self.cab_es = LabeledCombobox(tab, text="earthing", values=[""] + [e.value for e in EarthingSystem], width=8, label_width_px=170)
        self.cab_es.grid(row=8, column=0, sticky="we", pady=4)

        btns = ttk.Frame(tab)
        btns.grid(row=9, column=0, sticky="we", pady=(14, 0))

        ttk.Button(btns, text="Apply to connection", command=self._apply_cable).pack(side="left")
        ttk.Button(btns, text="Reset to defaults", command=self._reset_cable_to_defaults).pack(side="left", padx=6)
        ttk.Button(btns, text="Save as defaults", command=self._save_cable_as_defaults).pack(side="left")

        tab.columnconfigure(0, weight=1)

    def _reset_cable_to_defaults(self) -> None:
        d: CableInput = self.project.get_defaults(CableInput)
        self._fill_cable_form(d)

    def _save_cable_as_defaults(self) -> None:
        try:
            cab = self._read_cable_form()
        except Exception as e:
            show_error(self, "Cable defaults", str(e))
            return
        self.project.set_defaults(cab)
        show_info(self, "Defaults", "Cable defaults updated.")

    def _apply_cable(self) -> None:
        if not self.selected_conn_id:
            return
        try:
            cab = self._read_cable_form()
            self.project.set_component(self.selected_conn_id, cab)
        except Exception as e:
            show_error(self, "Cable", str(e))
            return
        show_info(self, "Applied", f"Cable stored for connection {self.selected_conn_id}.")
        self._refresh_connection_tree()
        self._select_conn_by_id(self.selected_conn_id)

    def _fill_cable_form(self, cab: CableInput) -> None:
        self.cab_name.set(cab.name)
        self.cab_L.set_quantity(cab.L)
        self.cab_mat.set(cab.conductor_material.value)
        self.cab_ins.set(cab.insulation_material.value)
        self.cab_im.set(cab.install_method.value)
        self.cab_mount.set(CableMountingMapper.get_label(cab.cable_mounting))
        self.cab_arr.set(CableArrangementMapper.get_label(cab.cable_arrangement))
        self.cab_phase.set(PhaseSystemMapper.get_label(cab.phase_system))
        self.cab_es.set("" if cab.earthing_system is None else cab.earthing_system.value)

    def _read_cable_form(self) -> CableInput:
        d: CableInput = self.project.get_defaults(CableInput)
        name = self.cab_name.get() or d.name
        L = self.cab_L.get_quantity()
        es_raw = self.cab_es.get()
        earthing = EarthingSystem(es_raw) if es_raw else None
        return replace(
            d,
            name=name,
            L=L,
            conductor_material=ConductorMaterial(self.cab_mat.get()),
            insulation_material=InsulationMaterial(self.cab_ins.get()),
            install_method=InstallMethod(self.cab_im.get()),
            cable_mounting=CableMountingMapper.get_enum(self.cab_mount.get()),
            cable_arrangement=CableArrangementMapper.get_enum(self.cab_arr.get()),
            phase_system=PhaseSystemMapper.get_enum(self.cab_phase.get()),
            earthing_system=earthing,
        )

    # -------- forms: Transformer --------
    def _build_transformer_tab(self, tab: ttk.Frame) -> None:
        self.tr_name = LabeledEntry(tab, text="name", label_width_px=170)
        self.tr_name.grid(row=0, column=0, sticky="we", pady=4)

        self.tr_Sn = QuantityInputField(tab, label="nominal power", units=["kVA", "MVA"], default_unit="kVA", Q_=Q_, label_width_px=170)
        self.tr_Sn.grid(row=1, column=0, sticky="we", pady=4)

        self.tr_Ulp = QuantityInputField(tab, label="primary voltage", units=["V", "kV"], default_unit="kV", Q_=Q_, label_width_px=170)
        self.tr_Ulp.grid(row=2, column=0, sticky="we", pady=4)

        self.tr_Uls = QuantityInputField(tab, label="secondary voltage", units=["V", "kV"], default_unit="V", Q_=Q_, label_width_px=170)
        self.tr_Uls.grid(row=3, column=0, sticky="we", pady=4)

        self.tr_ucc = QuantityInputField(tab, label="percent impedance voltage", units=["pct"], default_unit="pct", Q_=Q_, label_width_px=170)
        self.tr_ucc.grid(row=4, column=0, sticky="we", pady=4)

        self.tr_Pcu = QuantityInputField(tab, label="copper loss", units=["W", "kW"], default_unit="kW", Q_=Q_, label_width_px=170)
        self.tr_Pcu.grid(row=5, column=0, sticky="we", pady=4)

        wvals = [w.value for w in Transformer.WindingConn]
        self.tr_pri = LabeledCombobox(tab, text="primary connection", values=wvals, width=10, label_width_px=170)
        self.tr_sec = LabeledCombobox(tab, text="secondary connection", values=wvals, width=10, label_width_px=170)
        self.tr_pri.grid(row=6, column=0, sticky="we", pady=4)
        self.tr_sec.grid(row=7, column=0, sticky="we", pady=4)

        btns = ttk.Frame(tab)
        btns.grid(row=8, column=0, sticky="we", pady=(14, 0))
        ttk.Button(btns, text="Apply to connection", command=self._apply_transformer).pack(side="left")
        ttk.Button(btns, text="Reset to defaults", command=self._reset_transformer_to_defaults).pack(side="left", padx=6)
        ttk.Button(btns, text="Save as defaults", command=self._save_transformer_as_defaults).pack(side="left")

        tab.columnconfigure(0, weight=1)

    def _reset_transformer_to_defaults(self) -> None:
        d: TransformerInput = self.project.get_defaults(TransformerInput)
        self._fill_transformer_form(d)

    def _save_transformer_as_defaults(self) -> None:
        try:
            tr = self._read_transformer_form()
        except Exception as e:
            show_error(self, "Transformer defaults", str(e))
            return
        self.project.set_defaults(tr)
        show_info(self, "Defaults", "Transformer defaults updated.")

    def _apply_transformer(self) -> None:
        if not self.selected_conn_id:
            return
        try:
            tr = self._read_transformer_form()
            self.project.set_component(self.selected_conn_id, tr)
        except Exception as e:
            show_error(self, "Transformer", str(e))
            return
        show_info(self, "Applied", f"Transformer stored for connection {self.selected_conn_id}.")
        self._refresh_connection_tree()
        self._select_conn_by_id(self.selected_conn_id)

    def _fill_transformer_form(self, tr: TransformerInput) -> None:
        self.tr_name.set(tr.name)
        self.tr_Sn.set_quantity(tr.S_n)
        self.tr_Ulp.set_quantity(tr.U_lp)
        self.tr_Uls.set_quantity(tr.U_ls)
        self.tr_ucc.set_quantity(tr.u_cc)
        self.tr_Pcu.set_quantity(tr.P_Cu)
        self.tr_pri.set("" if tr.pri_conn is None else tr.pri_conn.value)
        self.tr_sec.set("" if tr.sec_conn is None else tr.sec_conn.value)

    def _read_transformer_form(self) -> TransformerInput:
        d: TransformerInput = self.project.get_defaults(TransformerInput)
        name = self.tr_name.get() or d.name
        pri_raw = self.tr_pri.get()
        sec_raw = self.tr_sec.get()
        pri = Transformer.WindingConn(pri_raw) if pri_raw else None
        sec = Transformer.WindingConn(sec_raw) if sec_raw else None
        return replace(
            d,
            name=name,
            S_n=self.tr_Sn.get_quantity(),
            U_lp=self.tr_Ulp.get_quantity(),
            U_ls=self.tr_Uls.get_quantity(),
            u_cc=self.tr_ucc.get_quantity(),
            P_Cu=self.tr_Pcu.get_quantity(),
            pri_conn=pri,
            sec_conn=sec,
        )

    # -------- forms: Load --------
    def _build_load_tab(self, tab: ttk.Frame) -> None:
        self.ld_Pe = QuantityInputField(tab, label="electrical power", units=["W", "kW"], default_unit="kW", Q_=Q_, label_width_px=170)
        self.ld_Pe.grid(row=0, column=0, sticky="we", pady=4)

        self.ld_cos = LabeledEntry(tab, text="cos phi", label_width_px=170)
        self.ld_cos.set("0.8")
        self.ld_cos.grid(row=1, column=0, sticky="we", pady=4)

        self.ld_Ul = QuantityInputField(tab, label="line-to-line voltage (opt.)", units=["V", "kV"], default_unit="V", Q_=Q_, label_width_px=170)
        self.ld_Ul.grid(row=2, column=0, sticky="we", pady=4)

        btns = ttk.Frame(tab)
        btns.grid(row=3, column=0, sticky="we", pady=(14, 0))
        ttk.Button(btns, text="Apply to connection", command=self._apply_load).pack(side="left")
        ttk.Button(btns, text="Clear load", command=self._clear_load).pack(side="left", padx=6)

        tab.columnconfigure(0, weight=1)

    def _apply_load(self) -> None:
        if not self.selected_conn_id:
            return
        try:
            Pe = self.ld_Pe.get_quantity()
            cos_phi = float(self.ld_cos.get().replace(",", "."))
            if not (0 < cos_phi <= 1.0):
                raise ValueError("cos phi must be in (0, 1].")
            # U_l optional: if empty, keep None so NetworkTopology sets it to U_n
            Ul_raw = self.ld_Ul.value_var.get().strip()
            Ul = None if Ul_raw == "" else self.ld_Ul.get_quantity()
            self.project.connections[self.selected_conn_id].load = Load(U_l=Ul, cos_phi=cos_phi, P_e=Pe)
        except Exception as e:
            show_error(self, "Load", str(e))
            return
        show_info(self, "Applied", f"Load stored for connection {self.selected_conn_id}.")
        self._refresh_connection_tree()
        self._select_conn_by_id(self.selected_conn_id)

    def _clear_load(self) -> None:
        if not self.selected_conn_id:
            return
        self.project.connections[self.selected_conn_id].load = None
        show_info(self, "Cleared", f"Load removed from connection {self.selected_conn_id}.")
        self._refresh_connection_tree()
        self._select_conn_by_id(self.selected_conn_id)

    # -------- load forms on selection --------
    def _load_forms_for_connection(self, conn_id: str) -> None:
        # Cable
        cab = self.project.get_component(conn_id, CableInput)
        if cab is None:
            cab = self.project.get_defaults(CableInput)
        self._fill_cable_form(cab)

        # Transformer
        tr = self.project.get_component(conn_id, TransformerInput)
        if tr is None:
            tr = self.project.get_defaults(TransformerInput)
        self._fill_transformer_form(tr)

        # Load (stored on connection)
        conn = self.project.connections[conn_id]
        if conn.load is None:
            # clear fields
            self.ld_Pe.set_quantity(Q_(0, "W"))
            self.ld_cos.set("0.8")
            self.ld_Ul.set_quantity(None)
        else:
            self.ld_Pe.set_quantity(conn.load.P_e)
            self.ld_cos.set(str(conn.load.cos_phi))
            self.ld_Ul.set_quantity(conn.load.U_l)
