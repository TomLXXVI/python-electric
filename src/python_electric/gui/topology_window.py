from __future__ import annotations

import tkinter as tk
from tkinter import ttk, filedialog

from .. import Q_
from ..protection.earthing_system import EarthingSystem
from .models import ProjectModel
from .widgets import (
    LabeledEntry,
    LabeledCombobox,
    QuantityInputField,
    show_error,
    show_info
)
from .components_window import ComponentsWindow


class TopologyWindow(tk.Tk):
    """Window 1: configure source and connections, save/load project."""

    def __init__(self, project: ProjectModel | None = None):
        super().__init__()
        self.title("python-electric – Topology Builder (MVP)")
        self.geometry("980x560")

        self.project: ProjectModel = project or ProjectModel()
        self._path: str | None = None

        self._build_menu()
        self._build_layout()
        self._refresh_topology_tree()

    # ---------- UI building ----------
    def _build_menu(self) -> None:
        menubar = tk.Menu(self)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="New", command=self._on_new)
        file_menu.add_command(label="Open…", command=self._on_open)
        file_menu.add_command(label="Save", command=self._on_save)
        file_menu.add_command(label="Save As…", command=self._on_save_as)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        tools_menu = tk.Menu(menubar, tearoff=0)
        tools_menu.add_command(label="Components…", command=self._open_components_window)
        menubar.add_cascade(label="Tools", menu=tools_menu)

        self.config(menu=menubar)

    def _build_layout(self) -> None:
        root = ttk.Frame(self, padding=10)
        root.pack(fill="both", expand=True)

        # Main area: left column (project settings + forms) and right column (topology tree)
        main = ttk.Panedwindow(root, orient="horizontal")
        main.pack(fill="both", expand=True)

        left = ttk.Frame(main)
        right = ttk.Frame(main)
        main.add(left, weight=1)
        main.add(right, weight=2)

        # --- Left column: Project settings (top) ---
        proj = ttk.LabelFrame(left, text="Project settings", padding=10)
        proj.pack(fill="x")

        self.name_field = LabeledEntry(proj, text="name", label_width_px=170)
        self.name_field.set(self.project.name)
        self.name_field.grid(row=0, column=0, sticky="we", pady=4)

        self.un_field = QuantityInputField(
            proj, label="line-to-line voltage", units=["V", "kV"], default_unit="V", Q_=Q_, label_width_px=170
        )
        self.un_field.set_quantity(self.project.U_n)
        self.un_field.grid(row=1, column=0, sticky="we", pady=4)

        self.es_field = LabeledCombobox(
            proj,
            text="earthing system",
            values=[e.value for e in EarthingSystem],
            width=10,
            label_width_px=170,
        )
        self.es_field.set(self.project.earthing_system.value)
        self.es_field.grid(row=2, column=0, sticky="we", pady=4)

        self.neutral_var = tk.BooleanVar(value=self.project.neutral_distributed)
        ttk.Checkbutton(proj, text="neutral distributed (IT)", variable=self.neutral_var).grid(
            row=3, column=0, sticky="we", pady=4
        )

        ttk.Button(proj, text="Apply project settings", command=self._apply_project_settings).grid(
            row=4, column=0, columnspan=2, sticky="e", pady=8
        )

        # Let the first field expand in the project settings row
        proj.columnconfigure(0, weight=1)

        # --- Left column: Notebook (source + connections) under project settings ---
        self.nb = ttk.Notebook(left)
        self.nb.pack(fill="both", expand=True, pady=(10, 0))

        self.source_tab = ttk.Frame(self.nb, padding=8)
        self.conn_tab = ttk.Frame(self.nb, padding=8)
        self.nb.add(self.source_tab, text="Source")
        self.nb.add(self.conn_tab, text="Connections")

        self._build_source_tab(self.source_tab)
        self._build_connections_tab(self.conn_tab)

        # --- Right column: unified topology tree (source + connections) ---
        cols = ("type", "route")
        self.conn_tree = ttk.Treeview(
            right, columns=cols, show="tree headings", selectmode="browse", height=18,
        )
        self.conn_tree.heading("#0", text="ID")
        self.conn_tree.heading("type", text="Type")
        self.conn_tree.heading("route", text="Route")

        self.conn_tree.column("#0", width=120, stretch=False, anchor="w")
        self.conn_tree.column("type", width=110, stretch=False, anchor="w")
        self.conn_tree.column("route", width=260, stretch=True, anchor="w")

        self.conn_tree.pack(fill="both", expand=True, pady=(8, 0))
        self.conn_tree.bind("<<TreeviewSelect>>", self._on_tree_select)

        ttk.Button(right, text="Components…", command=self._open_components_window).pack(
            anchor="e", pady=(8, 0)
        )


    def _build_source_tab(self, tab: ttk.Frame) -> None:
        # Source form only (overview is shown in the topology tree on the right).
        left = ttk.Frame(tab, padding=8)
        left.pack(fill="both", expand=True)

        s = self.project.source

        self.src_conn_id = LabeledEntry(left, text="connection ID", label_width_px=170)
        self.src_conn_id.set(s.conn_id)
        self.src_conn_id.grid(row=0, column=0, sticky="we", pady=4)

        self.src_end_id = LabeledEntry(left, text="end node ID", label_width_px=170)
        self.src_end_id.set(s.end_id)
        self.src_end_id.grid(row=1, column=0, sticky="we", pady=4)

        self.src_u_l = QuantityInputField(
            left, label="line-to-line voltage", units=["V", "kV"], default_unit="V", Q_=Q_, label_width_px=170
        )
        self.src_u_l.set_quantity(s.U_l)
        self.src_u_l.grid(row=2, column=0, sticky="we", pady=4)

        self.src_s_sc = QuantityInputField(
            left, label="short-circuit power", units=["kVA", "MVA"], default_unit="MVA", Q_=Q_, label_width_px=170
        )
        self.src_s_sc.set_quantity(s.S_sc)
        self.src_s_sc.grid(row=3, column=0, sticky="we", pady=4)

        self.src_r_to_x = LabeledEntry(left, text="ratio R/X", label_width_px=170)
        self.src_r_to_x.set(str(s.R_to_X))
        self.src_r_to_x.grid(row=4, column=0, sticky="we", pady=4)

        self.src_z0_r = LabeledEntry(left, text="Z0-factor R", label_width_px=170)
        self.src_z0_r.set(str(s.z0_r_factor))
        self.src_z0_r.grid(row=5, column=0, sticky="we", pady=4)

        self.src_z0_x = LabeledEntry(left, text="Z0-factor X", label_width_px=170)
        self.src_z0_x.set(str(s.z0_x_factor))
        self.src_z0_x.grid(row=6, column=0, sticky="we", pady=4)

        ttk.Button(left, text="Apply source", command=self._apply_source).grid(
            row=7, column=0, sticky="e", pady=8
        )

        left.columnconfigure(0, weight=1)


    def _build_connections_tab(self, tab: ttk.Frame) -> None:
        # Connections form only (overview is shown in the topology tree on the right).
        left = ttk.Frame(tab, padding=8)
        left.pack(fill="both", expand=True)

        self.c_conn_id = LabeledEntry(left, text="connection ID", label_width_px=170)
        self.c_start_id = LabeledEntry(left, text="start node ID", label_width_px=170)
        self.c_end_id = LabeledEntry(left, text="end node ID", label_width_px=170)
        self.c_conn_id.grid(row=0, column=0, sticky="we", pady=4)
        self.c_start_id.grid(row=1, column=0, sticky="we", pady=4)
        self.c_end_id.grid(row=2, column=0, sticky="we", pady=4)

        btns = ttk.Frame(left)
        btns.grid(row=3, column=0, sticky="we", pady=(10, 0))
        ttk.Button(btns, text="Add", command=self._add_connection).pack(side="left")
        ttk.Button(btns, text="Update selected", command=self._update_selected_connection).pack(
            side="left", padx=6
        )
        ttk.Button(btns, text="Delete selected", command=self._delete_selected_connection).pack(side="left")

        left.columnconfigure(0, weight=1)


    # ---------- actions ----------
    def _apply_project_settings(self) -> None:
        try:
            self.project.name = self.name_field.get() or "network"
            self.project.U_n = self.un_field.get_quantity()
            self.project.earthing_system = EarthingSystem(self.es_field.get())
            self.project.neutral_distributed = bool(self.neutral_var.get())
        except Exception as e:
            show_error(self, "Invalid project settings", str(e))
            return
        show_info(self, "OK", "Project settings applied.")

    def _apply_source(self) -> None:
        try:
            s = self.project.source
            s.conn_id = self.src_conn_id.get() or "C0"
            s.end_id = self.src_end_id.get() or "B1"
            s.U_l = self.src_u_l.get_quantity()
            s.S_sc = self.src_s_sc.get_quantity()
            s.R_to_X = float(self.src_r_to_x.get().replace(",", ".") or "0.1")
            s.z0_r_factor = float(self.src_z0_r.get().replace(",", ".") or "0.0")
            s.z0_x_factor = float(self.src_z0_x.get().replace(",", ".") or "0.0")
        except Exception as e:
            show_error(self, "Invalid source data", str(e))
            return
        self._refresh_topology_tree()
        show_info(self, "OK", "Source updated.")

    def _add_connection(self) -> None:
        conn_id = self.c_conn_id.get()
        start_id = self.c_start_id.get()
        end_id = self.c_end_id.get()
        try:
            if not conn_id:
                raise ValueError("connection ID is required.")
            if not start_id or not end_id:
                raise ValueError("start and end node ID are required.")
            if start_id == end_id:
                raise ValueError("start and end node ID must differ.")
            self.project.add_connection(conn_id, start_id, end_id)
        except Exception as e:
            show_error(self, "Cannot add connection", str(e))
            return
        self._refresh_topology_tree()
        self._clear_conn_form()

    def _update_selected_connection(self) -> None:
        sel = self._selected_conn_id_from_tree()
        if not sel:
            show_error(self, "Update", "No connection selected.")
            return
        conn_id = self.c_conn_id.get()
        start_id = self.c_start_id.get()
        end_id = self.c_end_id.get()
        try:
            if not conn_id:
                raise ValueError("connection ID is required.")
            if start_id == end_id:
                raise ValueError("start and end node ID must differ.")
            # if conn_id changed, re-key dict
            conn = self.project.connections.pop(sel)
            conn.conn_id = conn_id
            conn.start_id = start_id
            conn.end_id = end_id
            self.project.connections[conn_id] = conn
        except Exception as e:
            show_error(self, "Cannot update", str(e))
            return
        self._refresh_topology_tree()

    def _delete_selected_connection(self) -> None:
        sel = self._selected_conn_id_from_tree()
        if not sel:
            return
        self.project.delete_connection(sel)
        self._refresh_topology_tree()
        self._clear_conn_form()

    def _clear_conn_form(self) -> None:
        self.c_conn_id.set("")
        self.c_start_id.set("")
        self.c_end_id.set("")

    
    def _selected_conn_id_from_tree(self) -> str | None:
        """Return selected connection id from the topology tree (excluding the source row)."""
        sel = self.conn_tree.selection()
        if not sel:
            return None
        iid = sel[0]
        if iid == "__SOURCE__":
            return None
        parent = self.conn_tree.parent(iid)
        if parent:
            iid = parent
        return iid if iid in self.project.connections else None

    def _refresh_topology_tree(self) -> None:
        """Refresh the unified topology tree (source + all connections)."""
        for item in self.conn_tree.get_children():
            self.conn_tree.delete(item)

        # Source row
        s = self.project.source
        src_id = s.conn_id or "C0"
        end_id = s.end_id or "B1"
        self.conn_tree.insert(
            "",
            "end",
            iid="__SOURCE__",
            text=src_id,
            values=("Source", f"GRID - {end_id}"),
        )

        # Connection rows
        for conn_id in self.project.list_connection_ids():
            c = self.project.connections[conn_id]
            self.conn_tree.insert(
                "",
                "end",
                iid=conn_id,
                text=c.conn_id,
                values=("Connection", f"{c.start_id} - {c.end_id}"),
            )

    def _on_tree_select(self, _evt=None) -> None:
        sel = self.conn_tree.selection()
        if not sel:
            return
        iid = sel[0]
        parent = self.conn_tree.parent(iid)
        if parent:
            iid = parent

        if iid == "__SOURCE__":
            self.nb.select(self.source_tab)
            s = self.project.source
            self.src_conn_id.set(s.conn_id)
            self.src_end_id.set(s.end_id)
            self.src_u_l.set_quantity(s.U_l)
            self.src_s_sc.set_quantity(s.S_sc)
            self.src_r_to_x.set(str(s.R_to_X))
            self.src_z0_r.set(str(s.z0_r_factor))
            self.src_z0_x.set(str(s.z0_x_factor))
            return

        if iid in self.project.connections:
            self.nb.select(self.conn_tab)
            c = self.project.connections[iid]
            self.c_conn_id.set(c.conn_id)
            self.c_start_id.set(c.start_id)
            self.c_end_id.set(c.end_id)

    # ---------- file actions ----------
    def _on_new(self) -> None:
        self.project = ProjectModel()
        self._path = None
        self._rebind_project_to_widgets()
        show_info(self, "New", "New project created.")

    def _on_open(self) -> None:
        path = filedialog.askopenfilename(
            parent=self,
            title="Open project",
            filetypes=[("python-electric project", "*.pxe"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            self.project = ProjectModel.load(path)
            self._path = path
            self._rebind_project_to_widgets()
        except Exception as e:
            show_error(self, "Open failed", str(e))

    def _on_save(self) -> None:
        if not self._path:
            self._on_save_as()
            return
        try:
            self._apply_project_settings()
            self._apply_source()
            self.project.save(self._path)
            show_info(self, "Saved", f"Saved to {self._path}")
        except Exception as e:
            show_error(self, "Save failed", str(e))

    def _on_save_as(self) -> None:
        path = filedialog.asksaveasfilename(
            parent=self,
            title="Save project as",
            defaultextension=".pxe",
            filetypes=[("python-electric project", "*.pxe")]
        )
        if not path:
            return
        self._path = path
        self._on_save()

    def _rebind_project_to_widgets(self) -> None:
        # project settings
        self.name_field.set(self.project.name)
        self.un_field.set_quantity(self.project.U_n)
        self.es_field.set(self.project.earthing_system.value)
        self.neutral_var.set(self.project.neutral_distributed)

        # source
        s = self.project.source
        self.src_conn_id.set(s.conn_id)
        self.src_end_id.set(s.end_id)
        self.src_u_l.set_quantity(s.U_l)
        self.src_s_sc.set_quantity(s.S_sc)
        self.src_r_to_x.set(str(s.R_to_X))
        self.src_z0_r.set(str(s.z0_r_factor))
        self.src_z0_x.set(str(s.z0_x_factor))
        self._refresh_topology_tree()

        self._refresh_topology_tree()

    # ---------- other windows ----------
    def _open_components_window(self) -> None:
        if not self.project.connections:
            show_error(self, "Components", "Add at least one connection first.")
            return
        ComponentsWindow(self, self.project)
