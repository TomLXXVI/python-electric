from tkinter import ttk
from python_electric.gui import topology_window


if __name__ == '__main__':

    def main() -> None:
        root = topology_window.TopologyWindow()

        style = ttk.Style(root)
        style.configure("TButton", padding=10)

        root.mainloop()

    main()
