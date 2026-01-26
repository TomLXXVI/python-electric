import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from .. import Quantity, Q_
from ..utils.lookup_table import LookupTable
from ..materials import ConductorMaterial, InsulationMaterial, Soil, SOIL_RESISTIVITY

__all__ = [
    "LoopEarthElectrode",
    "RodEarthElectrode",
    "PEConductor",
    "EarthConductor"
]

PI = math.pi


def create_k_factor_table_1():
    row_header = [
        ConductorMaterial.COPPER.value,
        ConductorMaterial.ALUMINIUM.value,
        ConductorMaterial.STEEL.value
    ]
    col_header = [
        InsulationMaterial.PVC.value,
        InsulationMaterial.EPR.value,
        InsulationMaterial.B.value
    ]
    data = [
        [143, 176, 160],
        [95,  116, 110],
        [52,   64,  60]
    ]
    description = "k-factors for insulated PE-conductors not part of a cable"
    cols_description = "insulation type"
    rows_description = "conductor material"
    lookup_table = LookupTable(
        data, col_header, row_header,
        description,
        cols_description,
        rows_description
    )
    return lookup_table


def create_k_factor_table_2():
    row_header = [
        ConductorMaterial.COPPER.value,
        ConductorMaterial.ALUMINIUM.value,
    ]
    col_header = [
        InsulationMaterial.PVC.value,
        InsulationMaterial.EPR.value,
        InsulationMaterial.B.value
    ]
    data = [
        [115, 143, 134],
        [76,   94,  89],
    ]
    description = "k-factors for PE-conductors part of a multicore cable"
    cols_description = "insulation type"
    rows_description = "conductor material"
    lookup_table = LookupTable(
        data, col_header, row_header,
        description,
        cols_description,
        rows_description
    )
    return lookup_table


tbl_k_factor_1 = create_k_factor_table_1()
tbl_k_factor_2 = create_k_factor_table_2()


class EarthElectrode(ABC):

    @abstractmethod
    def earthing_resistance(self, soil: Soil) -> Quantity:
        pass

@dataclass
class LoopEarthElectrode(EarthElectrode):
    wire_diameter: Quantity
    area: Quantity

    def __post_init__(self):
        self._d_wire = self.wire_diameter.to('m').m
        self._D_loop = self._equivalent_loop_diameter(self.area.to('m ** 2').m)

    @staticmethod
    def _equivalent_loop_diameter(A: float) -> float:
        return math.sqrt(4 * A / PI)

    def earthing_resistance(self, soil: Soil) -> Quantity:
        rho_soil = SOIL_RESISTIVITY[soil]
        R_e = rho_soil / (PI ** 2 * self._D_loop)
        R_e *= math.log(2 * PI * self._D_loop / self._d_wire)
        return Q_(R_e, 'ohm')


@dataclass
class RodEarthElectrode(EarthElectrode):
    diameter: Quantity
    length: Quantity

    def __post_init__(self):
        self._d_rod = self.diameter.to('m').m
        self._L = self.length.to('m').m

    def earthing_resistance(self, soil: Soil) -> Quantity:
        rho_soil = SOIL_RESISTIVITY[soil]
        R_e = rho_soil / (2 * PI * self._d_rod)
        R_e *= math.log(2 * PI * self._L / self._d_rod)
        return Q_(R_e, 'ohm')


@dataclass
class PEConductor:
    STD_SIZES = [
        1.5, 2.5, 4, 6, 10,
        16, 25, 35, 50, 70, 95,
        120, 150, 185, 240, 300
    ]  # mm²
    conductor_material: ConductorMaterial
    insulation_material: InsulationMaterial
    separated: bool

    k: float = field(init=False, default=0.0)
    S: Quantity = field(init=False, default=None)

    def __post_init__(self):
        if self.insulation_material == InsulationMaterial.PRC:
            self._insul_mat = InsulationMaterial.EPR
        else:
            self._insul_mat = self.insulation_material

        if self.separated:
            self.k = tbl_k_factor_1.data_value(
                self.conductor_material,
                self._insul_mat
            )
        else:
            self.k = tbl_k_factor_2.data_value(
                self.conductor_material,
                self._insul_mat
            )

    @classmethod
    def get_std_csa(cls, S: Quantity) -> Quantity:
        S_mag = S.to('mm**2').m
        delta_S = [abs(S_mag - S_std) for S_std in cls.STD_SIZES]
        delta_S_min = min(delta_S)
        i_min = delta_S.index(delta_S_min)
        S_std = cls.STD_SIZES[i_min]
        if S_std < S_mag:
            S_std = cls.STD_SIZES[min(i_min + 1, len(cls.STD_SIZES) - 1)]
        return Q_(S_std, 'mm**2')

    def cross_section_area(
        self,
        If: Quantity,
        t_interrupt: Quantity,
        mech_protected: bool = True
    ) -> Quantity:
        """
        Returns the standardized cross-sectional area for the PE-conductor.

        Parameters
        ----------
        If: Quantity
            Fault current that may flow through the PE-conductor.
        t_interrupt: Quantity
            Time before the protective device has completely interrupted the
            fault current.
        mech_protected: bool
            Indicates whether the PE-conductor is mechanically protected, e.g.
            PE-conductor in a conduit.

        Returns
        -------
        Quantity
        """
        If = If.to('A').m
        t_interrupt = t_interrupt.to('s').m
        if t_interrupt > 5.0:
            raise ValueError("Maximum interruption time is 5 s.")

        S_cal = If * math.sqrt(t_interrupt) / self.k
        S_std = self.get_std_csa(Q_(S_cal, 'mm**2'))
        if mech_protected:
            S = max(S_std.m, 2.5)
        else:
            S = max(S_std.m, 4.0)
        self.S = Q_(S, 'mm ** 2')
        return self.S


class EarthConductor(PEConductor):

    def cross_section_area(
        self,
        If: Quantity,
        t_interrupt: Quantity,
        mech_protected: bool = True
    ) -> Quantity:
        S = super().cross_section_area(If, t_interrupt, False)
        if mech_protected:
            S = max(S.m, 16.0)
        elif self.conductor_material == ConductorMaterial.COPPER:
            S = max(S.m, 25.0)
        else:
            S = max(S.m, 50.0)
        return Q_(S, 'mm ** 2')
