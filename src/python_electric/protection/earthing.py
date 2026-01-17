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
    conductor_material: ConductorMaterial
    insulation_material: InsulationMaterial
    seperate: bool

    k: float = field(init=False, default=0.0)
    S: Quantity = field(init=False, default=None)

    def __post_init__(self):
        if self.insulation_material == InsulationMaterial.PRC:
            self._insulation_material = InsulationMaterial.EPR
        else:
            self._insulation_material = self.insulation_material

        if self.seperate:
            self.k = tbl_k_factor_1.data_value(
                self.conductor_material,
                self._insulation_material
            )
        else:
            self.k = tbl_k_factor_2.data_value(
                self.conductor_material,
                self._insulation_material
            )

    def cross_section_area(
        self,
        I_f: Quantity,
        t_u: Quantity,
        protected: bool = True
    ) -> Quantity:
        std_sizes = [1.5, 2.5, 4, 6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240, 300]
        I_f = I_f.to('A').m
        t_u = t_u.to('s').m
        if t_u > 5.0:
            raise ValueError("Maximum interruption time is 5 s.")
        S = I_f * math.sqrt(t_u) / self.k
        delta_S = [abs(S - S_std) for S_std in std_sizes]
        delta_S_min = min(delta_S)
        i_min = delta_S.index(delta_S_min)
        S = std_sizes[i_min]
        if protected:
            S = max(S, 2.5)
        else:
            S = max(S, 4.0)
        self.S = Q_(S, 'mm ** 2')
        return self.S


class EarthConductor(PEConductor):

    def cross_section_area(
        self,
        I_f: Quantity,
        t_u: Quantity,
        protected: bool = True
    ) -> Quantity:
        S = super().cross_section_area(I_f, t_u, False)
        if protected:
            S = max(S.m, 16.0)
        elif self.conductor_material == ConductorMaterial.COPPER:
            S = max(S.m, 25.0)
        else:
            S = max(S.m, 50.0)
        return Q_(S, 'mm ** 2')
