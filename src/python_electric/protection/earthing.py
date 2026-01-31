import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from .. import Quantity, Q_
from ..utils.lookup_table import LookupTable
from ..materials import ConductorMaterial, InsulationMaterial, Soil, SOIL_RESISTIVITY
from ..sizing import get_stdcsa


__all__ = [
    "EarthLoopElectrode",
    "EarthRodElectrode",
    "PEConductor",
    "EarthingConductor"
]


PI = math.pi

# ------------------------------------------------------------------------------
# Look-up table k-factor - PE-conductor separated.

def create_k_table_separated():
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
    data = [[143, 176, 160], [95,  116, 110], [52,   64,  60]]
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


# Look-up table k-factor - PE-conductor part of multicore table
def create_k_table_multicore():
    row_header = [
        ConductorMaterial.COPPER.value,
        ConductorMaterial.ALUMINIUM.value,
    ]
    col_header = [
        InsulationMaterial.PVC.value,
        InsulationMaterial.EPR.value,
        InsulationMaterial.B.value
    ]
    data = [[115, 143, 134], [76,   94,  89]]
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


ktbl_separated = create_k_table_separated()
ktbl_multicore = create_k_table_multicore()


class EarthElectrode(ABC):

    @abstractmethod
    def earthing_resistance(self, soil: Soil) -> Quantity:
        pass

@dataclass
class EarthLoopElectrode(EarthElectrode):
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
class EarthRodElectrode(EarthElectrode):
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
    separated: bool
    S_ph: Quantity | None = None

    k: float = field(init=False, default=0.0)
    S_pe: Quantity = field(init=False, default=None)

    def __post_init__(self):
        if self.insulation_material == InsulationMaterial.PRC:
            self.insulation_material = InsulationMaterial.EPR

        # Get k-factor.
        if self.separated:
            self.k = ktbl_separated.data_value(
                self.conductor_material,
                self.insulation_material
            )
        else:
            self.k = ktbl_multicore.data_value(
                self.conductor_material,
                self.insulation_material
            )

    def csa_adiabatic(self, If: Quantity, ti: Quantity) -> Quantity:
        """
        Calculates the standardized csa of the PE-conductor assuming adiabatic
        heating of the conductor during the fault.
        """
        If = If.to('A').m
        ti = min(ti.to('s').m, 5.0)
        S_pe = get_stdcsa(Q_(If * math.sqrt(ti) / self.k, 'mm**2'))
        return S_pe

    def csa(
        self,
        If: Quantity,
        ti: Quantity,
        mech_prot: bool = True
    ) -> Quantity:
        """
        Returns the standardized cross-sectional area for the PE-conductor.

        Parameters
        ----------
        If: Quantity
            Fault current that may flow through the PE-conductor.
        ti: Quantity
            Time the protective device needs to interrupt the fault current
            completely. This time is limited to 5 s.
        mech_prot: bool
            Indicates whether the PE-conductor is mechanically protected, e.g.
            a PE-conductor in a conduit.

        Returns
        -------
        Quantity
        """
        # First calculate the csa assuming adiabatic heating during the fault
        S_pe_cal = self.csa_adiabatic(If, ti)

        # Rule-based determination of S_pe
        if self.separated:
            if Q_(0.0, 'mm**2') < self.S_ph <= Q_(16, 'mm ** 2'):
                S_pe = Q_(max(self.S_ph.m, 2.5 if mech_prot else 4.0), 'mm**2')
            elif Q_(16, 'mm ** 2') < self.S_ph <= Q_(35, 'mm ** 2'):
                S_pe = Q_(16, 'mm ** 2')
            else:   # self.S_ph > Q_(35, 'mm ** 2'):
                S_pe = get_stdcsa(self.S_ph / 2)
            self.S_pe = max(S_pe_cal, S_pe)
        else:
            # PE = part of multicore cable -> same as S_ph
            self.S_pe = self.S_ph
        return self.S_pe


class EarthingConductor(PEConductor):

    def csa(
        self,
        If: Quantity,
        ti: Quantity,
        mech_prot: bool = True
    ) -> Quantity:
        """
        Returns the standardized cross-sectional area for the earthing
        conductor that connects with the earth electrode.
        """
        S_pe = self.csa_adiabatic(If, ti)
        if mech_prot:
            S_pe = max(S_pe, Q_(16.0, 'mm**2'))
        elif self.conductor_material == ConductorMaterial.COPPER:
            S_pe = max(S_pe, Q_(25.0, 'mm**2'))
        else:
            S_pe = max(S_pe, Q_(50.0, 'mm**2'))
        self.S_pe = S_pe
        return self.S_pe
