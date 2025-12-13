"""
Sizing a three-phase cable and then change the cross-sectional area of the
conductors.
"""
from python_electric import Quantity
from python_electric.materials import (
    ConductorMaterials,
    InsulationMaterials
)
from python_electric.sizing import (
    ThreePhaseCable,
    Ambient,
    InstallationMethods,
    CableMounting,
    CableArrangement
)

Q_ = Quantity


cable_W2 = ThreePhaseCable(
    P_mech=Q_(15, 'kW'),
    eta=0.84,
    cos_phi=0.81,
    U_line=Q_(400, 'V'),
    k_simul=1.0,
    k_ext=1.0,
    L=Q_(31, 'm'),
    conductor_material=ConductorMaterials.COPPER,
    insulation_material=InsulationMaterials.XLPE,
    ambient=Ambient.AIR,
    T_amb=Q_(30, 'degC'),
    install_method=InstallationMethods.E,
    cable_mounting=CableMounting.PERFORATED_TRAY,
    cable_arrangement=CableArrangement.MULTICORE,
    num_circuits=5,
    harmonic3_content=0.0,
    sizing_based_on_I_nom=False,
    name="W2"
)

print(f"{cable_W2.S:~P.0f}")
print(f"{cable_W2.I_z:~P.1f}")
print(f"{cable_W2.joule_integral.to('kA ** 2 * s'):~P.3f}")

cable_W2.change_csa_conductor(S=Q_(10, 'mm ** 2'))

print(f"\n{cable_W2.S:~P.0f}")
print(f"{cable_W2.I_z:~P.1f}")
print(f"{cable_W2.joule_integral.to('kA ** 2 * s'):~P.3f}")