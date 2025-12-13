"""
Sizing a three-phase cable with the high-level, user API and selecting a circuit
breaker for this cable.
"""
from python_electric import Quantity, VoltRef
from python_electric.materials import ConductorMaterials, InsulationMaterials
from python_electric.sizing import (
    ThreePhaseCable,
    InstallationMethods,
    CableMounting,
    CableArrangement,
    Ambient,
)
from python_electric.protection import CircuitBreaker

Q_ = Quantity


cable = ThreePhaseCable(
    I_load=Q_(100, 'A'),
    L=Q_(20, 'm'),
    conductor_material=ConductorMaterials.COPPER,
    insulation_material=InsulationMaterials.PVC,
    ambient=Ambient.AIR,
    T_amb=Q_(40, 'degC'),
    install_method=InstallationMethods.E,
    cable_mounting=CableMounting.LADDER,
    cable_arrangement=CableArrangement.MULTICORE,
    num_circuits=6,
    harmonic3_content=0.0
)

print(
    f"installation reduction factor = {cable.sizing_data.sizing_factor:.2f}",  # from low-level sizing API
    f"nominal current = {cable.I_nom:~P.1f}",
    f"standard current-carrying capacity = {cable.sizing_data.current_capacity_ref:.2f} A",  # from low-level sizing API
    f"cross-sectional area = {cable.S:~P.0f}",
    f"actual current-carrying capacity = {cable.I_z:~P.2f}",
    sep="\n"
)

cable.voltage_drop(
    U_line=Q_(230.0, 'V'),
    I_load=cable.I_load,
    volt_ref=VoltRef.PH3_GROUND_TO_LINE,
    cos_phi=0.8
)
print(f"voltage drop = {cable.U_drop:~P.2f} ({cable.U_drop_pct:~P.2f})")


# -> from short-circuit calculations:
Ik_max = Q_(11.2, 'kA').to('A')
Ik_min = Q_(8.76, 'kA').to('A')


cb = CircuitBreaker(
    standard=CircuitBreaker.Standard.INDUSTRIAL,
    category=CircuitBreaker.Category.ADJUSTABLE,
    load_current=cable.I_load,
    nom_current=cable.I_nom,
    ampacity=cable.I_z,
    joule_integral=cable.joule_integral,
    ultim_break_capacity=Q_(50, 'kA'),
    k_magn_trip=2.0
)
print(
    f"cable overload protection: "
    f"{cb.check_overload_protection()}"
)
print(
    f"cable short-circuit protection: "
    f"{cb.check_shortcircuit_protection(Ik_max, Ik_min)}"
)
print(
    f"conventional non-tripping current: {cb.I_nf}",
    f"conventional tripping current: {cb.I_f}",
    f"conventional tripping time: {cb.t_conv}",
    f"minimum magnetic trip current: {cb.I_m_min}",
    f"maximum magnetic trip current: {cb.I_m_max}",
    f"magnetic trip time limit: {cb.t_m_lim}",
    sep="\n"
)
