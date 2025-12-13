"""
Sizing a cable with the low-level API (without using Quantity)
"""
from python_electric.materials import (
    ConductorMaterials,
    InsulationMaterials
)
from python_electric.sizing.cable.sizing import (
    CableSizingData,
    Ambient,
    InstallationMethods,
    CableMounting,
    CableArrangement,
    get_cable_sizing,
    set_cable_sizing
)


# Collect cable data.
cable_data = CableSizingData(
    rated_load_current=31.8,
    length=31.0,
    conductor_material=ConductorMaterials.COPPER,
    insulation_material=InsulationMaterials.XLPE,
    ambient=Ambient.AIR,
    amb_temperature=30.0,
    install_method=InstallationMethods.E,
    cable_mounting=CableMounting.PERFORATED_TRAY,
    cable_arrangement=CableArrangement.MULTICORE,
    num_circuits=5
)

# Run sizing routine...
cable_data = get_cable_sizing(cable_data)
print(cable_data.cross_section_area)
print(cable_data.current_capacity)
print(cable_data.joule_integral)
print()


# Set conductor cross-sectional area of the cable...
cable_data = set_cable_sizing(cable_data, S=10.0)
print(cable_data.cross_section_area)
print(cable_data.current_capacity)
print(cable_data.joule_integral)