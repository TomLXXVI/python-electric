"""
Implementation of the ampacity table for buried cables in amperes.

References
----------
Electrical Installation Guide 2007, Merlin Gerin, Chapter G, Fig. G22.
"""
from ....materials import ConductorProperties, InsulationProperties
from ....utils.lookup_table import LookupTable


def _create_copper_table():
    col_header = ["PVC3", "PVC2", "XLPE3", "XLPE2"]
    row_header = [1.5, 2.5, 4, 6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240, 300]
    data = [
        [18.0,   22.0,  22.0,  26.0],
        [24.0,   29.0,  29.0,  34.0],
        [31.0,   38.0,  37.0,  44.0],
        [39.0,   47.0,  46.0,  56.0],
        [52.0,   63.0,  61.0,  73.0],
        [67.0,   81.0,  79.0,  95.0],
        [86.0,  104.0, 101.0, 121.0],
        [103.0, 125.0, 122.0, 146.0],
        [122.0, 148.0, 144.0, 173.0],
        [151.0, 183.0, 178.0, 213.0],
        [179.0, 216.0, 211.0, 252.0],
        [203.0, 246.0, 240.0, 287.0],
        [230.0, 278.0, 271.0, 324.0],
        [258.0, 312.0, 304.0, 363.0],
        [297.0, 361.0, 351.0, 419.0],
        [336.0, 408.0, 396.0, 474.0]
    ]
    description = "current-carrying capacities in amperes, buried, copper"
    cols_description = "number of loaded conductors and type of insulation"
    rows_description = "cross-sectional area"
    table = LookupTable(
        data, col_header, row_header,
        description,
        cols_description,
        rows_description
    )
    return table


def _create_aluminium_table():
    col_header = ["PVC3", "PVC2", "XLPE3", "XLPE2"]
    row_header = [2.5, 4, 6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240, 300]
    data = [
        [18.5,   22.0,  22.0,  26.0],
        [24.0,   29.0,  29.0,  34.0],
        [30.0,   36.0,  36.0,  42.0],
        [40.0,   48.0,  47.0,  56.0],
        [52.0,   62.0,  61.0,  73.0],
        [66.0,   80.0,  78.0,  93.0],
        [80.0,   96.0,  94.0, 112.0],
        [94.0,  113.0, 112.0, 132.0],
        [117.0, 140.0, 138.0, 163.0],
        [138.0, 166.0, 164.0, 193.0],
        [157.0, 189.0, 186.0, 220.0],
        [178.0, 213.0, 210.0, 249.0],
        [200.0, 240.0, 236.0, 279.0],
        [230.0, 277.0, 272.0, 322.0],
        [260.0, 313.0, 308.0, 364.0]
    ]
    description = "current-carrying capacities in amperes, buried, aluminium"
    cols_description = "number of loaded conductors and type of insulation"
    rows_description = "cross-sectional area"
    table = LookupTable(
        data, col_header, row_header,
        description,
        cols_description,
        rows_description
    )
    return table


copper_table = _create_copper_table()
alu_table = _create_aluminium_table()

def get_cross_sectional_area(
    conductor_props: ConductorProperties,
    insulation_props: InsulationProperties,
    num_loaded_conductors: int,
    current: float
) -> float:
    if conductor_props.type in ["copper", "aluminium"]:
        if insulation_props.type in ["PVC", "XLPE"]:
            if 1 < num_loaded_conductors <= 3:
                key = insulation_props.type + str(num_loaded_conductors)
                if conductor_props.type == "copper":
                    csa = copper_table.rowheader_value(key, current)
                else:
                    csa = alu_table.rowheader_value(key, current)
                return csa
            else:
                raise ValueError(
                    "The number of loaded conductors is limited to 2 or 3."
                )
        else:
            raise ValueError(f"{insulation_props.type} is not supported.")
    else:
        raise ValueError("Only copper and aluminium conductors are supported.")


def get_ampacity(
    conductor_props: ConductorProperties,
    insulation_props: InsulationProperties,
    num_loaded_conductors: int,
    cross_sectional_area: float
) -> float:
    if conductor_props.type in ["copper", "aluminium"]:
        if insulation_props.type in ["PVC", "XLPE"]:
            if 1 < num_loaded_conductors <= 3:
                key = insulation_props.type + str(num_loaded_conductors)
                if conductor_props.type == "copper":
                    amp = copper_table.data_value(cross_sectional_area, key)
                else:
                    amp = alu_table.data_value(cross_sectional_area, key)
                return amp
            else:
                raise ValueError(
                    "The number of loaded conductors is limited to 2 or 3."
                )
        else:
            raise ValueError(f"{insulation_props.type} is not supported.")
    else:
        raise ValueError("Only copper and aluminium conductors are supported.")
