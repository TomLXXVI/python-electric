"""
Implementation of the ampacity table for unburied cables in amperes.

References
----------
Electrical Installation Guide 2007, Merlin Gerin, Chapter G, Fig. G21a.
"""
from ....materials import ConductorProperties, InsulationProperties
from ....utils.lookup_table import LookupTable


def _create_index_table():

    def _create_row(
        PVC2: int,
        PVC3: int,
        XLPE2: int,
        XLPE3: int
    ) -> dict[str, int]:
        return {
            "PVC2": PVC2,
            "PVC3": PVC3,
            "XLPE2": XLPE2,
            "XLPE3": XLPE3
        }

    d = {
        "A1": _create_row(1, 2, 5, 4),
        "A2": _create_row(1, 0, 4, 3),
        "B1": _create_row(4, 3, 8, 6),
        "B2": _create_row(3, 2, 6, 5),
        "C": _create_row(6, 4, 9, 7),
        "E": _create_row(7, 5, 10, 8),
        "F": _create_row(8, 6, 11, 9)
    }
    return d


def _create_copper_table() -> LookupTable:
    col_header = [j for j in range(12)]
    row_header = [1.5, 2.5, 4, 6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240]
    nan = float("nan")
    data = [
        [13.0, 13.5, 14.5,  15.5,  17.0,  18.5,  20.0,  22.0,  23.0,  24.0,  26.0,   nan],
        [17.5, 18.0, 19.5,  21.0,  23.0,  25.0,  27.0,  30.0,  31.0,  33.0,  36.0,   nan],
        [23.0, 24.0, 26.0,  28.0,  31.0,  34.0,  36.0,  40.0,  42.0,  45.0,  49.0,   nan],
        [29.0, 31.0, 34.0,  36.0,  40.0,  43.0,  46.0,  51.0,  54.0,  58.0,  63.0,   nan],
        [39.0, 42.0, 46.0,  50.0,  54.0,  60.0,  63.0,  70.0,  75.0,  80.0,  86.0,   nan],
        [52.0, 56.0, 61.0,  68.0,  73.0,  80.0,  85.0,  94.0, 100.0, 107.0, 115.0,   nan],
        [68.0, 73.0, 80.0,  89.0,  95.0, 101.0, 110.0, 119.0, 127.0, 135.0, 149.0, 161.0],
        [ nan,  nan,  nan, 110.0, 117.0, 126.0, 137.0, 147.0, 158.0, 169.0, 185.0, 200.0],
        [ nan,  nan,  nan, 134.0, 141.0, 153.0, 167.0, 179.0, 192.0, 207.0, 225.0, 242.0],
        [ nan,  nan,  nan, 171.0, 179.0, 196.0, 213.0, 229.0, 246.0, 268.0, 289.0, 310.0],
        [ nan,  nan,  nan, 207.0, 216.0, 238.0, 258.0, 278.0, 298.0, 328.0, 352.0, 377.0],
        [ nan,  nan,  nan, 239.0, 249.0, 276.0, 299.0, 322.0, 346.0, 382.0, 410.0, 437.0],
        [ nan,  nan,  nan,   nan, 285.0, 318.0, 344.0, 371.0, 395.0, 441.0, 473.0, 504.0],
        [ nan,  nan,  nan,   nan, 324.0, 362.0, 392.0, 424.0, 450.0, 506.0, 542.0, 575.0],
        [ nan,  nan,  nan,   nan, 380.0, 424.0, 461.0, 500.0, 538.0, 599.0, 641.0, 679.0]
    ]
    description = "current-carrying capacities in amperes, unburied, copper"
    cols_description = "number of loaded conductors and type of insulation"
    rows_description = "cross-sectional area"
    table = LookupTable(
        data, col_header, row_header,
        description,
        cols_description,
        rows_description
    )
    return table


def _create_aluminium_table() -> LookupTable:
    col_header = [j for j in range(12)]
    row_header = [2.5, 4, 6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240]
    nan = float("nan")
    data = [
        [13.5, 14.0, 15.0,  16.5,  18.5,  19.5,  21.0,  23.0,  24.0,  26.0,  28.0,   nan],
        [17.5, 18.5, 20.0,  22.0,  25.0,  26.0,  28.0,  31.0,  32.0,  35.0,  38.0,   nan],
        [23.0, 24.0, 26.0,  28.0,  32.0,  33.0,  36.0,  39.0,  42.0,  45.0,  49.0,   nan],
        [31.0, 32.0, 36.0,  39.0,  44.0,  46.0,  49.0,  54.0,  58.0,  62.0,  67.0,   nan],
        [41.0, 43.0, 48.0,  53.0,  58.0,  61.0,  66.0,  73.0,  77.0,  84.0,  91.0,   nan],
        [53.0, 57.0, 63.0,  70.0,  73.0,  78.0,  83.0,  90.0,  97.0, 101.0, 108.0, 121.0],
        [ nan,  nan,  nan,  86.0,  90.0,  96.0, 103.0, 112.0, 120.0, 126.0, 135.0, 150.0],
        [ nan,  nan,  nan, 104.0, 110.0, 117.0, 125.0, 136.0, 146.0, 154.0, 164.0, 184.0],
        [ nan,  nan,  nan, 133.0, 140.0, 150.0, 160.0, 174.0, 187.0, 198.0, 211.0, 237.0],
        [ nan,  nan,  nan, 161.0, 170.0, 183.0, 195.0, 211.0, 227.0, 241.0, 257.0, 289.0],
        [ nan,  nan,  nan, 186.0, 197.0, 212.0, 226.0, 245.0, 263.0, 280.0, 300.0, 337.0],
        [ nan,  nan,  nan,   nan, 226.0, 245.0, 261.0, 283.0, 304.0, 324.0, 346.0, 389.0],
        [ nan,  nan,  nan,   nan, 256.0, 280.0, 298.0, 323.0, 347.0, 371.0, 397.0, 447.0],
        [ nan,  nan,  nan,   nan, 300.0, 330.0, 352.0, 382.0, 409.0, 439.0, 470.0, 530.0]
    ]
    description = "current-carrying capacities in amperes, unburied, aluminium"
    cols_description = "number of loaded conductors and type of insulation"
    rows_description = "cross-sectional area"
    table = LookupTable(
        data, col_header, row_header,
        description,
        cols_description,
        rows_description
    )
    return table


index_table = _create_index_table()
copper_table = _create_copper_table()
alu_table = _create_aluminium_table()

def get_cross_sectional_area(
    conductor_props: ConductorProperties,
    insulation_props: InsulationProperties,
    num_loaded_conductors: int,
    ref_method: str,
    current: float
) -> float:
    if conductor_props.type in ["copper", "aluminium"]:
        if insulation_props.type in ["PVC", "XLPE"]:
            if 1 < num_loaded_conductors <= 3:
                key = insulation_props.type + str(num_loaded_conductors)
                col_index = index_table[ref_method][key]
                if conductor_props.type == "copper":
                    csa = copper_table.rowheader_value(col_index, current)
                else:
                    csa = alu_table.rowheader_value(col_index, current)
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
    ref_method: str,
    cross_sectional_area: float
) -> float:
    if conductor_props.type == "copper" or conductor_props.type == "aluminium":
        if insulation_props.type in ["PVC", "XLPE"]:
            if 1 < num_loaded_conductors <= 3:
                key = insulation_props.type + str(num_loaded_conductors)
                col_index = index_table[ref_method][key]
                if conductor_props.type == "copper":
                    amp = copper_table.data_value(cross_sectional_area, col_index)
                else:
                    amp = alu_table.data_value(cross_sectional_area, col_index)
                return amp
            else:
                raise ValueError(
                    "The number of loaded conductors is limited to 2 or 3."
                )
        else:
            raise ValueError(f"{insulation_props.type} is not supported.")
    else:
        raise ValueError("Only copper and aluminium conductors are supported.")
