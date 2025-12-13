from enum import StrEnum

__all__ = ["Soil", "SOIL_RESISTIVITY"]


class Soil(StrEnum):
    PEAT = "peat"
    LOAM = "loam"
    CLAY = "clay"
    HUMUS = "humus"
    MOIST_SAND = "moist_sand"
    MOIST_GRAVEL = "moist_gravel"
    DRY_SAND = "dry_sand"
    DRY_GRAVEL = "dry_gravel"
    STONEY_SOIL = "stoney_soil"


SOIL_RESISTIVITY = {  # ohm.m²/m
    "peat": 50,
    "loam": 100,
    "clay": 100,
    "humus": 100,
    "moist_sand": 200,
    "moist_gravel": 500,
    "dry_sand": 1000,
    "dry_gravel": 1000,
    "stoney_soil": 3000
}
