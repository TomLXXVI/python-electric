__all__ = [
    "DataLowerBoundError",
    "DataUpperBoundError",
    "DataNotFoundError",
    "RowheaderNotFoundError",
    "ColumnheaderNotFoundError",
]


class DataLowerBoundError(Exception):
    pass


class DataUpperBoundError(Exception):
    pass


class DataNotFoundError(Exception):
    pass


class RowheaderNotFoundError(Exception):
    pass


class ColumnheaderNotFoundError(Exception):
    pass
