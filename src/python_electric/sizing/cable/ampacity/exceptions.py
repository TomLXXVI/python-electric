__all__ = [
    "CurrentOverflowError",
    "CurrentUnderflowError",
    "CSANotFoundError",
    "AmpacityError"
]

class CurrentOverflowError(Exception):
    pass


class CurrentUnderflowError(Exception):
    pass


class CSANotFoundError(Exception):
    pass


class AmpacityError(Exception):
    pass
