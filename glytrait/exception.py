class GlyTraitError(Exception):
    """Base class for exceptions in this module."""


class StructureParseError(GlyTraitError):
    """Raised when a structure cannot be parsed."""


class InputError(GlyTraitError):
    """The input file format error."""


class FormulaError(GlyTraitError):
    """Raised if a formula is invalid."""


class SiaLinkageError(GlyTraitError):
    """Raised if a sialic acid linkage is not specified but used."""
