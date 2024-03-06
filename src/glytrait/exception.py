"""This module contains the set of GlyTrait's exceptions."""


class GlyTraitError(Exception):
    """Base class for exceptions in this module."""


class GlycanParseError(GlyTraitError):
    """Raised when a glycan cannot be parsed."""


class StructureParseError(GlycanParseError):
    """Raised when a structure cannot be parsed."""


class CompositionParseError(GlycanParseError):
    """Raised when a composition cannot be parsed."""


class DataInputError(GlyTraitError):
    """The input file format error."""


class FileTypeError(GlyTraitError):
    """Raised when a file type is not supported."""


class NotEnoughGroupsError(GlyTraitError):
    """Raised when there are not enough groups."""


class FormulaError(GlyTraitError):
    """Raised if a formula is invalid."""


class SiaLinkageError(GlyTraitError):
    """Raised if a sialic acid linkage is not specified but used."""
