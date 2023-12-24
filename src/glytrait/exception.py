"""This module contains the set of GlyTrait's exceptions."""


class GlyTraitError(Exception):
    """Base class for exceptions in this module."""


class StructureParseError(GlyTraitError):
    """Raised when a structure cannot be parsed."""


class CompositionParseError(GlyTraitError):
    """Raised when a composition cannot be parsed."""


class InputError(GlyTraitError):
    """The input file format error."""


class FormulaError(GlyTraitError):
    """Raised if a formula is invalid."""


class SiaLinkageError(GlyTraitError):
    """Raised if a sialic acid linkage is not specified but used."""


class HypothesisTestingError(GlyTraitError):
    """Raised if a hypothesis test is not possible."""


class ConfigError(GlyTraitError):
    """Raised if a parameter is invalid."""
