"""This module implements classes about representing glycans.

Classes:
    NGlycan: A class representing an N-glycan structure.
    Composition: A class representing the composition of a glycan.

Functions:
    load_glycans: Load glycan structures from a list of structure strings.
    load_compositions: Load glycan compositions from a list of composition strings.
"""
from __future__ import annotations

import re
from collections.abc import Generator
from enum import Enum, auto
from typing import Literal, NoReturn, Iterable, Optional

from attrs import frozen, field
from glypy.io.glycoct import loads as glycoct_loads, GlycoCTError
from glypy.structure.glycan import Glycan as GlypyGlycan
from glypy.structure.glycan_composition import (
    GlycanComposition,
    MonosaccharideResidue,
)

from glytrait.exception import *


class GlycanType(Enum):
    """The type of glycan."""

    COMPLEX = auto()
    HIGH_MANNOSE = auto()
    HYBRID = auto()


@frozen
class Structure:
    """The structure of a glycan.

    Attributes:
        name (str): The name of the glycan.
        composition (dict): The composition of the glycan.

    Methods:
        from_string: classmethod for building a `Structure` instance from a string.
        from_glycoct: classmethod for building a `Structure` instance from a glycoCT string.
        breadth_first_traversal: traverse the structure in a "bfs" manner.
        depth_first_traversal: traverse the structure in a "dfs" manner.
    """

    name: str = field()
    _glypy_glycan: GlypyGlycan = field(repr=False)
    composition: dict[str, int] = field(init=False, repr=False, hash=False)

    def __attrs_post_init__(self):
        self._init_composition()

    def _init_composition(self):
        glypy_comp = GlycanComposition.from_glycan(self._glypy_glycan)
        comp = {str(k): v for k, v in glypy_comp.items()}
        object.__setattr__(self, "composition", comp)

    @classmethod
    def from_string(
        cls, name: str, string: str, *, format: Literal["glycoct"] = "glycoct"
    ) -> Structure:
        """Build a glycan from a string representation.

        Args:
            name (str): The name of the glycan.
            string (str): The string representation of the glycan.
            format (Literal["glycoct"], optional): The format of the string.
                Defaults to "glycoct".

        Returns:
            NGlycan: The glycan.

        Raises:
            StructureParseError: When the string cannot be parsed.
        """
        if format == "glycoct":
            try:
                return cls(name, glycoct_loads(string))
            except GlycoCTError:
                raise StructureParseError(f"Could not parse string: {string}")
        else:
            raise StructureParseError(f"Unknown format: {format}")

    @classmethod
    def from_glycoct(cls, name: str, glycoct: str) -> Structure:
        """Build a glycan from a GlycoCT string.

        Args:
            name (str): The name of the glycan.
            glycoct (str): The GlycoCT string.

        Returns:
            NGlycan: The glycan.

        Raises:
            StructureParseError: When the string cannot be parsed.
        """
        return cls.from_string(name, glycoct, format="glycoct")

    def _traversal(
        self,
        method: Literal["bfs", "dfs"],
        *,
        skip: Optional[Iterable[str]] = None,
        only: Optional[Iterable[str]] = None,
    ) -> Generator[MonosaccharideResidue, None, None]:
        # set the traversal method
        if method == "bfs":
            traversal_func = self._glypy_glycan.breadth_first_traversal
        elif method == "dfs":
            traversal_func = self._glypy_glycan.depth_first_traversal
        else:
            raise ValueError(f"Unknown traversal method: {method}")

        # check the validation of `skip` and `only`
        if skip and only:
            raise ValueError("Cannot specify both `skip` and `only`.")

        # traverse the glycan
        if skip is None and only is None:
            yield from traversal_func()
        elif skip is not None:
            for node in traversal_func():
                if get_mono_str(node) not in skip:
                    yield node
        else:  # only is not None
            for node in traversal_func():
                if get_mono_str(node) in only:
                    yield node

    def breadth_first_traversal(
        self,
        *,
        skip: Optional[Iterable[str]] = None,
        only: Optional[Iterable[str]] = None,
    ) -> Generator[MonosaccharideResidue, None, None]:
        """Traverse the structure in a "bfs" manner.

        Args:
            skip (Optional[Iterable[str]], optional): The monosaccharides to skip.
                Defaults to None.
            only (Optional[Iterable[str]], optional): The monosaccharides to traverse.
                Defaults to None.

        Yields:
            glypy.MonosaccharideResidue: The monosaccharide residues.
        """
        yield from self._traversal("bfs", skip=skip, only=only)

    def depth_first_traversal(
        self,
        *,
        skip: Optional[Iterable[str]] = None,
        only: Optional[Iterable[str]] = None,
    ) -> Generator[MonosaccharideResidue, None, None]:
        """Traverse the structure in a "dfs" manner.

        Args:
            skip (Optional[Iterable[str]], optional): The monosaccharides to skip.
                Defaults to None.
            only (Optional[Iterable[str]], optional): The monosaccharides to traverse.
                Defaults to None.

        Yields:
            glypy.MonosaccharideResidue: The monosaccharide residues.
        """
        yield from self._traversal("dfs", skip=skip, only=only)

    def get(self, key: str, default=0):
        """Get the number of the given monosaccharide."""
        return self.composition.get(key, default)


@frozen
class Composition:
    """A glycan composition.

    Attributes:
        name (str): The name of the glycan.

    Methods:
        from_string: classmethod to build a `Composition` instance from a string.
        asdict: convert the `Composition` instance into a dict.
        get: get the number of the given monosaccharide.
    """

    name: str = field()
    _comp: dict[str, int] = field(converter=dict, hash=False)

    @_comp.validator
    def _check_comp(self, attribute, value):
        """Check the composition."""
        for k, v in value.items():
            if k not in {"H", "N", "F", "S", "L", "E"}:
                raise CompositionParseError(f"Unknown monosaccharide: {k}.")
            if v < 0:
                raise CompositionParseError(f"Monosacharride must be above 0: {k}={v}.")

    @classmethod
    def from_string(cls, s: str) -> Composition:
        """Create a composition from a string.

        Args:
            s (str): The string, e.g. H5N4F1S1.

        Returns:
            Composition: The composition.

        Raises:
            CompositionParseError: When the string cannot be parsed.
        """
        cls._validate_string(s)
        # noinspection PyUnreachableCode
        mono_comp: dict[str, int] = {}
        pattern = r"([A-Z])(\d+)"
        for m in re.finditer(pattern, s):
            mono_comp[m.group(1)] = int(m.group(2))
        return cls(s, mono_comp)

    @staticmethod
    def _validate_string(s: str) -> NoReturn:
        if s == "":
            raise CompositionParseError("Empty string.")
        pattern = r"^([A-Z]\d+)*$"
        if not re.fullmatch(pattern, s):
            raise CompositionParseError(f"Invalid composition: {s}.")

    def asdict(self) -> dict[str, int]:
        """Return the composition as a dict."""
        return self._comp.copy()

    def get(self, key: str, default: int = 0):
        """Get the number of the given monosaccharide."""
        return self._comp.get(key, default)


def get_mono_str(mono: MonosaccharideResidue) -> str:
    """Get the string representation of a monosaccharide residue.

    This is a helper function for checking the identity of a `glypy.MonosaccarideResidue`
    instance.

    Args:
        mono (glypy.MonosaccharideResidue): The monosaccharide residue.
    """
    return MonosaccharideResidue.from_monosaccharide(mono).name()


def load_structures(names: Iterable[str], structures: Iterable[str]) -> list[Structure]:
    """Load structures from a list of structures.

    Args:
        names (Iterable[str]): The names of the glycans.
        structures (Iterable[str]): The structure strings of the glycans.

    Returns:
        list[Structure]: The glycan structures.
    """
    return [Structure.from_string(n, s) for n, s in zip(names, structures)]


def load_compositions(compositions: Iterable[str]) -> list[Composition]:
    """Load compositions from a list of strings.

    Args:
        compositions (Iterable[str]): The compositions.

    Returns:
        list[Composition]: The compositions.
    """
    return [Composition.from_string(s) for s in compositions]
