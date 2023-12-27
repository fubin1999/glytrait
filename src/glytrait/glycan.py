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
from collections.abc import Generator, Mapping, Iterator, Callable
from typing import Literal, Iterable, Optional, Final, Union, TypeVar, NewType

from attrs import frozen, field
from glypy.io.glycoct import loads as glycoct_loads, GlycoCTError  # type: ignore
from glypy.structure.glycan import Glycan as GlypyGlycan  # type: ignore
from glypy.structure.glycan_composition import (  # type: ignore
    GlycanComposition,
    MonosaccharideResidue,
)

from glytrait.exception import (
    GlycanParseError,
    StructureParseError,
    CompositionParseError,
)

StructureDict = NewType("StructureDict", dict[str, "Structure"])
"""A type hint for a dictionary of `Structure` instances.
The keys are the names of the structures, and the values are the structures.
"""

CompositionDict = NewType("CompositionDict", dict[str, "Composition"])
"""A type hint for a dictionary of `Composition` instances.
The keys are the names of the compositions, and the values are the compositions.
"""


def parse_structures(__iter: Iterable[tuple[str, str]]) -> StructureDict:
    """Parse glycan structures from a list of structure strings.

    Args:
        An iterable of tuples of (name, structure_string).

    Returns:
        A dictionary of `Structure` instances, with the keys being the names of the structures.

    Raises:
        StructureParseError: When the string cannot be parsed.
            The names of the structures that cannot be parsed are included
            in the error message.
    """
    try:
        result: dict[str, Structure] = _load_glycans(__iter, Structure.from_string)
    except GlycanParseError as exc:
        raise StructureParseError(f"Could not parse structures for: {exc}.")
    else:
        return StructureDict(result)


def parse_compositions(__iter: Iterable[tuple[str, str]]) -> CompositionDict:
    """Parse glycan compositions from a list of composition strings.

    Args:
        An iterable of tuples of (name, composition_string).

    Returns:
        A dictionary of `Composition` instances, with the keys being the names of the compositions.

    Raises:
        CompositionParseError: When the string cannot be parsed.
            The names of the compositions that cannot be parsed are included
            in the error message.
    """
    try:
        result: dict[str, Composition] = _load_glycans(__iter, Composition.from_string)
    except GlycanParseError as exc:
        raise CompositionParseError(f"Could not parse compositions for: {exc}.")
    else:
        return CompositionDict(result)


Glycan = TypeVar("Glycan", bound=Union["Structure", "Composition"])
GlycanBuilder = Callable[[str, str], Glycan]


def _load_glycans(
    __iter: Iterable[tuple[str, str]], builder: GlycanBuilder
) -> dict[str, Glycan]:
    failed_names: list[str] = []
    glycans: dict[str, Glycan] = {}
    for name, string in __iter:
        try:
            glycan = builder(name, string)
        except GlycanParseError:
            failed_names.append(name)
        else:
            glycans[name] = glycan
    if failed_names:
        msg = ", ".join(f"'{name}'" for name in failed_names)
        raise StructureParseError(msg)
    return glycans


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
    _composition: dict[str, int] = field(init=False, repr=False, hash=False)

    def __attrs_post_init__(self):
        self._init_composition()

    def _init_composition(self):
        glypy_comp = GlycanComposition.from_glycan(self._glypy_glycan)
        comp = {str(k): v for k, v in glypy_comp.items()}
        object.__setattr__(self, "_composition", comp)

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
            Structure: The glycan.

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
                if get_mono_str(node) in only:  # type: ignore
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

    @property
    def composition(self) -> dict[str, int]:
        """The composition of the glycan."""
        return self._composition.copy()


VALID_MONOS: Final = ["H", "N", "F", "S", "L", "E"]


# The order here are used in the string representation of a composition.


@frozen
class Composition(Mapping[str, int]):
    """A glycan composition.

    Valid monosaccharides are: H, N, F, S, L, E.
    S must not be with L or E.
    Numbers of monosaccharides must be above 0.

    Attributes:
        name (str): The name of the glycan.

    Methods:
        from_string: classmethod to build a `Composition` instance from a string.

    Examples:
        >>> comp = Composition("G", {"H": 5, "N": 4, "F": 1})
        >>> comp
        Composition(name='G', comp={'H': 5, 'N': 4, 'F': 1})
        >>> str(comp)
        "H5N4F1"
        >>> comp["H"]
        5
        >>> comp["S"]
        0  # note that 'S' is not used when initialization
        >>> comp["P"]
        # raise KeyError
        >>> comp["H"] = 6
        # raise TypeError (unchangeable)
        >>> comp.get("P", default=0)
        0  # default value
        >>> len(comp)
        10  # 5 + 4 + 1
        >>> for mono, num in comp.items():
        ...     print(mono, num)
        H 5
        N 4
        F 1
        >>> Composition("G", {"H": 5, "N":4, "S": 1, "L": 1})
        # raise CompositionParseError (S must not be with L or E)
    """

    name: str = field()
    _comp: dict[str, int] = field(converter=dict, hash=False)

    def __attrs_post_init__(self):
        self._validate_comp()
        self._remove_zero()
        self._reorder()

    def _validate_comp(self) -> None:
        """Validate the composition."""
        self._validate_mono()
        self._validate_num()
        self._validate_sia()

    def _validate_mono(self) -> None:
        """Make sure that the monosaccharides are valid."""
        for k in self._comp:
            if k not in VALID_MONOS:
                raise CompositionParseError(f"Unknown monosaccharide: {k}.")

    def _validate_num(self) -> None:
        """Make sure that the numbers are valid."""
        for k, v in self._comp.items():
            if v < 0:
                raise CompositionParseError(f"Monosacharride must be above 0: {k}={v}.")

    def _validate_sia(self) -> None:
        """Make sure that S is not with L or E."""
        has_S = self._comp.get("S", 0) > 0
        has_L = self._comp.get("L", 0) > 0
        has_E = self._comp.get("E", 0) > 0
        if has_S and (has_L or has_E):
            raise CompositionParseError("S must not be with L or E.")

    def _remove_zero(self) -> None:
        to_delete: set[str] = set()
        for k, v in self._comp.items():
            if v == 0:
                to_delete.add(k)
        for k in to_delete:
            del self._comp[k]

    def _reorder(self) -> None:
        new_comp: dict[str, int] = {}
        for mono in VALID_MONOS:
            if mono in self._comp:
                new_comp[mono] = self._comp[mono]
        object.__setattr__(self, "_comp", new_comp)

    @classmethod
    def from_string(cls, name: str, string: str) -> Composition:
        """Create a composition from a string.

        Args:
            name (str): The name of the glycan.
            string (str): The string representation of the composition.

        Returns:
            Composition: The composition.

        Raises:
            CompositionParseError: When the string cannot be parsed.
        """
        if string == "":
            raise CompositionParseError("Empty string.")
        pattern = r"^([A-Z]\d+)*$"
        if not re.fullmatch(pattern, string):
            raise CompositionParseError(f"Invalid composition: {string}.")
        mono_comp: dict[str, int] = {}
        pattern = r"([A-Z])(\d+)"
        for m in re.finditer(pattern, string):
            mono_comp[m.group(1)] = int(m.group(2))
        return cls(name, mono_comp)  # type: ignore

    def __getitem__(self, __key: str) -> int:
        try:
            return self._comp[__key]
        except KeyError:
            if __key in VALID_MONOS:
                return 0
            else:
                raise

    def __len__(self) -> int:
        return sum(self._comp.values())

    def __iter__(self) -> Iterator[str]:
        return iter(self._comp)

    def __str__(self) -> str:
        mono_strings: list[str] = []
        for mono, num in self._comp.items():
            mono_strings.append(f"{mono}{num}")
        return "".join(mono_strings)


def get_mono_str(mono: MonosaccharideResidue) -> str:
    """Get the string representation of a monosaccharide residue.

    This is a helper function for checking the identity of a `glypy.MonosaccarideResidue`
    instance.

    Args:
        mono (glypy.MonosaccharideResidue): The monosaccharide residue.
    """
    return MonosaccharideResidue.from_monosaccharide(mono).name()
