from __future__ import annotations

from collections.abc import Generator, Iterable
from enum import Enum, auto
from typing import Literal

from attrs import frozen, field
from glypy.io.glycoct import loads as glycoct_loads, GlycoCTError
from glypy.structure.glycan import Glycan as GlypyGlycan
from glypy.structure.glycan_composition import (
    GlycanComposition,
    to_iupac_lite,
    MonosaccharideResidue,
)

from .utils import get_mono_comp


N_glycan_core = GlycanComposition.parse("{Man:3; Glc2NAc:2}")
Glc2NAc = MonosaccharideResidue.from_iupac_lite("Glc2NAc")
Man = MonosaccharideResidue.from_iupac_lite("Man")
Gal = MonosaccharideResidue.from_iupac_lite("Gal")
Neu5Ac = MonosaccharideResidue.from_iupac_lite("Neu5Ac")


class StructureParseError(Exception):
    """Raised when a structure cannot be parsed."""


class GlycanType(Enum):
    """The type of glycan."""

    COMPLEX = auto()
    HIGH_MANNOSE = auto()
    HYBRID = auto()


@frozen
class NGlycan:
    """A glycan."""

    _glypy_glycan: GlypyGlycan = field(repr=False)
    _composition: GlycanComposition = field(init=False)

    def __attrs_post_init__(self):
        object.__setattr__(
            self,
            "_composition",
            GlycanComposition.from_glycan(self._glypy_glycan),
        )

    @classmethod
    def from_string(cls, string: str, format: Literal["glycoct"] = "glycoct") -> NGlycan:
        """Build a glycan from a string representation.

        Args:
            string (str): The string representation of the glycan.
            format (Literal["glycoct"], optional): The format of the string.
                Defaults to "glycoct".

        Returns:
            NGlycan: The glycan.
        """
        if format == "glycoct":
            try:
                return cls(glycoct_loads(string))
            except GlycoCTError:
                raise StructureParseError(f"Could not parse string: {string}")
        else:
            raise ValueError(f"Unknown format: {format}")

    @classmethod
    def from_glycoct(cls, glycoct: str) -> NGlycan:
        """Build a glycan from a GlycoCT string.

        Args:
            glycoct (str): The GlycoCT string.

        Returns:
            NGlycan: The glycan.
        """
        return cls.from_string(glycoct, format="glycoct")

    @property
    def type(self) -> GlycanType:
        """The type of the glycan. Either 'complex', 'high-mannose' and 'hybrid'."""
        # N-glycan core is defined as the complex type.
        if self._composition == N_glycan_core:
            return GlycanType.COMPLEX

        # Bisecting could only be found in complex type.
        if self.is_bisecting:
            return GlycanType.COMPLEX

        if self._composition[Glc2NAc] == 2:
            return GlycanType.HIGH_MANNOSE

        bft_iter = self._breadth_first_traversal(skip=["Fuc"])
        for i in range(3):
            next(bft_iter)
        node1 = next(bft_iter)
        node2 = next(bft_iter)
        if any((len(node1.links) == 1, len(node2.links) == 1)):
            return GlycanType.COMPLEX

        if self._composition[Glc2NAc] == 3:
            return GlycanType.HYBRID

        return GlycanType.COMPLEX

    def is_complex(self) -> bool:
        """Whether the glycan is complex."""
        return self.type == GlycanType.COMPLEX

    def is_high_mannose(self) -> bool:
        """Whether the glycan is high-mannose."""
        return self.type == GlycanType.HIGH_MANNOSE

    def is_hybrid(self) -> bool:
        """Whether the glycan is hybrid."""
        return self.type == GlycanType.HYBRID

    @property
    def is_bisecting(self) -> bool:
        """Whether the glycan has a bisection."""
        bft_iter = self._breadth_first_traversal(skip=["Fuc"])
        for i in range(2):
            next(bft_iter)
        next_node = next(bft_iter)
        return len(next_node.links) == 4

    def _breadth_first_traversal(
        self, skip: Iterable[str] = None
    ) -> Generator[MonosaccharideResidue, None, None]:
        if skip is None:
            skip = []
        for node in self._glypy_glycan.breadth_first_traversal():
            if get_mono_comp(node) not in skip:
                yield node

    def __repr__(self) -> str:
        return f"NGlycan({to_iupac_lite(self._composition)})"
