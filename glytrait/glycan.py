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
from glypy.structure.monosaccharide import Monosaccharide

from .utils import get_mono_comp

N_glycan_core = GlycanComposition.parse("{Man:3; Glc2NAc:2}")
Glc2NAc = MonosaccharideResidue.from_iupac_lite("Glc2NAc")
Man = MonosaccharideResidue.from_iupac_lite("Man")
Gal = MonosaccharideResidue.from_iupac_lite("Gal")
Neu5Ac = MonosaccharideResidue.from_iupac_lite("Neu5Ac")
Neu5Gc = MonosaccharideResidue.from_iupac_lite("Neu5Gc")
Fuc = MonosaccharideResidue.from_iupac_lite("Fuc")


class StructureParseError(Exception):
    """Raised when a structure cannot be parsed."""


class BranchError(Exception):
    """Raised when `count_branch` is called on non-complex glycans."""


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
    _cores: list[int] = field(init=False)

    def __attrs_post_init__(self):
        self._init_composition()
        self._init_cores()

    def _init_composition(self):
        object.__setattr__(
            self,
            "_composition",
            GlycanComposition.from_glycan(self._glypy_glycan),
        )

    def _init_cores(self):
        should_be = ["Glc2NAc", "Glc2NAc", "Man", "Man", "Man"]
        cores: list[int] = []
        for node in self._breadth_first_traversal(skip=["Fuc"]):
            # The first two monosaacharides could only be "Glc2NAc".
            # However, when the glycan is bisecting, the rest of the three monosaccharides
            # might not all be "Man". So we look for the monosaacharides in the order of
            # `should_be`, and skip the ones that are not in the order.
            if get_mono_comp(node) == should_be[len(cores)]:
                cores.append(node.id)
            if len(cores) == 5:
                break
        self._check_cores(cores)
        object.__setattr__(self, "_cores", cores)

    def _check_cores(self, cores):
        """Check the validation of the cores."""
        core_nodes = [self._glypy_glycan.get(i) for i in cores]
        core_residues = [
            MonosaccharideResidue.from_monosaccharide(n) for n in core_nodes
        ]
        core_comps = [get_mono_comp(n) for n in core_residues]
        N_glycan_core_comps = ["Glc2NAc", "Glc2NAc", "Man", "Man", "Man"]
        if sorted(core_comps) != N_glycan_core_comps:
            raise ValueError(f"Invalid core: {core_comps}")

    @classmethod
    def from_string(
        cls, string: str, format: Literal["glycoct"] = "glycoct"
    ) -> NGlycan:
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

    def _breadth_first_traversal(
        self, skip: Iterable[str] = None
    ) -> Generator[MonosaccharideResidue, None, None]:
        if skip is None:
            skip = []
        for node in self._glypy_glycan.breadth_first_traversal():
            if get_mono_comp(node) not in skip:
                yield node

    def _get_branch_core_man(self) -> tuple[Monosaccharide, Monosaccharide]:
        # Branch core mannose is defined as:
        #
        # X - Man  <- This one
        #        \
        #         Man - Glc2NAc - Glc2NAc
        #        /
        # X - Man  <- And this one
        bft_iter = self._breadth_first_traversal(skip=["Fuc"])
        for i in range(3):
            next(bft_iter)
        while True:
            node1 = next(bft_iter)
            if get_mono_comp(node1) == "Man":
                break
        while True:
            node2 = next(bft_iter)
            if get_mono_comp(node2) == "Man":
                break
        return node1, node2

    @property
    def type(self) -> GlycanType:
        """The type of the glycan. Either 'complex', 'high-mannose' and 'hybrid'."""
        # N-glycan core is defined as the complex type.
        if self._composition == N_glycan_core:
            return GlycanType.COMPLEX

        # Bisecting could only be found in complex type.
        if self.is_bisecting():
            return GlycanType.COMPLEX

        # If the glycan is not core, and it only has 2 "GlcNAc", it is high-mannose.
        if self._composition[Glc2NAc] == 2:
            return GlycanType.HIGH_MANNOSE

        # If the glycan is mono-antennary and not high-monnose, it is complex.
        node1, node2 = self._get_branch_core_man()
        if any((len(node1.links) == 1, len(node2.links) == 1)):
            return GlycanType.COMPLEX

        # Then, if it has 3 "Glc2NAc", it must be hybrid.
        if self._composition[Glc2NAc] == 3:
            return GlycanType.HYBRID

        # All rest cases are complex.
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

    def is_bisecting(self) -> bool:
        """Whether the glycan has a bisection."""
        bft_iter = self._breadth_first_traversal(skip=["Fuc"])
        for i in range(2):
            next(bft_iter)
        next_node = next(bft_iter)
        return len(next_node.links) == 4

    def count_antenna(self) -> int:
        """The number of branches in the glycan."""
        if not self.is_complex():
            raise BranchError("Cannot count branches on non-complex glycans.")
        node1, node2 = self._get_branch_core_man()
        return len(node1.links) + len(node2.links) - 2

    def count_fuc(self) -> int:
        """The number of fucoses."""
        return self._composition[Fuc]

    def count_core_fuc(self) -> int:
        """The number of core fucoses."""
        n = 0
        for node in self._breadth_first_traversal():
            # node.parents()[0] is the nearest parent, and is a tuple of (Link, Monosaccharide)
            # node.parents()[0][1] is the Monosaccharide
            if get_mono_comp(node) == "Fuc" and node.parents()[0][1].id in self._cores:
                n = n + 1
        return n

    def count_antennary_fuc(self) -> int:
        """The number of antennary fucoses."""
        return self.count_fuc() - self.count_core_fuc()

    def count_sia(self) -> int:
        """The number of sialic acids."""
        return self._composition[Neu5Ac] + self._composition[Neu5Gc]

    def count_man(self) -> int:
        """The number of mannoses."""
        return self._composition[Man]

    def count_gal(self) -> int:
        """The number of galactoses."""
        return self._composition[Gal]

    def __repr__(self) -> str:
        return f"NGlycan({to_iupac_lite(self._composition)})"
