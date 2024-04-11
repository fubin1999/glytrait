"""This module provides functions for calculating meta-properties of glycans.

Functions:
    build_meta_property_table: Build a table of meta-properties for glycans.

Currently, the following meta-properties are supported:
    - type: The type of glycan. (structure only)
    - B: Whether the glycan has bisection. (structure only)
    - nAnt: The number of antennas. (structure only)
    - nF: The number of fucoses.
    - nFc: The number of fucoses on the core. (structure only)
    - nFa: The number of fucoses on the antenna. (structure only)
    - nS: The number of sialic acids.
    - nM: The number of mannoses.
    - nG: The number of galactoses.
    - nN: The number of GlcNAcs.
    - Pl: Whether the glycan has any poly-LacNAc. (structure only)
    - nL: The number of sialic acids with an alpha-2,3 linkage.
    - nE: The number of sialic acids with an alpha-2,6 linkage.
"""

from __future__ import annotations

from enum import Enum, auto
from functools import singledispatch, cache
from typing import Literal, Any, ClassVar, Type

import pandas as pd
from attrs import define
from glypy import Monosaccharide  # type: ignore

from glytrait.data_type import MetaPropertyTable
from glytrait.exception import SiaLinkageError
from glytrait.glycan import Structure, Composition, get_mono_str, GlycanDict

__all__ = [
    "build_meta_property_table",
    "available_meta_properties",
    "get_meta_property",
]


def build_meta_property_table(
    glycans: GlycanDict,
    mode: Literal["composition", "structure"],
    sia_linkage: bool = False,
) -> MetaPropertyTable:
    """Build a table of meta-properties for glycans.

    Args:
        glycans: (GlycanDict): A dict of glycans, with glycan IDs as keys and
            Compositions or Structures as values.
        mode (Literal["composition", "structure"]): The calculation mode.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta-properties. Defaults to False.

    Returns:
        MetaPropertyTable: The table of meta-properties.
    """
    mp_series_list: list[pd.Series] = []
    for mp_name in available_meta_properties(mode, sia_linkage):
        mp = get_meta_property(mp_name)
        s = mp.calculate_many(glycans)
        mp_series_list.append(s)
    mp_table_df = pd.concat(mp_series_list, axis=1)
    return MetaPropertyTable(mp_table_df)


class GlycanType(Enum):
    """The type of glycan."""

    COMPLEX = auto()
    HIGH_MANNOSE = auto()
    HYBRID = auto()


# ===== Registering meta-properties =====
Glycan = Composition | Structure

_mp_objects: dict[str, MetaProperty] = {}


def register(mp_type: Type[MetaProperty]) -> Type[MetaProperty]:
    """Decorator for registering meta-properties.

    Args:
        mp_type: The meta-property class.

    Returns:
        Type[MetaProperty]: The meta-property class.
    """
    if mp_type.name in _mp_objects:
        raise ValueError(f"Meta-property '{mp_type.name}' already exists")
    _mp_objects[mp_type.name] = mp_type()
    return mp_type


# ===== Get available meta-properties =====
def available_meta_properties(
    mode: Literal["composition", "structure"], sia_linkage: bool
) -> list[str]:
    """Get the available meta-property names.

    Args:
        mode: The calculation mode.
        sia_linkage: Whether to include the sialic acid linkage meta-properties.
            If True, all meta-properties will be included.
            Otherwise, "nE" and "nL" will be excluded.

    Returns:
        list[str]: The available meta-property names.
    """
    to_return: list[str] = []
    for name, mp in _mp_objects.items():
        mode_correct = mp.supported_mode in ("both", mode)  # type: ignore
        if sia_linkage:
            sia_linkage_correct = True
        else:
            sia_linkage_correct = name not in ("nL", "nE")
        if mode_correct and sia_linkage_correct:
            to_return.append(name)
    return to_return


def get_meta_property(name: str) -> MetaProperty:
    """Get a meta-property by its name.

    Args:
        name: The name of the meta-property.

    Returns:
        MetaProperty: The meta-property object.

    Raises:
        KeyError: If the meta-property does not exist.
    """
    return _mp_objects[name]


# ===== Base class for meta-properties =====
@define
class MetaProperty:
    """The base class for meta-properties.

    A meta-property is a basic property of a glycan, such as the number of fucoses,
    the number of sialic acids, etc.
    It is calculated from the glycan structure or composition.

    Besides three class attributes,
    subclasses should also implement the `calculate_one` method for calculating the
    meta-property for a single glycan.
    `calculate_many` is ready to use and should not be overridden.

    For consistency, all subclasses should be named with the suffix "MP".
    e.g. `GlycanTypeMP`, `BisectionMP`, etc.

    Attributes:
        name: The name of the meta-property.
        supported_mode: The supported calculation mode, which can be "composition",
            "structure", or "both".
        return_type: The type of the series returned by `calculate_many`, which can be
            "UInt8", "boolean", or "category".
    """

    name: ClassVar[str]
    supported_mode: ClassVar[Literal["composition", "structure", "both"]]
    return_type: ClassVar[Literal["UInt8", "boolean", "category"]]
    # "UInt8" is used mainly for the number of monosaccharides,
    # which is always non-negative and less than 256.

    def calculate_one(self, glycan: Glycan) -> int | bool | str:
        """Calculate the meta-property for a single glycan.

        Subclasses should always override this method.

        Args:
            glycan: The glycan.

        Returns:
            Any: The value of the meta-property.
        """
        raise NotImplementedError

    def calculate_many(self, glycans: GlycanDict) -> pd.Series:
        """Calculate the meta-property for multiple glycans.

        Args:
            glycans: A dict of glycans, with glycan IDs as keys and
                Compositions or Structures as values.

        Returns:
            pd.Series: A series of the meta-property values, with glycan IDs as the index.
        """
        values = [self.calculate_one(glycan) for glycan in glycans.values()]
        return pd.Series(
            data=values,
            index=pd.Index(glycans.keys()),
            name=self.name,
            dtype=self.return_type,
        )


# ===== Concrete meta-properties =====
@register
@define
class GlycanTypeMP(MetaProperty):
    """The type of glycan."""

    name: ClassVar = "type"
    supported_mode: ClassVar = "structure"
    return_type: ClassVar = "category"

    def calculate_one(self, glycan: Glycan) -> str:
        """Calculate the type of glycan."""
        return _glycan_type(glycan).name.lower()


@register
@define
class BisectionMP(MetaProperty):
    """Whether the glycan has bisection."""

    name: ClassVar = "B"
    supported_mode: ClassVar = "structure"
    return_type: ClassVar = "boolean"

    def calculate_one(self, glycan: Glycan) -> bool:
        """Calculate whether the glycan has bisection."""
        return _is_bisecting(glycan)


@register
@define
class CountAntennaMP(MetaProperty):
    """The number of antennas."""

    name: ClassVar = "nAnt"
    supported_mode: ClassVar = "structure"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of antennas."""
        return _count_antenna(glycan)


@register
@define
class CountFucMP(MetaProperty):
    """The number of fucoses."""

    name: ClassVar = "nF"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of fucoses."""
        return _count_fuc(glycan)


@register
@define
class CountCoreFucMP(MetaProperty):
    """The number of fucoses on the core."""

    name: ClassVar = "nFc"
    supported_mode: ClassVar = "structure"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of fucoses on the core."""
        return _count_core_fuc(glycan)


@register
@define
class CountAntennaryFucMP(MetaProperty):
    """The number of fucoses on the antenna."""

    name: ClassVar = "nFa"
    supported_mode: ClassVar = "structure"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of fucoses on the antenna."""
        return _count_fuc(glycan) - _count_core_fuc(glycan)


@register
@define
class CountSiaMP(MetaProperty):
    """The number of sialic acids."""

    name: ClassVar = "nS"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of sialic acids."""
        return _count_sia(glycan)


@register
@define
class CountManMP(MetaProperty):
    """The number of mannoses."""

    name: ClassVar = "nM"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of mannoses."""
        return _count_man(glycan)


@register
@define
class CountGalMP(MetaProperty):
    """The number of galactoses."""

    name: ClassVar = "nG"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of galactoses."""
        return _count_gal(glycan)


@register
@define
class CountGlcNAcMP(MetaProperty):
    """The number of GlcNAcs."""

    name: ClassVar = "nN"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of GlcNAcs."""
        return _count_glcnac(glycan)


@register
@define
class HasPolyLacNAcMP(MetaProperty):
    """Whether the glycan has any poly-LacNAc."""

    name: ClassVar = "Pl"
    supported_mode: ClassVar = "structure"
    return_type: ClassVar = "boolean"

    def calculate_one(self, glycan: Glycan) -> bool:
        """Calculate whether the glycan has any poly-LacNAc."""
        return _has_poly_lacnac(glycan)


@register
@define
class CountA23SiaMP(MetaProperty):
    """The number of sialic acids with an alpha-2,3 linkage."""

    name: ClassVar = "nL"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of sialic acids with an alpha-2,3 linkage."""
        return _count_a23_sia(glycan)


@register
@define
class CountA26SiaMP(MetaProperty):
    """The number of sialic acids with an alpha-2,6 linkage."""

    name: ClassVar = "nE"
    supported_mode: ClassVar = "both"
    return_type: ClassVar = "UInt8"

    def calculate_one(self, glycan: Glycan) -> int:
        """Calculate the number of sialic acids with an alpha-2,6 linkage."""
        return _count_a26_sia(glycan)


# ===== Cacheable functions for calculating meta-properties =====
@cache
def _get_branch_core_man(glycan: Structure) -> tuple[Monosaccharide, Monosaccharide]:
    # Branch core mannose is defined as:
    #
    # X - Man  <- This one
    #        \
    #         Man - Glc2NAc - Glc2NAc
    #        /
    # X - Man  <- And this one
    bft_iter = glycan.breadth_first_traversal(skip=["Fuc"])
    for i in range(3):
        next(bft_iter)
    while True:
        node1 = next(bft_iter)
        if get_mono_str(node1) == "Man":
            break
    while True:
        node2 = next(bft_iter)
        if get_mono_str(node2) == "Man":
            break
    return node1, node2


@cache
def _glycan_type(glycan: Structure) -> GlycanType:
    """Decide whether a glycan is a 'complex', 'hybrid', or 'high-mannose' type."""
    # N-glycan core is defined as the complex type.
    if glycan.composition == {"Glc2NAc": 2, "Man": 3}:
        # This type of glycan is actually pausimannose type.
        # However, it is currently regarded as a complex type right now.
        return GlycanType.COMPLEX

    # Bisecting could only be found in complex type.
    if _is_bisecting(glycan):
        return GlycanType.COMPLEX

    # If the glycan is not core, and it only has 2 "GlcNAc", it is high-mannose.
    if glycan.composition["Glc2NAc"] == 2:
        return GlycanType.HIGH_MANNOSE

    # If the glycan is mono-antennary and not high-monnose, it is complex.
    node1, node2 = _get_branch_core_man(glycan)
    if any((len(node1.links) == 1, len(node2.links) == 1)):
        return GlycanType.COMPLEX

    # Then, if it has 3 "Glc2NAc", it must be hybrid.
    if glycan.composition["Glc2NAc"] == 3:
        return GlycanType.HYBRID

    # All rest cases are complex.
    return GlycanType.COMPLEX


@cache
def _is_bisecting(glycan: Structure) -> bool:
    """Decide whether a glycan has bisection."""
    bft_iter = glycan.breadth_first_traversal(skip=["Fuc"])
    for i in range(2):
        next(bft_iter)
    next_node = next(bft_iter)
    return len(next_node.links) == 4


@cache
def _count_antenna(glycan: Structure) -> int:
    """Count the number of branches in a glycan."""
    if _glycan_type(glycan) != GlycanType.COMPLEX:
        return 0
    node1, node2 = _get_branch_core_man(glycan)
    return len(node1.links) + len(node2.links) - 2


@singledispatch
def _count_fuc(glycan) -> int:
    """Count the number of fucoses."""
    raise TypeError


@_count_fuc.register
@cache
def _(glycan: Structure) -> int:
    return glycan.composition.get("Fuc", 0)


@_count_fuc.register
@cache
def _(glycan: Composition) -> int:
    return glycan.get("F", 0)


@cache
def _get_core_ids(glycan: Structure) -> list[int]:
    """Get the IDs of the monosaccharides of the core."""
    should_be = ["Glc2NAc", "Glc2NAc", "Man", "Man", "Man"]
    cores: list[int] = []
    for node in glycan.breadth_first_traversal(skip=["Fuc"]):
        # The first two monosaacharides could only be "Glc2NAc".
        # However, when the glycan is bisecting, the rest of the three monosaccharides
        # might not all be "Man". So we look for the monosaacharides in the order of
        # `should_be`, and skip the ones that are not in the order.
        if get_mono_str(node) == should_be[len(cores)]:
            cores.append(node.id)
        if len(cores) == 5:
            break
    return cores


@cache
def _count_core_fuc(glycan: Structure) -> int:
    cores = _get_core_ids(glycan)
    n = 0
    for node in glycan.breadth_first_traversal():
        # node.parents()[0] is the nearest parent, and is a tuple of (Link, Monosaccharide)
        # node.parents()[0][1] is the Monosaccharide
        if get_mono_str(node) == "Fuc" and node.parents()[0][1].id in cores:
            n = n + 1
    return n


@singledispatch
def _count_sia(glycan) -> int:
    """Count the number of sialic acids."""
    raise TypeError


@_count_sia.register
@cache
def _(glycan: Structure) -> int:
    return glycan.composition.get("Neu5Ac", 0) + glycan.composition.get("Neu5Gc", 0)


@_count_sia.register
@cache
def _(glycan: Composition) -> int:
    return glycan.get("S", 0) + glycan.get("E", 0) + glycan.get("L", 0)


@singledispatch
def _count_gal(glycan) -> int:
    """Count the number of Gals."""
    raise TypeError


@_count_gal.register
@cache
def _(glycan: Structure) -> int:
    return glycan.composition.get("Gal", 0)


@_count_gal.register
@cache
def _(glycan: Composition) -> int:
    n_H = glycan.get("H", 0)
    n_N = glycan.get("N", 0)
    if n_H >= 4 and n_N >= n_H - 1:
        return n_H - 3
    return 0


@singledispatch
def _count_man(glycan) -> int:
    """Count the number of mannoses."""
    raise TypeError


@_count_man.register
@cache
def _(glycan: Structure) -> int:
    return glycan.composition.get("Man", 0)


@_count_man.register
@cache
def _(glycan: Composition) -> int:
    return glycan.get("H", 0) - _count_gal(glycan)


@singledispatch
def _count_glcnac(glycan) -> int:
    """Count the number of GlcNAcs."""
    raise TypeError


@_count_glcnac.register
@cache
def _(glycan: Structure) -> int:
    return glycan.composition.get("Glc2NAc", 0)


@_count_glcnac.register
@cache
def _(glycan: Composition) -> int:
    return glycan.get("N", 0)


@cache
def _has_poly_lacnac(glycan: Structure) -> bool:
    """Decide whether a glycan has any poly-LacNAc."""
    # Iterate all Gal residues.
    # If a Gal residue has a GlcNAc child, and the GlcNAc residue has a Gal child,
    # then the glycan has poly-LacNAc.
    gals = glycan.breadth_first_traversal(only=["Gal"])
    for gal in gals:
        for _, child in gal.children():
            if get_mono_str(child) == "Glc2NAc":
                children_of_glcnac = child.children()
                for _, child_of_glcnac in children_of_glcnac:
                    if get_mono_str(child_of_glcnac) == "Gal":
                        return True
    return False


@singledispatch
def _count_a23_sia(glycan) -> int:
    """Count the number of sialic acids with an alpha-2,3 linkage."""
    raise TypeError


@_count_a23_sia.register
@cache
def _(glycan: Structure) -> int:
    n = 0
    for node in glycan.breadth_first_traversal():
        if get_mono_str(node) in ("Neu5Ac", "Neu5Gc"):
            if node.links[2][0].parent_position == -1:
                raise SiaLinkageError("Sialic acid linkage not specified")
            elif node.links[2][0].parent_position == 3:
                n = n + 1
    return n


@_count_a23_sia.register
@cache
def _(glycan: Composition) -> int:
    return glycan.get("L", 0)


@singledispatch
def _count_a26_sia(glycan) -> int:
    """Count the number of sialic acids with an alpha-2,6 linkage."""
    raise TypeError


@_count_a26_sia.register
@cache
def _(glycan: Structure) -> int:
    n = 0
    for node in glycan.breadth_first_traversal():
        if get_mono_str(node) in ("Neu5Ac", "Neu5Gc"):
            if node.links[2][0].parent_position == -1:
                raise SiaLinkageError("Sialic acid linkage not specified")
            elif node.links[2][0].parent_position == 6:
                n = n + 1
    return n


@_count_a26_sia.register
@cache
def _(glycan: Composition) -> int:
    return glycan.get("E", 0)
