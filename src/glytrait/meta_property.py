"""This module provides functions for calculating meta-properties of glycans.

Functions:
    build_meta_property_table: Build a table of meta-properties for glycans.
"""

from __future__ import annotations

from enum import Enum, auto
from functools import singledispatch, cache
from typing import Literal, ClassVar, Type

import pandas as pd
from attrs import define
from glypy import Monosaccharide  # type: ignore

from glytrait.exception import SiaLinkageError
from glytrait.glycan import Structure, Composition, get_mono_str, GlycanDict
from glytrait.data_type import MetaPropertyTable

__all__ = [
    "build_meta_property_table",
    "available_meta_properties",
]


# To be deleted
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
    for mp_name, MpClass in available_meta_properties(mode, sia_linkage).items():
        mp = MpClass()
        values = [mp(glycan) for glycan in glycans.values()]
        s = pd.Series(values, index=pd.Index(glycans.keys()), name=mp_name)
        mp_series_list.append(s)
    mp_table_df = pd.concat(mp_series_list, axis=1)
    return MetaPropertyTable(mp_table_df)


class GlycanType(Enum):
    """The type of glycan."""

    COMPLEX = auto()
    HIGH_MANNOSE = auto()
    HYBRID = auto()


# ===== Base classes for meta-properties =====
@define
class MetaProperty:
    """The base class for meta-properties.

    A meta-property is a function that takes a glycan as input and returns a value.
    The glycan can be either a `Composition` or a `Structure`.
    The value can be an integer, a boolean, or a string.

    For explicitness, all subclasses should be named with the suffix
    "MetaProperty" or "MP".
    """

    name: ClassVar[str]
    supported_mode: ClassVar[Literal["composition", "structure", "both"]]

    def __call__(self, glycan: Structure | Composition) -> int | bool | str:
        raise NotImplementedError


@define
class IntegerMetaProperty(MetaProperty):
    """Meta-property that returns an integer value."""

    def __call__(self, glycan: Structure | Composition) -> int:
        raise NotImplementedError


@define
class BooleanMetaProperty(MetaProperty):
    """Meta-property that returns a boolean value."""

    def __call__(self, glycan: Structure | Composition) -> bool:
        raise NotImplementedError


@define
class StringMetaProperty(MetaProperty):
    """Meta-property that returns a string value."""

    def __call__(self, glycan: Structure | Composition) -> str:
        raise NotImplementedError


# ===== Concrete meta-properties =====
@define
class GlycanTypeMP(StringMetaProperty):
    """The type of glycan."""

    name: ClassVar = "glycanType"
    supported_mode: ClassVar = "structure"

    def __call__(self, glycan: Structure) -> str:
        return _glycan_type(glycan).name.lower()


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


@define
class BisectionMP(BooleanMetaProperty):
    """Whether the glycan has bisection."""

    name: ClassVar = "bisection"
    supported_mode: ClassVar = "structure"

    def __call__(self, glycan: Structure) -> bool:
        return _is_bisecting(glycan)


@cache
def _is_bisecting(glycan: Structure) -> bool:
    """Decide whether a glycan has bisection."""
    bft_iter = glycan.breadth_first_traversal(skip=["Fuc"])
    for i in range(2):
        next(bft_iter)
    next_node = next(bft_iter)
    return len(next_node.links) == 4


@define
class CountAntennaMP(IntegerMetaProperty):
    """The number of antennas."""

    name: ClassVar = "nAnt"
    supported_mode: ClassVar = "structure"

    def __call__(self, glycan: Structure) -> int:
        # Calling `_count_antenna` instead of putting the implementation here is for the purpose
        # of caching.
        # Methods decorated with `@cache` will not work properly due to the `self` argument.
        # This is a workaround and a common practice in this module.
        return _count_antenna(glycan)


@cache
def _count_antenna(glycan: Structure) -> int:
    """Count the number of branches in a glycan."""
    if _glycan_type(glycan) != GlycanType.COMPLEX:
        return 0
    node1, node2 = _get_branch_core_man(glycan)
    return len(node1.links) + len(node2.links) - 2


@define
class CountFucMP(IntegerMetaProperty):
    """The number of fucoses."""

    name: ClassVar = "nF"
    supported_mode: ClassVar = "both"

    def __call__(self, glycan: Structure | Composition) -> int:
        return _count_fuc(glycan)


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


@define
class CountCoreFucMP(IntegerMetaProperty):
    """The number of fucoses on the core."""

    name: ClassVar = "nFc"
    supported_mode: ClassVar = "structure"

    def __call__(self, glycan: Structure) -> int:
        return _count_core_fuc(glycan)


@define
class CountAntennaryFucMP(IntegerMetaProperty):
    """The number of fucoses on the antenna."""

    name: ClassVar = "nFa"
    supported_mode: ClassVar = "structure"

    def __call__(self, glycan: Structure) -> int:
        return _count_fuc(glycan) - _count_core_fuc(glycan)


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


@define
class CountSiaMP(IntegerMetaProperty):
    """The number of sialic acids."""

    name: ClassVar = "nS"
    supported_mode: ClassVar = "both"

    def __call__(self, glycan: Structure | Composition) -> int:
        return _count_sia(glycan)


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


@define
class CountManMP(IntegerMetaProperty):
    """The number of mannoses."""

    name: ClassVar = "nM"
    supported_mode: ClassVar = "both"

    def __call__(self, glycan: Structure | Composition) -> int:
        return _count_man(glycan)


@define
class CountGalMP(IntegerMetaProperty):
    """The number of galactoses."""

    name: ClassVar = "nG"
    supported_mode: ClassVar = "both"

    def __call__(self, glycan: Structure | Composition) -> int:
        return _count_gal(glycan)


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


# ===== The old module =====


struc_meta_properties: list[Type[MetaProperty]] = []
comp_meta_properties: list[Type[MetaProperty]] = []


def register_struc(meta_property: Type[MetaProperty]):
    """Register a structure meta-property."""
    struc_meta_properties.append(meta_property)
    return meta_property


def register_comp(meta_property: Type[MetaProperty]):
    """Register a composition meta-property."""
    comp_meta_properties.append(meta_property)
    return meta_property


def available_meta_properties(
    mode: Literal["composition", "structure"],
    sia_linkage: bool,
    only_sia_linkage: bool = False,
) -> dict[str, Type[MetaProperty]]:
    """Get the names of the available meta properties.

    Args:
        mode (Literal["composition", "structure"]): The calculation mode.
        sia_linkage (bool): Whether to include the sialic acid linkage meta-properties.
        only_sia_linkage (bool, optional): Whether to only include the sialic acid linkage.
            meta-properties. Defaults to False.
            `sia_linkage` must be True if this is True.

    Returns:
        dict[str, Type[MetaProperty]]: The mapping from the names to the meta-property classes.
    """
    if only_sia_linkage and not sia_linkage:
        raise ValueError("sia_linkage must be True if only_sia_linkage is True")

    match mode:
        case "composition":
            mp_list = comp_meta_properties
        case "structure":
            mp_list = struc_meta_properties
        case _:
            raise ValueError(f"Invalid mode: {mode}")

    if only_sia_linkage:
        return {mp.name: mp for mp in mp_list if mp.sia_linkage}
    elif sia_linkage:
        return {mp.name: mp for mp in mp_list}
    else:
        return {mp.name: mp for mp in mp_list if not mp.sia_linkage}


@register_comp
@register_struc
@define
class Wildcard(MetaProperty):
    """This meta-property has value 1 for all glycans."""

    name: ClassVar = "."
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Composition | Structure) -> float:
        return 1.0


@register_struc
@define
class IsComplex(MetaProperty):
    """Whether the glycan is a complex type."""

    name: ClassVar = "isComplex"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _glycan_type(glycan) == GlycanType.COMPLEX


@register_struc
@define
class IsHighMannose(MetaProperty):
    """Whether the glycan is a high-mannose type."""

    name: ClassVar = "isHighMannose"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _glycan_type(glycan) == GlycanType.HIGH_MANNOSE


@register_struc
@define
class IsHybrid(MetaProperty):
    """Whether the glycan is a hybrid type."""

    name: ClassVar = "isHybrid"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _glycan_type(glycan) == GlycanType.HYBRID


@register_struc
@define
class IsBisecting(MetaProperty):
    """Whether the glycan has bisection."""

    name: ClassVar = "isBisecting"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _is_bisecting(glycan)


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


@register_struc
@define
class Is1Antennary(MetaProperty):
    """Whether the glycan has one antenna."""

    name: ClassVar = "is1Antennary"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _count_antenna(glycan) == 1


@register_struc
@define
class Is2Antennary(MetaProperty):
    """Whether the glycan has two antennas."""

    name: ClassVar = "is2Antennary"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _count_antenna(glycan) == 2


@register_struc
@define
class Is3Antennary(MetaProperty):
    """Whether the glycan has three antennas."""

    name: ClassVar = "is3Antennary"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _count_antenna(glycan) == 3


@register_struc
@define
class Is4Antennary(MetaProperty):
    """Whether the glycan has four antennas."""

    name: ClassVar = "is4Antennary"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _count_antenna(glycan) == 4


@register_struc
@define
class TotalAntenna(MetaProperty):
    """The total number of antennas."""

    name: ClassVar = "totalAntenna"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _count_antenna(glycan)


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


@register_struc
@define
class CoreFuc(MetaProperty):
    """The number of fucoses on the core."""

    name: ClassVar = "coreFuc"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return _count_core_fuc(glycan)


@register_struc
@define
class AntennaryFuc(MetaProperty):
    """The number of fucoses on the antenna."""

    name: ClassVar = "antennaryFuc"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
        return glycan._composition.get("Fuc", 0) - _count_core_fuc(glycan)


@register_struc
@register_comp
@define
class TotalFuc(MetaProperty):
    """Total number of fucoses."""

    name: ClassVar = "totalFuc"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_fuc(glycan)


@register_struc
@register_comp
@define
class HasFuc(MetaProperty):
    """Whether the glycan has any fucoses."""

    name: ClassVar = "hasFuc"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_fuc(glycan) > 0


@register_struc
@register_comp
@define
class NoFuc(MetaProperty):
    """Whether the glycan has no fucoses."""

    name: ClassVar = "noFuc"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_fuc(glycan) == 0


@register_struc
@register_comp
@define
class TotalSia(MetaProperty):
    """The total number of sialic acids."""

    name: ClassVar = "totalSia"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_sia(glycan)


@register_struc
@register_comp
@define
class HasSia(MetaProperty):
    """Whether the glycan has any sialic acids."""

    name: ClassVar = "hasSia"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_sia(glycan) > 0


@register_struc
@register_comp
@define
class NoSia(MetaProperty):
    """Whether the glycan has no sialic acids."""

    name: ClassVar = "noSia"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_sia(glycan) == 0


@register_struc
@register_comp
@define
class TotalMan(MetaProperty):
    """The total number of mannoses."""

    name: ClassVar = "totalMan"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_man(glycan)


@register_struc
@register_comp
@define
class TotalGal(MetaProperty):
    """The total number of galactoses."""

    name: ClassVar = "totalGal"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_gal(glycan)


@register_struc
@define
class HasPolyLacNAc(MetaProperty):
    """Whether the glycan has any poly-LacNAc."""

    name: ClassVar = "hasPolyLacNAc"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Structure) -> float:
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


@register_struc
@register_comp
@define
class A23Sia(MetaProperty):
    """The number of sialic acids with an alpha-2,3 linkage."""

    name: ClassVar = "a23Sia"
    sia_linkage: ClassVar = True

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_a23_sia(glycan)


@register_struc
@register_comp
@define
class A26Sia(MetaProperty):
    """The number of sialic acids with an alpha-2,6 linkage."""

    name: ClassVar = "a26Sia"
    sia_linkage: ClassVar = True

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_a26_sia(glycan)


@register_struc
@register_comp
@define
class HasA23Sia(MetaProperty):
    """Whether the glycan has any sialic acids with an alpha-2,3 linkage."""

    name: ClassVar = "hasa23Sia"
    sia_linkage: ClassVar = True

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_a23_sia(glycan) > 0


@register_struc
@register_comp
@define
class HasA26Sia(MetaProperty):
    """Whether the glycan has any sialic acids with an alpha-2,6 linkage."""

    name: ClassVar = "hasa26Sia"
    sia_linkage: ClassVar = True

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_a26_sia(glycan) > 0


@register_struc
@register_comp
@define
class NoA23Sia(MetaProperty):
    """Whether the glycan has no sialic acids with an alpha-2,3 linkage."""

    name: ClassVar = "noa23Sia"
    sia_linkage: ClassVar = True

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_a23_sia(glycan) == 0


@register_struc
@register_comp
@define
class NoA26Sia(MetaProperty):
    """Whether the glycan has no sialic acids with an alpha-2,6 linkage."""

    name: ClassVar = "noa26Sia"
    sia_linkage: ClassVar = True

    def __call__(self, glycan: Structure | Composition) -> float:
        return _count_a26_sia(glycan) == 0


@register_comp
@define
class IsHighBranching(MetaProperty):
    """Whether the glycan has a high branching."""

    name: ClassVar = "isHighBranching"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Composition) -> float:
        return glycan.get("N", 0) > 4


@register_comp
@define
class IsLowBranching(MetaProperty):
    """Whether the glycan has a low branching."""

    name: ClassVar = "isLowBranching"
    sia_linkage: ClassVar = False

    def __call__(self, glycan: Composition) -> float:
        return glycan.get("N", 0) <= 4
