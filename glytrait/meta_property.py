from typing import Iterable, Literal

import pandas as pd

from glytrait.glycan import NGlycan, Composition

basic_struc_meta_properties = {
    ".",
    "isComplex",
    "isHighMannose",
    "isHybrid",
    "isBisecting",
    "is1Antennary",
    "is2Antennary",
    "is3Antennary",
    "is4Antennary",
    "totalAntenna",
    "coreFuc",
    "antennaryFuc",
    "hasAntennaryFuc",
    "totalFuc",
    "hasFuc",
    "noFuc",
    "totalSia",
    "hasSia",
    "noSia",
    "totalMan",
    "totalGal",
}
sia_struc_meta_properties = {
    "a23Sia",
    "a26Sia",
    "hasa23Sia",
    "hasa26Sia",
    "noa23Sia",
    "noa26Sia",
}
struc_meta_properties = basic_struc_meta_properties | sia_struc_meta_properties
basic_comp_meta_properties = {
    ".",
    "isHighBranching",
    "isLowBranching",
    "totalSia",
    "totalFuc",
    "totalGal",
    "hasSia",
    "hasFuc",
    "hasGal",
    "noSia",
    "noFuc",
    "noGal",
}
sia_comp_meta_properties = {
    "a23Sia",
    "a26Sia",
    "hasa23Sia",
    "hasa26Sia",
    "noa23Sia",
    "noa26Sia",
}
comp_meta_properties = basic_comp_meta_properties | sia_comp_meta_properties


def build_meta_property_table(
    glycan_ids: Iterable[str],
    glycans: Iterable[NGlycan | Composition],
    mode: Literal["composition", "structure"],
    sia_linkage: bool = False,
) -> pd.DataFrame:
    """Build a table of meta properties for glycans.

    Args:
        glycan_ids (Iterable[str]): The IDs of the glycans.
        glycans (Iterable[NGlycan]): The glycans.
        mode (Literal["composition", "structure"]): The calculation mode.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta properties. Defaults to False.

    Returns:
        pd.DataFrame: The table of meta properties, with Compositions as the index,
    """
    if mode == "composition":
        return _build_comp_meta_property_table(glycan_ids, glycans, sia_linkage)
    elif mode == "structure":
        return _build_struc_meta_property_table(glycan_ids, glycans, sia_linkage)
    else:
        raise ValueError(f"Invalid type: {mode}")


def _build_struc_meta_property_table(
    glycan_ids: Iterable[str], glycans: Iterable[NGlycan], sia_linkage: bool = False
) -> pd.DataFrame:
    """Build a table of structural meta properties for glycans.

    The following meta properties are included:
        - isComplex: Whether the glycan is a complex type.
        - isHighMannose: Whether the glycan is a high-mannose type.
        - isHybrid: Whether the glycan is a hybrid type.
        - isBisecting: Whether the glycan has a bisection.
        - is1Antennary: Whether the glycan has a 1-antenna.
        - is2Antennary: Whether the glycan has a 2-antenna.
        - is3Antennary: Whether the glycan has a 3-antenna.
        - is4Antennary: Whether the glycan has a 4-antenna.
        - totalAntenna: The total number of antennae.
        - coreFuc: The number of fucoses on the core.
        - antennaryFuc: The number of fucoses on the antenna.
        - hasAntennaryFuc: Whether the glycan has any fucoses on the antenna.
        - totalFuc: The total number of fucoses.
        - hasFuc: Whether the glycan has any fucoses.
        - noFuc: Whether the glycan has no fucoses.
        - totalSia: The total number of sialic acids.
        - hasSia: Whether the glycan has any sialic acids.
        - noSia: Whether the glycan has no sialic acids.
        - totalMan: The total number of mannoses.
        - totalGal: The total number of galactoses.

    If `sia_linkage` is True, the following meta properties are also included:
        - a23Sia: The number of sialic acids with an alpha-2,3 linkage.
        - a26Sia: The number of sialic acids with an alpha-2,6 linkage.
        - hasa23Sia: Whether the glycan has any sialic acids with an alpha-2,3 linkage.
        - hasa26Sia: Whether the glycan has any sialic acids with an alpha-2,6 linkage.
        - noa23Sia: Whether the glycan has no sialic acids with an alpha-2,3 linkage.
        - noa26Sia: Whether the glycan has no sialic acids with an alpha-2,6 linkage.

    Args:
        glycan_ids (Iterable[str]): The IDs of the glycans. Could be the
            accession, or the composition.
        glycans (Iterable[NGlycan]): The glycans.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta properties. Defaults to False.

    Returns:
        pd.DataFrame: The table of meta properties, with Compositions as the index,
            and the property names as the columns.
    """
    meta_property_table = pd.DataFrame(index=list(glycan_ids))

    meta_property_table["isComplex"] = [g.is_complex() for g in glycans]
    meta_property_table["isHighMannose"] = [g.is_high_mannose() for g in glycans]
    meta_property_table["isHybrid"] = [g.is_hybrid() for g in glycans]
    meta_property_table["isBisecting"] = [g.is_bisecting() for g in glycans]
    meta_property_table["is1Antennary"] = [g.count_antenna() == 1 for g in glycans]
    meta_property_table["is2Antennary"] = [g.count_antenna() == 2 for g in glycans]
    meta_property_table["is3Antennary"] = [g.count_antenna() == 3 for g in glycans]
    meta_property_table["is4Antennary"] = [g.count_antenna() == 4 for g in glycans]
    meta_property_table["totalAntenna"] = [g.count_antenna() for g in glycans]
    meta_property_table["coreFuc"] = [g.count_core_fuc() for g in glycans]
    meta_property_table["antennaryFuc"] = [g.count_antennary_fuc() for g in glycans]
    meta_property_table["hasAntennaryFuc"] = [
        g.count_antennary_fuc() > 0 for g in glycans
    ]
    meta_property_table["totalFuc"] = [g.count_fuc() for g in glycans]
    meta_property_table["hasFuc"] = [g.count_fuc() > 0 for g in glycans]
    meta_property_table["noFuc"] = [g.count_fuc() == 0 for g in glycans]
    meta_property_table["totalSia"] = [g.count_sia() for g in glycans]
    meta_property_table["hasSia"] = [g.count_sia() > 0 for g in glycans]
    meta_property_table["noSia"] = [g.count_sia() == 0 for g in glycans]
    meta_property_table["totalMan"] = [g.count_man() for g in glycans]
    meta_property_table["totalGal"] = [g.count_gal() for g in glycans]

    if sia_linkage:
        meta_property_table["a23Sia"] = [g.count_a23_sia() for g in glycans]
        meta_property_table["a26Sia"] = [g.count_a26_sia() for g in glycans]
        meta_property_table["hasa23Sia"] = [g.count_a23_sia() > 0 for g in glycans]
        meta_property_table["hasa26Sia"] = [g.count_a26_sia() > 0 for g in glycans]
        meta_property_table["noa23Sia"] = [g.count_a23_sia() == 0 for g in glycans]
        meta_property_table["noa26Sia"] = [g.count_a26_sia() == 0 for g in glycans]

    return meta_property_table


def _build_comp_meta_property_table(
    glycan_ids: Iterable[str], glycans: Iterable[Composition], sia_linkage: bool = False
) -> pd.DataFrame:
    """Build a table of compositional meta properties for glycans.

    The following meta properties are included:
        - isHighBranching: Whether the glycan has a high branching.
        - isLowBranching: Whether the glycan has a low branching.
        - totalSia: The total number of sialic acids.
        - totalFuc: The total number of fucoses.
        - totalGal: The total number of galactoses.
        - hasSia: Whether the glycan has any sialic acids.
        - hasFuc: Whether the glycan has any fucoses.
        - hasGal: Whether the glycan has any galactoses.
        - noSia: Whether the glycan has no sialic acids.
        - noFuc: Whether the glycan has no fucoses.
        - noGal: Whether the glycan has no galactoses.

    If `sia_linkage` is True, the following meta properties are also included:
        - a23Sia: The number of sialic acids with an alpha-2,3 linkage.
        - a26Sia: The number of sialic acids with an alpha-2,6 linkage.
        - hasa23Sia: Whether the glycan has any sialic acids with an alpha-2,3 linkage.
        - hasa26Sia: Whether the glycan has any sialic acids with an alpha-2,6 linkage.
        - noa23Sia: Whether the glycan has no sialic acids with an alpha-2,3 linkage.
        - noa26Sia: Whether the glycan has no sialic acids with an alpha-2,6 linkage.

    Args:
        glycan_ids (Iterable[str]): The IDs of the glycans. Could be the
            accession, or the composition.
        glycans (Iterable[Composition]): The glycans.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta properties. Defaults to False.

    Returns:
        pd.DataFrame: The table of meta properties, with Compositions as the index,
            and the property names as the columns.
    """
    meta_property_table = pd.DataFrame(index=list(glycan_ids))

    meta_property_table["isHighBranching"] = [g.is_high_branching() for g in glycans]
    meta_property_table["isLowBranching"] = [g.is_low_branching() for g in glycans]
    meta_property_table["totalSia"] = [g.count_sia() for g in glycans]
    meta_property_table["totalFuc"] = [g.count_fuc() for g in glycans]
    meta_property_table["totalGal"] = [g.count_gal() for g in glycans]
    meta_property_table["hasSia"] = [g.count_sia() > 0 for g in glycans]
    meta_property_table["hasFuc"] = [g.count_fuc() > 0 for g in glycans]
    meta_property_table["hasGal"] = [g.count_gal() > 0 for g in glycans]
    meta_property_table["noSia"] = [g.count_sia() == 0 for g in glycans]
    meta_property_table["noFuc"] = [g.count_fuc() == 0 for g in glycans]
    meta_property_table["noGal"] = [g.count_gal() == 0 for g in glycans]

    if sia_linkage:
        meta_property_table["a23Sia"] = [g.count_a23_sia() for g in glycans]
        meta_property_table["a26Sia"] = [g.count_a26_sia() for g in glycans]
        meta_property_table["hasa23Sia"] = [g.count_a23_sia() > 0 for g in glycans]
        meta_property_table["hasa26Sia"] = [g.count_a26_sia() > 0 for g in glycans]
        meta_property_table["noa23Sia"] = [g.count_a23_sia() == 0 for g in glycans]
        meta_property_table["noa26Sia"] = [g.count_a26_sia() == 0 for g in glycans]

    return meta_property_table
