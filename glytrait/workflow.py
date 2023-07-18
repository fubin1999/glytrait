import pandas as pd

from glytrait.analysis import auto_hypothesis_test, calcu_roc_auc
from glytrait.config import Config
from glytrait.formula import TraitFormula, load_formulas
from glytrait.glycan import load_glycans, load_compositions
from glytrait.io import (
    read_input,
    write_output,
    read_group,
    read_structure,
    load_default_structures,
)
from glytrait.meta_property import build_meta_property_table
from glytrait.preprocessing import preprocess_pipeline
from glytrait.trait import (
    calcu_derived_trait,
    filter_invalid, filter_colinearity,
)

__all__ = ["run_workflow"]


def run_workflow(config: Config) -> None:
    """Run the whole GlyTrait workflow.

    Args:
        config: Configuration for the workflow.
    """
    glycans, abund_df = _load_and_preprocess_data(config)
    formulas = _load_formulas(config)
    meta_prop_df, derived_traits = _calcu_derived_traits(
        config, abund_df, glycans, formulas
    )
    derived_traits, formulas = _filter_invalid_traits(config, derived_traits, formulas)
    important_traits = _important_trait_table(config, derived_traits, formulas)
    trait_table = _combine_data(abund_df, derived_traits)
    groups, hypo_test_result, roc_result = _statistical_analysis(config, trait_table)
    write_output(
        config,
        derived_traits,
        abund_df,
        important_traits,
        meta_prop_df,
        formulas,
        groups,
        hypo_test_result,
        roc_result,
    )


def _load_and_preprocess_data(config: Config) -> tuple[list, pd.DataFrame]:
    """Load glycans and the abundance table."""
    comps, strucs, abund_df = read_input(config.get("input_file"))
    config.update({"_has_struc_col": strucs is not None})
    if config.get("mode") == "structure":
        if strucs is None and config.get("structure_file") is not None:
            glycans = read_structure(config.get("structure_file"), comps)
        elif strucs is None and config.get("database") is not None:
            glycans = load_default_structures(config.get("database"), comps)
        elif strucs is not None:
            glycans = load_glycans(strucs)
        else:
            raise ValueError("Config error: no structure information.")
    else:
        glycans = load_compositions(comps, sia_linkage=config.get("sia_linkage"))

    abund_df = preprocess_pipeline(
        abund_df, config.get("filter_glycan_max_na"), config.get("impute_method")
    )
    return glycans, abund_df


def _load_formulas(config: Config) -> list[TraitFormula]:
    """Load formulas."""
    formulas = load_formulas(config.get("mode"), config.get("formula_file"))
    if not config.get("sia_linkage"):
        formulas = [f for f in formulas if f.sia_linkage is False]
    return formulas


def _calcu_derived_traits(
    config: Config, abund_df: pd.DataFrame, glycans: list, formulas: list[TraitFormula]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate derived traits."""
    meta_prop_df = build_meta_property_table(
        abund_df.columns,
        glycans,
        config.get("mode"),
        config.get("sia_linkage"),
    )
    derived_traits = calcu_derived_trait(abund_df, meta_prop_df, formulas)
    return meta_prop_df, derived_traits


def _filter_invalid_traits(
    config: Config, derived_traits: pd.DataFrame, formulas: list[TraitFormula]
) -> tuple[pd.DataFrame, list[TraitFormula]]:
    """Filter invalid traits if needed."""
    if config.get("filter_invalid_traits"):
        derived_traits = filter_invalid(derived_traits)
        formulas = [f for f in formulas if f.name in derived_traits.columns]
    return derived_traits, formulas


def _get_groups(config: Config) -> pd.Series:
    """Get group information."""
    if config.get("group_file") is not None:
        groups = read_group(config.get("group_file"))
    else:
        groups = None
    return groups


def _combine_data(abund_df: pd.DataFrame, derived_traits: pd.DataFrame) -> pd.DataFrame:
    """Combine abundance table and derived traits."""
    return pd.concat([abund_df, derived_traits], axis=1)


def _important_trait_table(
    config: Config, derived_trait_df: pd.DataFrame, formulas: list[TraitFormula]
) -> pd.DataFrame:
    """Get important trait table."""
    important = filter_colinearity(
        formulas, derived_trait_df, threshold=config.get("correlation_threshold")
    )
    important_trait = derived_trait_df.columns[important]
    important_trait_df = derived_trait_df[important_trait]
    return important_trait_df


def _statistical_analysis(
    config: Config, trait_table: pd.DataFrame
) -> tuple[pd.Series | None, pd.DataFrame | None, pd.DataFrame | None]:
    """Perform statistical analysis."""
    if (groups := _get_groups(config)) is not None:
        hypo_result = auto_hypothesis_test(trait_table, groups)
        if groups.nunique() == 2:
            roc_result = calcu_roc_auc(trait_table, groups)
        else:
            roc_result = None
    else:
        hypo_result = None
        roc_result = None
    return groups, hypo_result, roc_result
