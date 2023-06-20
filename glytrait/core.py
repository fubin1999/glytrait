import pandas as pd

from glytrait.config import Config
from glytrait.database import load_default
from glytrait.glycan import load_glycans, load_compositions
from glytrait.io import read_input, write_output, read_group, read_structure
from glytrait.preprocessing import preprocess_pipeline
from glytrait.stats import auto_hypothesis_test
from glytrait.trait import (
    load_formulas,
    build_meta_property_table,
    calcu_derived_trait,
    filter_derived_trait,
    TraitFormula,
)

__all__ = ["run_workflow"]


# --------------------------------------------------
# Helper functions directly called by `run_workflow`
# --------------------------------------------------


def _load_and_preprocess_data(config: Config) -> tuple[list, pd.DataFrame]:
    """Load glycans and the abundance table."""
    comps, strucs, abund_df = read_input(config.get("input_file"))
    config.update({"_has_struc_col": strucs is not None})
    if config.get("mode") == "structure":
        if strucs is None and config.get("structure_file") is not None:
            glycans = read_structure(config.get("structure_file"), comps)
        elif strucs is None and config.get("database") is not None:
            glycans = load_default(config.get("database"), comps)
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
        derived_traits = filter_derived_trait(derived_traits)
        formulas = [f for f in formulas if f.name in derived_traits.columns]
    return derived_traits, formulas


def _perform_hypothesis_test(
    config: Config,
    abund_df: pd.DataFrame,
    derived_traits: pd.DataFrame,
) -> tuple[pd.DataFrame | None, pd.Series | None]:
    """Perform hypothesis test."""
    if config.get("group_file") is not None:
        groups = read_group(config.get("group_file"))
        combined_traits = pd.concat([abund_df, derived_traits], axis=1)
        hypo_result = auto_hypothesis_test(combined_traits, groups)
    else:
        groups = None
        hypo_result = None
    return hypo_result, groups


# ---------------------------------------------------------
# End of helper functions directly called by `run_workflow`
# ---------------------------------------------------------


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
    hypo_result, groups = _perform_hypothesis_test(config, abund_df, derived_traits)
    write_output(
        config, derived_traits, abund_df, meta_prop_df, formulas, groups, hypo_result
    )
