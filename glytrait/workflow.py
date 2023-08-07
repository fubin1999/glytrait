"""The main workflow of GlyTrait.

Functions:
    run_workflow: Run the whole GlyTrait workflow.
"""
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
    filter_invalid,
    filter_colinearity,
)

__all__ = ["Workflow"]


class Workflow:
    """The main workflow of GlyTrait.

    Attributes:
        config: Configuration for the workflow.

    Methods:
        run: Run the whole GlyTrait workflow.

    Examples:
        >>> config = Config()
        >>> config.update({"input_file": "data/abundance.csv"})
        >>> workflow = Workflow(config)
        >>> workflow.run()
    """

    def __init__(self, config: Config):
        self.config = config

    def run(self) -> None:
        """Run the whole GlyTrait workflow."""
        print('check')
        # Load data
        comps, strucs, abund_df = self._load_data()
        groups = self._load_groups()
        self.config.update({"_has_struc_col": strucs is not None})
        if len(abund_df.index) < 3:
            self.config.update({"post_filtering": False})

        # Preprocess and prepare data
        comps, strucs, abund_df = self._preprocess(comps, strucs, abund_df)
        glycans = self._load_glycans(comps, strucs)
        formulas = self._load_formulas()

        # Calculate meta properties and derived traits
        meta_prop_df = self._calcu_mate_properties(abund_df, glycans)
        derived_traits = self._calcu_derived_traits(abund_df, meta_prop_df, formulas)
        if self.config.get("post_filtering"):
            formulas, derived_traits = self._filter_traits(formulas, derived_traits)
        trait_table = self._combine_data(abund_df, derived_traits)

        # Statistical analysis
        if groups is not None:
            hypo_test_result, roc_result = self._statistical_analysis(trait_table, groups)
        else:
            hypo_test_result, roc_result = None, None

        # Write output
        write_output(
            self.config,
            derived_traits,
            abund_df,
            meta_prop_df,
            formulas,
            groups,
            hypo_test_result,
            roc_result,
        )

    def _load_data(self):
        """Load compositions, structures, and the abundance table."""
        return read_input(self.config.get("input_file"))

    def _preprocess(
        self, comps, strucs, abund_df
    ) -> tuple[list[str], list[str] | None, pd.DataFrame]:
        """Preprocess the data."""
        abund_df = preprocess_pipeline(
            abund_df,
            self.config.get("filter_glycan_max_na"),
            self.config.get("impute_method"),
        )
        to_keep = [c in abund_df.columns for c in comps]
        comps = [c for c, k in zip(comps, to_keep) if k]
        if strucs is not None:
            strucs = [s for s, k in zip(strucs, to_keep) if k]
        return comps, strucs, abund_df

    def _load_glycans(self, comps, strucs) -> list:
        """Load glycans."""
        if self.config.get("mode") == "structure":
            if strucs is None and self.config.get("structure_file") is not None:
                glycans = read_structure(self.config.get("structure_file"), comps)
            elif strucs is None and self.config.get("database") is not None:
                glycans = load_default_structures(self.config.get("database"), comps)
            elif strucs is not None:
                glycans = load_glycans(strucs)
            else:
                raise ValueError("Config error: no structure information.")
        else:
            glycans = load_compositions(
                comps, sia_linkage=self.config.get("sia_linkage")
            )
        return glycans

    def _load_formulas(self) -> list[TraitFormula]:
        """Load formulas."""
        formulas = load_formulas(self.config.get("mode"), self.config.get("formula_file"))
        if not self.config.get("sia_linkage"):
            formulas = [f for f in formulas if f.sia_linkage is False]
        return formulas

    def _calcu_mate_properties(self, abund_df, glycans) -> pd.DataFrame:
        """Calculate meta properties."""
        return build_meta_property_table(
            abund_df.columns,
            glycans,
            self.config.get("mode"),
            self.config.get("sia_linkage"),
        )

    def _calcu_derived_traits(self, abund_df, meta_prop_df, formulas) -> pd.DataFrame:
        """Calculate derived traits."""
        return calcu_derived_trait(abund_df, meta_prop_df, formulas)

    def _filter_traits(
        self, formulas: list[TraitFormula], derived_traits: pd.DataFrame,
    ) -> tuple[list[TraitFormula], pd.DataFrame]:
        """Filter invalid traits if needed."""
        formulas, derived_traits = filter_invalid(formulas, derived_traits)
        formulas, derived_traits = filter_colinearity(
            formulas,
            derived_traits,
            self.config.get("corr_threshold"),
            self.config.get("corr_method"),
        )
        return formulas, derived_traits

    def _combine_data(self, abund_df: pd.DataFrame, derived_traits: pd.DataFrame) -> pd.DataFrame:
        """Combine abundance table and derived traits."""
        return pd.concat([abund_df, derived_traits], axis=1)

    def _load_groups(self) -> pd.Series:
        """Get group information."""
        if self.config.get("group_file") is not None:
            groups = read_group(self.config.get("group_file"))
        else:
            groups = None
        return groups

    def _statistical_analysis(
        self, trait_table: pd.DataFrame, groups: pd.Series
    ) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
        """Perform statistical analysis."""
        if self.config.get("post_filtering") is False:
            print(
                "Warning: post-filtering is not enabled, "
                "so statistical analysis is not performed."
            )
            return None, None
        hypo_result = auto_hypothesis_test(trait_table, groups)
        if groups.nunique() == 2:
            roc_result = calcu_roc_auc(trait_table, groups)
        else:
            roc_result = None
        return hypo_result, roc_result
