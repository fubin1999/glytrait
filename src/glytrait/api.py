"""The main API of GlyTrait.

This module provides the `Experiment` class, which is the main entry point
to perform the GlyTrait workflow.

Example:
    >>> from glytrait.api import Experiment
    >>> experiment = Experiment(
    ...    abundance_file="glycan_abundance.csv",
    ...    glycan_file="glycan_structure.csv",
    ...    group_file="group.csv",
    ...    mode="structure",
    ... )
    >>> experiment.run_workflow()
    >>> experiment.filtered_derived_trait_table
    # the filtered derived trait table
"""

from collections import defaultdict
from collections.abc import Iterable
from typing import cast, Literal, Optional

import pandas as pd
from attrs import define, field

from glytrait.data_type import (
    MetaPropertyTable,
    DerivedTraitTable,
    GroupSeries,
    AbundanceTable,
)
from glytrait.exception import GlyTraitError
from glytrait.formula import load_default_formulas, TraitFormula, parse_formulas
from glytrait.data_input import (
    load_glycans,
    load_groups,
    load_abundance,
    check_all_glycans_have_struct_or_comp,
    check_same_samples_in_abund_and_groups,
)
from glytrait.meta_property import build_meta_property_table
from glytrait.post_filtering import post_filter
from glytrait.preprocessing import preprocess
from glytrait.trait import calcu_derived_trait
from glytrait.stat import auto_test
from glytrait.glycan import Structure, Composition


class ExperimentError(GlyTraitError):
    """Base class for exceptions in the Experiment class."""


class InvalidOperationOrderError(ExperimentError):
    """Raised when calling a method in an invalid order."""


class MissingDataError(ExperimentError):
    """Raised when some data is missing for the operation."""


GlycanDict = dict[str, Structure] | dict[str, Composition]


@define
class Experiment:
    """GlyTrait experiment.

    Create an instance of this class to perform all the steps in the GlyTrait workflow.
    After the instance is created, the simplest way to run the workflow is to call
    `run_workflow` with the desired parameters.
    This will call the methods in correct order and store the results as attributes.
    If you want a more flexible way to run the workflow, you can call the methods
    mannually in the correct order.
    The detailed steps are described below.

    Firstly, call the `preprocessed` method to filter glycans, impute missing values,
    and normalize the abundance.
    Calling this method makes the `processed_abundance_table` attribute available.

    Secondly, call the `derive_traits` method to calculate all the derived traits.
    The result table is stored as the `derived_trait_table` attribute.

    Thirdly, call the `post_filter` method
    to remove invalid traits and highly correlated traits.
    The result table is stored as the `filtered_derived_trait_table` attribute.

    Finally, and optionally, call the `diff_analysis` method with a group series
    to perform differential analysis.
    The result is stored as the `diff_results` attribute.

    Note that you could always re-run any workflow step to update the intermediate results
    and clean up the following results.
    For example, if you call `preprocess` after `derive_traits`,
    the `processed_abundance_table` will be updated,
    and the `derived_trait_table` will be cleared.
    This is convenient when you want to try different parameters or formulas.

    This class also provides a method `try_formulas` to try new formula expressions.
    First, call the `preprocess`, then call `try_formulas` with the formula expressions to try.
    Instead of saving the result as an attribute,
    this method returns the result DataFrame or Series directly
    (depending on the number of formula expressions).

    Methods:
        preprocess: Preprocess the data.
        derive_traits: Calculate derived traits.
        post_filter: Post-filter the derived traits.
        diff_analysis: Perform differential analysis.
        derive_one_trait: Try a single trait formula.
        reset: Reset the workflow to the beginning (before `preprocess`).
        get_data: Get the data by name.

    Attributes:
        abundance_file: The path to the abundance file.
        glycan_file: The path to the glycan file.
        group_file: The path to the group file.
        mode: "structure" or "composition".
        sia_linkage: Whether to consider the linkage of sialic acid.
        abundance_table: The original abundance table.
        groups: The group series. Available if `group_file` is provided.
        processed_abundance_table: The processed abundance table. Available after `preprocess`.
        meta_property_table: The meta property table. Either provided or extracted from
            `glycan_file`.
        derived_trait_table: The derived trait table. Available after `derive_traits`.
        filtered_derived_trait_table: The filtered derived trait table.
            Available after `post_filter`.
        diff_results: The differential analysis results. Available after `diff_analysis`.

    Examples:
        >>> from glytrait.api import Experiment
        >>> experiment = Experiment(
        ...    abundance_file="glycan_abundance.csv",
        ...    glycan_file="glycan_structure.csv",
        ...    group_file="group.csv",
        ...    mode="structure",
        ... )

        # Run the entire workflow
        >>> experiment.run_workflow()

        # Or run the workflow step by step
        >>> experiment.preprocess(filter_max_na=0.5, impute_method="min")
        >>> experiment.derive_traits()  # with default formulas
        >>> experiment.post_filter(corr_threshold=0.9)
        >>> experiment.diff_analysis()

        # Get results
        >>> experiment.filtered_derived_trait_table
        >>> experiment.diff_results

        # Try new formulas
        >>> experiment.try_formulas("TC = [type == 'complex'] / [1]")
    """

    abundance_file: Optional[str] = field(repr=False)
    glycan_file: Optional[str] = field(repr=False)
    group_file: Optional[str] = field(default=None, kw_only=True, repr=False)
    mode: Literal["structure", "composition"] = field(default="structure", kw_only=True)
    sia_linkage: bool = field(default=False, kw_only=True)

    # The following attributes are generated during the workflow
    _abundance_table: AbundanceTable = field(init=False, default=None)
    _groups: GroupSeries | None = field(init=False, default=None)
    _meta_property_table: MetaPropertyTable = field(init=False, default=None)
    _processed_abundance_table: AbundanceTable = field(init=False, default=None)
    _formulas: list[TraitFormula] = field(init=False, default=None)
    _derived_trait_table: DerivedTraitTable = field(init=False, default=None)
    _filtered_derived_trait_table: DerivedTraitTable = field(init=False, default=None)
    _diff_results: dict[str, pd.DataFrame] = field(init=False, default=None)

    def __attrs_post_init__(self):
        abundance_type = defaultdict(lambda: "float64")
        abundance_type.update({"Sample": "O"})
        abund_df = load_abundance(
            pd.read_csv(self.abundance_file, dtype=abundance_type)
        )
        self._abundance_table = abund_df

        glycans = load_glycans(pd.read_csv(self.glycan_file), mode=self.mode)
        check_all_glycans_have_struct_or_comp(abund_df, glycans)
        self._meta_property_table = self._extract_meta_properties(glycans)

        if self.group_file:
            groups = load_groups(pd.read_csv(self.group_file))
            check_same_samples_in_abund_and_groups(abund_df, groups)
            self._groups = groups

    def _extract_meta_properties(self, glycans: GlycanDict) -> MetaPropertyTable:
        """Extract meta-properties from glycan structures or compositions."""
        glycan_names: list[str] = self._abundance_table.columns.tolist()
        glycan_dict = cast(GlycanDict, {g: glycans[g] for g in glycan_names})
        return build_meta_property_table(glycan_dict, self.mode, self.sia_linkage)

    # ===== Properties of the data ===== START
    @property
    def abundance_table(self) -> AbundanceTable:
        """The original abundance table."""
        return AbundanceTable(self._abundance_table.copy())

    @property
    def groups(self) -> GroupSeries | None:
        """The group series."""
        return GroupSeries(self._groups.copy()) if self._groups is not None else None

    @property
    def meta_property_table(self) -> MetaPropertyTable:
        """The meta property table."""
        if self._meta_property_table is not None:
            return MetaPropertyTable(self._meta_property_table)
        msg = (
            "Meta property table is not available. "
            "Please call the `preprocess` method first."
        )
        raise MissingDataError(msg)

    @property
    def processed_abundance_table(self) -> AbundanceTable:
        """The processed abundance table."""
        if self._processed_abundance_table is not None:
            return self._processed_abundance_table
        msg = "Please call the `preprocess` method first."
        raise MissingDataError(msg)

    @property
    def derived_trait_table(self) -> DerivedTraitTable:
        """The derived trait table."""
        if self._derived_trait_table is not None:
            return self._derived_trait_table
        msg = "Please call the `derive_traits` method to get derived traits table."
        raise MissingDataError(msg)

    @property
    def filtered_derived_trait_table(self) -> DerivedTraitTable:
        """The filtered derived trait table."""
        if self._filtered_derived_trait_table is not None:
            return self._filtered_derived_trait_table
        msg = "Please call the `post_filter` method first."
        raise MissingDataError(msg)

    @property
    def diff_results(self) -> dict[str, pd.DataFrame]:
        """The differential analysis results."""
        if self._diff_results is not None:
            return self._diff_results
        if self.groups is None:
            msg = (
                "Grouping information is missing. Please create a new `Experiment` "
                "object with `group_file` first."
            )
            raise MissingDataError(msg)
        else:
            msg = "Please call the `diff_analysis` method first."
            raise MissingDataError(msg)

    # ===== Properties of the data ===== END

    # ===== Workflow methods ===== START
    # === STEP 1 ===
    def preprocess(
        self,
        filter_max_na: float = 1.0,
        impute_method: Literal["zero", "min", "lod", "mean", "median"] = "zero",
    ) -> None:
        """Preprocess the data.

        Calling this method will make the `processed_abundance_table` attribute available.

        Args:
            filter_max_na (float): The maximum ratio of missing values in a sample.
                Range: [0, 1].
                If the ratio of missing values in a sample is greater than this value,
                the sample will be removed.
                Setting to 1.0 means no filtering.
                Setting to 0.0 means only keeping glycans with no missing values.
                Default: 1.0.
            impute_method (str): The method to impute missing values.
                "zero": fill with 0.
                "min": fill with the minimum value of the glycan.
                "lod": fill with the limit of detection of the glycan.
                "mean": fill with the mean value of the glycan.
                "median": fill with the median value of the glycan.
                Default: "zero".
        """
        self._processed_abundance_table = self._do_preprocess(
            abundance_table=self._abundance_table,
            filter_max_na=filter_max_na,
            impute_method=impute_method,
        )
        self._clear_after("preprocess")

    # === STEP 2 ===
    def derive_traits(self, formulas: Optional[list[TraitFormula]] = None) -> None:
        """Calculate derived traits.

        Calling this method will make the `derived_trait_table` attribute available.

        Args:
            formulas: The formulas to calculate the derived traits.
                If not provided, the default formulas will be used.
                Default: None.
        """
        if self._processed_abundance_table is None:
            msg = "Please call the `preprocess` method first."
            raise InvalidOperationOrderError(msg)
        if self._meta_property_table is None:
            msg = "Please call the `extract_meta_properties` method first."
            raise InvalidOperationOrderError(msg)
        self._formulas = self._do_get_formulas(formulas)
        self._derived_trait_table = self._do_derive_traits(
            self._processed_abundance_table, self._meta_property_table, self._formulas
        )
        self._clear_after("derive_traits")

    # === STEP 3 ===
    def post_filter(self, corr_threshold: float = 1.0) -> None:
        """Post-filter the derived traits.

        Calling this method will make the `filtered_derived_trait_table` attribute available.

        Args:
            corr_threshold: The correlation threshold for post filtering.
                If the correlation between two traits is greater than this value,
                one of them will be removed.
                Setting to -1.0 means no correlation filtering.
                Default: 1.0.
        """
        if self._derived_trait_table is None:
            raise InvalidOperationOrderError(
                "Please call the `derive_traits` method first."
            )
        self._filtered_derived_trait_table = self._do_post_filter(
            self._formulas, self._derived_trait_table, corr_threshold
        )
        self._clear_after("post_filter")

    # === STEP 4 ===
    def diff_analysis(self) -> None:
        """Perform differential analysis.

        Calling this method will make the `diff_results` attribute available.
        """
        if self._groups is None:
            raise MissingDataError(
                "Group information is required for differential analysis."
            )
        if self._filtered_derived_trait_table is None:
            raise InvalidOperationOrderError(
                "Please call the `post_filter` method first."
            )
        self._diff_results = self._do_diff_analysis(
            self._filtered_derived_trait_table, self._groups
        )

    # ===== Workflow methods ===== END

    # ===== Business logic methods ===== START
    def _do_preprocess(
        self,
        abundance_table: AbundanceTable,
        filter_max_na: float,
        impute_method: Literal["zero", "min", "lod", "mean", "median"],
    ) -> AbundanceTable:
        """Business logic for the preprocess method."""
        return preprocess(
            data=abundance_table,
            filter_max_na=filter_max_na,
            impute_method=impute_method,
        )

    def _do_get_formulas(self, formulas: list[TraitFormula] | None):
        """Get the formulas to calculate the derived traits."""
        if formulas is None:
            return load_default_formulas(self.mode, self.sia_linkage)
        else:
            if not self.sia_linkage and any(f.sia_linkage for f in formulas):
                raise ValueError(
                    "Could not use SIA linkage formulas with current settings. "
                    "Please create a new Experiment instance with sia_linkage=True."
                )
        return formulas

    def _do_derive_traits(
        self,
        processed_abundance: AbundanceTable,
        meta_properties: MetaPropertyTable,
        formulas: list[TraitFormula],
    ) -> DerivedTraitTable:
        """Business logic for the derive_traits method."""
        return calcu_derived_trait(
            abund_df=processed_abundance,
            meta_prop_df=meta_properties,
            formulas=formulas,
        )

    def _do_post_filter(
        self,
        formulas: list[TraitFormula],
        trait_df: DerivedTraitTable,
        threshold: float,
    ) -> DerivedTraitTable:
        """Business logic for the post_filter method."""
        return post_filter(formulas, trait_df, threshold, method="pearson")

    def _do_diff_analysis(
        self, trait_df: DerivedTraitTable, groups: GroupSeries
    ) -> dict[str, pd.DataFrame]:
        """Business logic for the diff_analysis method."""
        return auto_test(trait_df, groups)

    def _clear_after(self, method_name: str) -> None:
        """Clear the attributes after a certain method."""
        # The order of the keys is the order of the methods in the workflow
        # The values are the attributes generated by the methods
        method_outputs = {
            "preprocess": ["_processed_abundance_table"],
            "derive_traits": ["_formulas", "_derived_trait_table"],
            "post_filter": ["_filtered_derived_trait_table"],
            "diff_analysis": ["_diff_results"],
        }
        method_order = list(method_outputs.keys())

        # Find the index of the current method
        try:
            idx = method_order.index(method_name)
        except ValueError:
            raise ValueError(f"Unknown method name '{method_name}'.")

        # Collect attributes from methods that come after the current one
        attributes_to_clear = []
        for subsequent_method in method_order[idx + 1 :]:
            attributes_to_clear.extend(method_outputs.get(subsequent_method, []))

        # Clear the collected attributes
        for attr in attributes_to_clear:
            setattr(self, attr, None)

    # ===== Business logic methods ===== END

    def run_workflow(
        self,
        *,
        filter_max_na: float = 1.0,
        impute_method: Literal["zero", "min", "lod", "mean", "median"] = "zero",
        formulas: Optional[list[TraitFormula]] = None,
        corr_threshold: float = 1.0,
    ) -> None:
        """Run the entire workflow.

        Call `preprocess`, `derive_traits`, and `post_filter` sequentially.
        If group information is provided, call `diff_analysis` at the end.

        Args:
            filter_max_na: The maximum ratio of missing values in a sample.
                Range: [0, 1].
                If the ratio of missing values in a sample is greater than this value,
                the sample will be removed.
                Setting to 1.0 means no filtering.
                Setting to 0.0 means only keeping glycans with no missing values.
                Default: 1.0.
            impute_method: The method to impute missing values.
                "zero": fill with 0.
                "min": fill with the minimum value of the glycan.
                "lod": fill with the limit of detection of the glycan.
                "mean": fill with the mean value of the glycan.
                "median": fill with the median value of the glycan.
                Default: "zero".
            formulas: The formulas to calculate the derived traits.
                If not provided, the default formulas will be used.
                Default: None.
            corr_threshold: The correlation threshold for post filtering.
                If the correlation between two traits is greater than this value,
                one of them will be removed.
                Setting to -1.0 means no correlation filtering.
                Default: 1.0.
        """
        self.preprocess(filter_max_na, impute_method)
        self.derive_traits(formulas)
        self.post_filter(corr_threshold)
        if self.groups is not None:
            self.diff_analysis()

    def try_formulas(
        self, __expr: Iterable[str] | str, /, squeeze: bool = True
    ) -> pd.DataFrame | pd.Series:
        """Calculate derived traits with the given formulas

        Args:
            __expr: The formulas to calculate the derived traits.
                Could be a list of formula expressions (str)
                or a single formula expression (str).
            squeeze: Whether to squeeze the result if only one formula is provided.
                Default: True.

        Returns:
            If `squeeze=True` and only one formula is provided, return a Series.
            Otherwise, return a DataFrame with derived trait names as columns.
            Index of the DataFrame or Series is the sample names
            (same as the processed abundance table).

        Raises:
            FormulaParseError: If any formula expression is invalid.

        Examples:
            >>> exp = Experiment(...)
            >>> exp.preprocess()
            >>> exp.try_formulas("TC = [type == 'complex'] / [1]")
            pd.Series(...)
            >>> exp.try_formulas(["TC = [type == 'complex'] / [1]", "TS = [nS > 0] / [1]"])
            pd.DataFrame(...)
        """
        if self._processed_abundance_table is None:
            msg = "Please call the `preprocess` method first."
            raise InvalidOperationOrderError(msg)
        if self._meta_property_table is None:
            msg = "Please call the `extract_meta_properties` method first."
            raise InvalidOperationOrderError(msg)

        if isinstance(__expr, str):
            exprs = [__expr]
        else:
            exprs = list(__expr)
        formulas = parse_formulas(exprs)

        trait_table = calcu_derived_trait(
            abund_df=self.processed_abundance_table,
            meta_prop_df=self.meta_property_table,
            formulas=formulas,
        )

        if squeeze and len(formulas) == 1:
            return trait_table.iloc[:, 0]
        else:
            return trait_table
