from collections import defaultdict
from collections.abc import Iterable
from functools import wraps
from typing import cast, Literal, Optional, ClassVar, Any

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
from glytrait.data_input import GlyTraitInputData, load_data
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


@define
class _Workflow:
    """Manage the method calling order for the `Experiment` class.

    To use this class, create a subclass of it,
    set `_all_steps` and `_data_dict` class variables,
    and decorate each method with the `_step` decorator.

    The `_all_steps` attribute is a list of method names in the correct order.
    The `_data_dict` attribute is a dictionary where the keys are method names,
    and the values are the data keys each method will generate.

    All methods being a workflow step should return a dictionary,
    where the keys are the data keys, and the values are the data.
    Remember to decorate each step method with the `_step` decorator.
    This ensures the methods are converted to the correct form.

    Following this pattern will ensure:

    1. If a method is called out of order, an `InvalidOperationOrderError` will be raised.
    2. If a method is called twice, an `InvalidOperationOrderError` will be raised,
         unless `force=True` is passed.
    3. The data after a method will be reset if `force=True` is passed.
    4. The data generated by each method will be stored in the `_data` attribute.

    Also, a `reset` method is provided to reset the workflow,
    and the current step name can be accessed via the `current_step` property.
    The current step is "__START__" at the beginning,
    or after calling the `reset` method.

    Finally, the data generated by step methods can be accessed via the `get_data` method.
    """

    # `_all_steps` contains all the method names in the correct order.
    _all_steps: ClassVar[list[str]] = list()

    # `_data_dict` is the data keys each method will generate.
    _data_dict: ClassVar[dict[str, list[str]]] = dict()

    _data: dict[str, Any] = field(init=False, factory=dict)
    _current_step: str = field(init=False, default="__START__")

    def __attrs_post_init__(self) -> None:
        # Check if all methods in `_all_steps` are implemented.
        for step in self._all_steps:
            if not hasattr(self, step):
                raise NotImplementedError(f"Method {step} is not implemented.")

        # Check if the keys of `_data_dict` are the same as `_all_steps`.
        if set(self._all_steps) != set(self._data_dict.keys()):
            raise ValueError("`_all_steps` and `_data_dict` should have the same keys.")

    @property
    def current_step(self) -> str:
        """The current step name."""
        return self._current_step

    def reset(self) -> None:
        """Reset the workflow."""
        self._data.clear()
        self._current_step = "__START__"

    def _reset_data_after(self, step_name: str) -> None:
        """Reset all data after the given step."""
        for key in self._all_steps[self._all_steps.index(step_name) + 1 :]:
            for data_key in self._data_dict[key]:
                self._data.pop(data_key, None)

    def get_data(self, data_name: str) -> Any:
        """Get the data by name."""
        try:
            return self._data[data_name]
        except KeyError:
            if method_name := self._get_method_name(data_name):
                raise MissingDataError(
                    f"Data '{data_name}' is not available. "
                    f"Call the corresponding method ('{method_name}') to generate it."
                )
            else:
                raise KeyError(
                    f"Data '{data_name}' is not generated by any method. "
                    f"You may have misspelled the data name."
                )

    def _get_method_name(self, data_name: str) -> str | None:
        """Get the method name that generates the data."""
        for method, data_keys in self._data_dict.items():
            if data_name in data_keys:
                return method
        return None


def _step(func):
    """This function should be used along with `Workflow` class.

    It adds the following features to the decorated method:

    1. Check if the method is called in the correct order.
    2. Reset all data after the method if `force=True` is passed.
    3. Update the current step name.
    4. Store the result of the method in the `_data` attribute.

    See Also:
        `Workflow` class.
    """

    @wraps(func)
    def wrapper(self, *args, force=False, **kwargs):
        # This is the index of the last called step method
        if self._current_step == "__START__":
            last_pos = -1
        else:
            last_pos = self._all_steps.index(self._current_step)

        # This is the index of this step method
        this_pos = self._all_steps.index(func.__name__)

        # The current step is exactly the next step of the last called step
        if this_pos == last_pos + 1:
            self._current_step = func.__name__
            self._data.update(func(self, *args, **kwargs))

        # The current step missed some prerequisites
        elif this_pos > last_pos + 1:
            missing_steps = self._all_steps[last_pos + 1 : this_pos]
            raise InvalidOperationOrderError(
                f"Missing steps: {', '.join(missing_steps)}. "
                f"Call them sequentially before calling {func.__name__}."
            )

        # The current step is already called
        elif this_pos <= last_pos:
            if force:
                self._current_step = func.__name__
                self._reset_data_after(func.__name__)
                self._data.update(func(self, *args, **kwargs))
            else:
                raise InvalidOperationOrderError(
                    f"Step {func.__name__} is already called. "
                    f"Use `force=True` to call it again. "
                    f"Note that this will reset all data after this step."
                )

    return wrapper


GlycanDict = dict[str, Structure] | dict[str, Composition]


@define
class Experiment(_Workflow):
    """GlyTrait experiment.

    Create an instance of this class to perform all the steps in the GlyTrait workflow.

    Firstly, call the `preprocessed` method to filter glycans, impute missing values,
    and normalize the abundance.
    Calling this method makes the `processed_abundance_table` and
    `meta_property_table` attributes available.
    The latter is a table of structural or compositional properties,
    which is used to calculate derived traits.

    Secondly, call the `derive_traits` method to calculate all the derived traits.
    The result table is stored as the `derived_trait_table` attribute.

    Thirdly, call the `post_filter` method
    to remove invalid traits and highly correlated traits.
    The result table is stored as the `filtered_derived_trait_table` attribute.

    Finally, and optionally, call the `diff_analysis` method with a group series
    to perform differential analysis.
    The result is stored as the `diff_results` attribute.

    Alternatively, you can call the `run_workflow` method to run the entire workflow.

    If you want to try new formula expressions, use the `try_formulas` method.
    Instead of saving the result as an attribute,
    this method returns the result DataFrame or Series directly.

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
        input_data: The input data (a `GlyTraitInputData` instance).
        abundance_table: The original abundance table.
        glycans: The glycans.
        groups: The groups.
        processed_abundance_table: The processed abundance table.
        meta_property_table: The meta property table.
        derived_trait_table: The derived trait table.
        filtered_derived_trait_table: The filtered derived trait table.
        diff_results: The differential analysis results.

    Examples:
        >>> from glytrait.api import Experiment
        >>> experiment = Experiment(
        ...    abundance_file="glycan_abundance.csv",
        ...    glycan_file="glycan_structure.csv",
        ...    group_file="group.csv",
        ...    mode="structure",
        ... )
        >>> experiment.preprocess(filter_max_na=0.5, impute_method="min")
        >>> experiment.derive_traits()  # with default formulas
        >>> experiment.post_filter(corr_threshold=0.9)
        >>> experiment.diff_analysis()
    """

    _all_steps = [
        "preprocess",
        "derive_traits",
        "post_filter",
        "diff_analysis",
    ]
    _data_dict = {
        "preprocess": ["processed_abundance_table", "meta_property_table"],
        "derive_traits": ["derived_trait_table", "formulas"],
        "post_filter": ["filtered_derived_trait_table"],
        "diff_analysis": ["diff_results"],
    }

    abundance_file: Optional[str] = field(default=None, kw_only=True, repr=False)
    glycan_file: Optional[str] = field(default=None, kw_only=True, repr=False)
    group_file: Optional[str] = field(default=None, kw_only=True, repr=False)
    input_data: GlyTraitInputData = field(default=None, kw_only=True)
    mode: Literal["structure", "composition"] = field(default="structure", kw_only=True)
    sia_linkage: bool = field(default=False, kw_only=True)

    def __attrs_post_init__(self):
        if self.abundance_file is None and self.glycan_file is None:
            if self.input_data is None:
                msg = "Either `input_data` or `abundance_file` and `glycan_file` are required."
                raise ValueError(msg)
        if self.abundance_file and self.glycan_file is None:
            raise ValueError(
                "`glycan_file` is required if `abundance_file` is provided."
            )
        if self.abundance_file is None and self.glycan_file:
            raise ValueError(
                "`abundance_file` is required if `glycan_file` is provided."
            )
        if self.input_data and (self.abundance_file or self.glycan_file):
            raise ValueError("Do not provide both `input_data` and files.")

        abundance_type = defaultdict(lambda: "float64")
        abundance_type.update({"Sample": "O"})
        if self.input_data is None:
            self.input_data = load_data(
                abundance_df=pd.read_csv(self.abundance_file, dtype=abundance_type),
                glycan_df=pd.read_csv(self.glycan_file),
                group_df=pd.read_csv(self.group_file) if self.group_file else None,
                mode=self.mode,
            )

    @property
    def abundance_table(self) -> pd.DataFrame:
        """The original abundance table."""
        return self.input_data.abundance_table

    @abundance_table.setter
    def abundance_table(self, value: pd.DataFrame) -> None:
        self.input_data.abundance_table = value  # type: ignore
        self.reset()

    @property
    def glycans(self) -> GlycanDict:
        """The glycans."""
        return self.input_data.glycans

    @glycans.setter
    def glycans(self, value: GlycanDict) -> None:
        self.input_data.glycans = value
        self.reset()

    @property
    def groups(self) -> GroupSeries | None:
        """The groups."""
        return self.input_data.groups

    @groups.setter
    def groups(self, value: pd.Series | None) -> None:
        self.input_data.groups = value  # type: ignore
        self.reset()

    @property
    def processed_abundance_table(self) -> AbundanceTable:
        """The processed abundance table."""
        return self.get_data("processed_abundance_table")

    @property
    def meta_property_table(self) -> MetaPropertyTable:
        """The meta property table."""
        return self.get_data("meta_property_table")

    @property
    def derived_trait_table(self) -> DerivedTraitTable:
        """The derived trait table."""
        return self.get_data("derived_trait_table")

    @property
    def filtered_derived_trait_table(self) -> DerivedTraitTable:
        """The filtered derived trait table."""
        return self.get_data("filtered_derived_trait_table")

    @property
    def diff_results(self) -> dict[str, pd.DataFrame]:
        """The differential analysis results."""
        return self.get_data("diff_results")

    @_step
    def preprocess(
        self,
        filter_max_na: float = 1.0,
        impute_method: Literal["zero", "min", "lod", "mean", "median"] = "zero",
    ) -> None:
        """Preprocess the data.

        Calling this method will make the `processed_abundance_table` and
        `meta_property_table` attributes available.

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
        processed = preprocess(
            data=self.input_data.abundance_table,
            filter_max_na=filter_max_na,
            impute_method=impute_method,
        )
        mp_table = self._extract_meta_properties(processed)
        return {  # type: ignore
            "processed_abundance_table": AbundanceTable(processed),
            "meta_property_table": mp_table,
        }

    def _extract_meta_properties(
        self, processed_abund_df: AbundanceTable
    ) -> MetaPropertyTable:
        """Extract meta-properties.

        Calling this method will make the `meta_property_table` attribute available.
        """
        glycans: list[str] = processed_abund_df.columns.tolist()
        glycan_dict = cast(GlycanDict, {g: self.input_data.glycans[g] for g in glycans})
        mp_table = build_meta_property_table(glycan_dict, self.mode, self.sia_linkage)
        return mp_table

    @_step
    def derive_traits(self, formulas: Optional[list[TraitFormula]] = None) -> None:
        """Calculate derived traits.

        Calling this method will make the `derived_trait_table` attribute available.

        Args:
            formulas: The formulas to calculate the derived traits.
                If not provided, the default formulas will be used.
                Default: None.
        """
        if formulas is None:
            formulas = load_default_formulas(self.mode, self.sia_linkage)
        else:
            if not self.sia_linkage and any(f.sia_linkage for f in formulas):
                raise ValueError(
                    "Could not use SIA linkage formulas with current settings. "
                    "Please create a new Experiment instance with sia_linkage=True."
                )
        trait_table = calcu_derived_trait(
            abund_df=self.processed_abundance_table,
            meta_prop_df=self.meta_property_table,
            formulas=formulas,
        )
        return {"derived_trait_table": trait_table, "formulas": formulas}  # type: ignore

    @_step
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
        filtered_table = post_filter(
            formulas=self.get_data("formulas"),
            trait_df=self.derived_trait_table,
            threshold=corr_threshold,
            method="pearson",
        )
        return {"filtered_derived_trait_table": filtered_table}  # type: ignore

    @_step
    def diff_analysis(self) -> None:
        """Perform differential analysis.

        Calling this method will make the `diff_results` attribute available.
        """
        groups = self.input_data.groups
        trait_table = self.filtered_derived_trait_table
        if groups is None:
            raise MissingDataError(
                "Group information is required for differential analysis."
            )
        result = auto_test(trait_table, groups)
        return {"diff_results": result}  # type: ignore

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
        if self.input_data.groups is not None:
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
            >>> exp = Experiment(input_data)
            >>> exp.preprocess()
            >>> exp.try_formulas("TC = [type == 'complex'] / [1]")
            pd.Series(...)
            >>> exp.try_formulas(["TC = [type == 'complex'] / [1]", "TS = [nS > 0] / [1]"])
            pd.DataFrame(...)
        """

        if self._current_step == "__START__":
            raise InvalidOperationOrderError("Call `preprocess` first.")

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
