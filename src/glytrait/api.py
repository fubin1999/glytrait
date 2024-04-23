from functools import wraps
from typing import cast, Literal, Optional, ClassVar, Any

import pandas as pd
from attrs import define, field

from glytrait.data_export import export_all
from glytrait.data_type import (
    MetaPropertyTable,
    DerivedTraitTable,
    GroupSeries,
    AbundanceTable,
)
from glytrait.exception import GlyTraitError
from glytrait.formula import (
    TraitFormula,
    load_formulas_from_file,
    load_default_formulas,
    parse_formulas,
)
from glytrait.load_data import GlyTraitInputData, load_data
from glytrait.meta_property import build_meta_property_table
from glytrait.post_filtering import post_filter
from glytrait.preprocessing import preprocess
from glytrait.trait import calcu_derived_trait
from glytrait.stat import t_test, anova
from glytrait.glycan import Structure, Composition


class GlyTrait:
    """GlyTrait API.

    Args:
        mode (str): "structure" or "composition". Default: "structure".
        filter_max_na (float): The maximum ratio of missing values in a sample.
            If the ratio of missing values in a sample is greater than this value,
            the sample will be removed.
            Setting to 1.0 means no filtering.
            Setting to 0.0 means only keeping glycans with no missing values.
            Default: 1.0.
        impute_method (str): The method to impute missing values.
            "zero": fill with 0.
            "min": fill with the minimum value of the column.
            "lod": fill with the limit of detection of the column.
            "mean": fill with the mean value of the column.
            "median": fill with the median value of the column.
            Default: "zero".
        post_filtering (bool): Whether to perform post filtering.
            If True, the invalid traits and the highly correlated traits will be removed.
            Default: True.
        correlation_threshold (float): The correlation threshold for post filtering.
            If the correlation between two traits is greater than this value,
            one of them will be removed.
            Setting to -1.0 means no correlation filtering.
            Default: 1.0.
        sia_linkage (bool): Whether to consider the linkage of sialic acid.
            If True, the linkage of sialic acid will be considered in the calculation of
            meta properties. Default: False.
        custom_formula_file (str): Path to the custom formula file.

    Examples:
        >>> from glytrait import GlyTrait
        >>> glytrait = GlyTrait()
        >>> glytrait.run(
        ...     output_dir="output",
        ...     abundance_file="glycan_abundance.csv",
        ...     glycan_file="glycan_structure.csv",
        ...     group_file="group.csv",
        ... )
    """

    def __init__(self, **kwargs):
        self._config = _Config(**kwargs)
        self._formulas = self._init_formulas()

    def _init_formulas(self) -> list[TraitFormula]:
        if self._config.custom_formula_file is not None:
            return load_formulas_from_file(
                self._config.custom_formula_file, self._config.sia_linkage
            )
        else:
            return load_default_formulas(self._config.mode, self._config.sia_linkage)

    def run(
        self,
        output_dir: str,
        abundance_file: str,
        glycan_file: str,
        group_file: Optional[str] = None,
    ):
        """Run GlyTrait.

        Three files are required: abundance_file, glycan_file, and group_file.

        The abundance file is a csv file with the first column as the sample names, and the
        remaining columns as the abundance values for different glycans.
        The first column should be named "Sample".

        Two choices for the glycan file:
        1. A csv file with the first column "GlycanID" as the glycan names,
        and the second column "Structure" or "Composition" as the glycan
        structure or composition.
        2. A directory with each file named as the glycan names, and the file
        content as the glycan structure.
        The second choice is only available when mode is "structure".
        At this moment, only "glycoct" format is supported.
        So all files in this directory should have the extension ".glycoct".

        The group file is a csv file with the first column "Sample" as the sample names,
        and the second column "Group" as the group names.

        Args:
            output_dir (str): Path to the output directory.
            abundance_file (str): Path to the abundance file.
            glycan_file (str): Path to the glycan file.
            group_file (str): Path to the group file. Optional.
        """
        input_data = self._load_input_data(abundance_file, glycan_file, group_file)
        self._preprocess(input_data)
        meta_property_table = self._calcu_meta_property(input_data)
        derived_trait_table = self._calcu_derived_trait(input_data, meta_property_table)
        if self._config.post_filtering:
            filtered_derived_trait_table = self._post_filtering(derived_trait_table)
        else:
            filtered_derived_trait_table = None
        if group_file is not None and self._config.post_filtering:
            diff_result = self._diff_analysis(filtered_derived_trait_table, input_data)
        else:
            diff_result = {}
        self._export_data(
            output_dir,
            input_data,
            meta_property_table,
            derived_trait_table,
            filtered_derived_trait_table,
            diff_result,
        )

    def _load_input_data(
        self, abundance_file: str, glycan_file: str, group_file: Optional[str] = None
    ) -> GlyTraitInputData:
        abundance_df = pd.read_csv(abundance_file)
        glycan_df = pd.read_csv(glycan_file)
        if group_file is not None:
            group_df = pd.read_csv(group_file)
        else:
            group_df = None
        return load_data(abundance_df, glycan_df, group_df, mode=self._config.mode)

    def _preprocess(self, input_data: GlyTraitInputData) -> None:
        processed_abund_df = preprocess(
            data=input_data.abundance_table,
            filter_max_na=self._config.filter_max_na,
            impute_method=self._config.impute_method,
        )
        input_data._abundance_table = processed_abund_df
        for glycan in set(input_data.glycans).difference(processed_abund_df.columns):
            del input_data.glycans[glycan]

    def _calcu_meta_property(self, input_data: GlyTraitInputData) -> MetaPropertyTable:
        return build_meta_property_table(
            input_data.glycans, self._config.mode, self._config.sia_linkage
        )

    def _calcu_derived_trait(
        self, input_data: GlyTraitInputData, meta_property_table: MetaPropertyTable
    ) -> DerivedTraitTable:
        return calcu_derived_trait(
            input_data.abundance_table, meta_property_table, self._formulas
        )

    def _post_filtering(
        self, derived_trait_table: DerivedTraitTable
    ) -> DerivedTraitTable:
        return post_filter(
            formulas=self._formulas,
            trait_df=derived_trait_table,
            threshold=self._config.correlation_threshold,
            method="pearson",
        )

    def _diff_analysis(
        self,
        derived_trait_table: DerivedTraitTable | None,
        input_data: GlyTraitInputData,
    ) -> dict[str, pd.DataFrame]:
        groups = cast(GroupSeries, input_data.groups)
        trait_table = cast(DerivedTraitTable, derived_trait_table)
        if groups.unique().size == 2:
            return {"t_test.csv": t_test(trait_table, groups)}
        else:  # groups size > 2
            anova_df, post_hoc_df = anova(trait_table, groups)
            return {"anova.csv": anova_df, "post_hoc.csv": post_hoc_df}

    def _export_data(
        self,
        output_dir: str,
        input_data: GlyTraitInputData,
        meta_property_table: MetaPropertyTable,
        derived_trait_table: DerivedTraitTable,
        filtered_derived_trait_table: Optional[DerivedTraitTable] = None,
        diff_results: Optional[dict[str, pd.DataFrame]] = None,
    ) -> None:
        data_to_export = [
            ("meta_properties.csv", meta_property_table),
            ("derived_traits.csv", derived_trait_table),
            ("glycan_abundance_processed.csv", input_data.abundance_table),
        ]
        if filtered_derived_trait_table is not None:
            data_to_export.append(
                ("derived_traits_filtered.csv", filtered_derived_trait_table)
            )
        if diff_results is not None:
            data_to_export.extend(diff_results.items())
        export_all(data_to_export, output_dir)


@define
class _Config:
    """GlyTrait configuration.

    This class encapsulates the configuration validation logic.
    """

    mode: Literal["structure", "composition"] = field(default="structure")
    filter_max_na: float = field(default=1.0)
    impute_method: Literal["zero", "min", "lod", "mean", "median"] = field(
        default="zero"
    )
    post_filtering: bool = field(default=True)
    correlation_threshold: float = field(default=1.0)
    sia_linkage: bool = field(default=False)
    custom_formula_file: Optional[str] = field(default=None)

    @mode.validator
    def _validate_mode(self, attribute, value):  # type: ignore
        if not isinstance(value, str):
            raise ValueError("mode must be a string.")
        if value not in {"structure", "composition"}:
            raise ValueError("mode must be one of: structure, composition.")

    @filter_max_na.validator
    def _validate_filter_max_na(self, attribute, value):  # type: ignore
        if not isinstance(value, (float, int)):
            raise ValueError("filter_max_na must be a float.")
        if not 0 <= value <= 1:
            raise ValueError("filter_max_na must be between 0 and 1.")

    @impute_method.validator
    def _validate_impute_method(self, attribute, value):  # type: ignore
        if not isinstance(value, str):
            raise ValueError("impute_method must be a string.")
        if value not in {"zero", "min", "lod", "mean", "median"}:
            raise ValueError(
                "impute_method must be one of: zero, min, lod, mean, median."
            )

    @post_filtering.validator
    def _validate_post_filtering(self, attribute, value):  # type: ignore
        if not isinstance(value, bool):
            raise ValueError("post_filtering must be a boolean.")

    @correlation_threshold.validator
    def _validate_correlation_threshold(self, attribute, value):  # type: ignore
        if not isinstance(value, (float, int)):
            raise ValueError("correlation_threshold must be a float.")
        if not (0 <= value <= 1 or value == -1):
            raise ValueError("correlation_threshold must be between 0 and 1, or -1.")

    @sia_linkage.validator
    def _validate_sia_linkage(self, attribute, value):  # type: ignore
        if not isinstance(value, bool):
            raise ValueError("sia_linkage must be a boolean.")

    @custom_formula_file.validator
    def _validate_custom_formula_file(self, attribute, value):  # type: ignore
        if value is not None:
            if not isinstance(value, str):
                raise ValueError("custom_formula_file must be a string.")


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

    Create an instance of this class with an abundance table and a structure series
    to perform all the steps in the GlyTrait workflow.

    Firstly, call the `preprocessed` method to filter glycans, impute missing values,
    and normalize the abundance.
    Calling this method makes the `processed_abundance_table` attribute available.

    Secondly, call the `extract_meta_properties` method.
    This method generates a table of meta-properties for each glycan,
    stored as the `meta_property_table` attribute.

    Thirdly, call the `derive_traits` method to calculate all the derived traits.
    The result table is stored as the `derived_trait_table` attribute.

    Fourthly, call the `post_filter` method
    to remove invalid traits and highly correlated traits.
    The result table is stored as the `filtered_derived_trait_table` attribute.

    Finally, and optionally, call the `diff_analysis` method with a group series
    to perform differential analysis.
    The result is stored as the `diff_results` attribute.

    Besides, there is a `derive_one_trait` method to try a single trait formula.
    It returns a Series instead of storing the result.
    This method is useful to test custom trait formulas.

    Methods:
        preprocess: Preprocess the data.
        extract_meta_properties: Extract meta properties.
        derive_traits: Calculate derived traits.
        post_filter: Post-filter the derived traits.
        diff_analysis: Perform differential analysis.
        derive_one_trait: Try a single trait formula.
        reset: Reset the workflow to the beginning (before `preprocess`).
        get_data: Get the data by name.

    Attributes:
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
        >>> from glytrait.load_data import load_data
        # This is the safe way to load data.
        # It ensures all data are in correct formats.
        >>> input_data = load_data(
        ...    abundance_df=pd.read_csv("glycan_abundance.csv"),
        ...    glycan_df=pd.read_csv("glycan_structure.csv"),
        ...    group_df=pd.read_csv("group.csv"),
        ...    mode="structure",
        )
        # `input_data` is a GlyTraitInputData instance.
        >>> experiment = Experiment(input_data)
        >>> experiment.preprocess(filter_max_na=0.5, impute_method="min")
        >>> experiment.extract_meta_properties()
        >>> experiment.derive_traits()  # with default formulas
        >>> experiment.post_filter(corr_threshold=0.9)
        >>> experiment.diff_analysis()
    """

    _all_steps = [
        "preprocess",
        "extract_meta_properties",
        "derive_traits",
        "post_filter",
        "diff_analysis",
    ]
    _data_dict = {
        "preprocess": ["processed_abundance_table"],
        "extract_meta_properties": ["meta_property_table"],
        "derive_traits": ["derived_trait_table", "formulas"],
        "post_filter": ["filtered_derived_trait_table"],
        "diff_analysis": ["diff_results"],
    }

    input_data: GlyTraitInputData
    mode: Literal["structure", "composition"] = field(default="structure", kw_only=True)
    sia_linkage: bool = field(default=False, kw_only=True)

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
        processed = preprocess(
            data=self.input_data.abundance_table,
            filter_max_na=filter_max_na,
            impute_method=impute_method,
        )
        return {"processed_abundance_table": AbundanceTable(processed)}  # type: ignore

    @_step
    def extract_meta_properties(self) -> None:
        """Extract meta-properties.

        Calling this method will make the `meta_property_table` attribute available.
        """
        glycans: list[str] = self.processed_abundance_table.columns.tolist()
        glycan_dict = cast(GlycanDict, {g: self.input_data.glycans[g] for g in glycans})
        mp_table = build_meta_property_table(glycan_dict, self.mode, self.sia_linkage)
        return {"meta_property_table": mp_table}  # type: ignore

    @_step
    def derive_traits(self, exprs: Optional[list[str]] = None) -> None:
        """Calculate derived traits.

        Calling this method will make the `derived_trait_table` attribute available.

        Args:
            exprs: If provided, these formula expressions will be used to calculate
                the derived traits.
                Otherwise, the default formulas will be used.
        """
        if exprs is None:
            formulas = load_default_formulas(self.mode, self.sia_linkage)
        else:
            formulas = parse_formulas(exprs)
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
        if groups.unique().size == 2:
            result = {"t_test.csv": t_test(trait_table, groups)}
        else:  # groups size > 2
            anova_df, post_hoc_df = anova(trait_table, groups)
            result = {"anova.csv": anova_df, "post_hoc.csv": post_hoc_df}
        return {"diff_results": result}  # type: ignore
