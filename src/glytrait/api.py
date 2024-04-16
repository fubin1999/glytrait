from typing import Literal, Optional

import pandas as pd
from attrs import define, field

from glytrait.data_export import export_all
from glytrait.data_type import MetaPropertyTable, DerivedTraitTable
from glytrait.formula import TraitFormula, load_formulas
from glytrait.load_data import load_input_data_from_csv, GlyTraitInputData
from glytrait.meta_property import build_meta_property_table
from glytrait.post_filtering import post_filter
from glytrait.preprocessing import preprocess
from glytrait.trait import calcu_derived_trait
from glytrait.stat import t_test, anova


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
        self._formulas: list[TraitFormula] = load_formulas(
            self._config.mode,
            self._config.custom_formula_file,
            self._config.sia_linkage,
        )

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
        return load_input_data_from_csv(
            abundance_file=abundance_file,
            glycan_file=glycan_file,
            group_file=group_file,
            mode=self._config.mode,
        )

    def _preprocess(self, input_data: GlyTraitInputData) -> None:
        preprocess(input_data, self._config.filter_max_na, self._config.impute_method)

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
        derived_trait_table: DerivedTraitTable,
        input_data: GlyTraitInputData,
    ) -> dict[str, pd.DataFrame]:
        if input_data.groups.unique().size == 2:
            return {"t_test.csv": t_test(derived_trait_table, input_data.groups)}
        else:  # groups size > 2
            anova_df, post_hoc_df = anova(derived_trait_table, input_data.groups)
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
