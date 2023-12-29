from typing import Literal, Optional

from attrs import define, field

from glytrait.formula import TraitFormula, load_formulas
from glytrait.load_data import load_input_data, GlyTraitInputData
from glytrait.preprocessing import preprocess
from glytrait.meta_property import build_meta_property_table
from glytrait.trait import calcu_derived_trait
from glytrait.post_filtering import filter_invalid, filter_colinearity
from glytrait.data_export import export_all


@define
class GlyTrait:
    """GlyTrait API."""

    mode: Literal["structure", "composition"] = field()
    output_dir: str = field()
    filter_max_na: float = field(default=0.0)
    impute_method: Literal["zero", "min", "lod", "mean", "median"] = field(
        default="zero"
    )
    post_filtering: bool = field(default=True)
    correlation_threshold: bool = field(default=0.9)
    sia_linkage: bool = field(default=False)
    custom_formula_file: Optional[str] = field(default=None)

    _formulas: list[TraitFormula] = field(init=False)

    @mode.validator
    def _validate_mode(self, attribute, value):  # type: ignore
        if not isinstance(value, str):
            raise ValueError("mode must be a string.")
        if value not in {"structure", "composition"}:
            raise ValueError("mode must be one of: structure, composition.")

    @output_dir.validator
    def _validate_output_dir(self, attribute, value):
        if not isinstance(value, str):
            raise ValueError("output_dir must be a string.")

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
            raise ValueError("impute_method must be one of: zero, min, lod, mean, median.")

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

    def __attrs_post_init__(self):
        self._formulas = load_formulas(self.mode, self.custom_formula_file, self.sia_linkage)

    def run(
        self,
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
            abundance_file (str): Path to the abundance file.
            glycan_file (str): Path to the glycan file.
            group_file (str): Path to the group file. Optional.
        """
        input_data: GlyTraitInputData = load_input_data(
            abundance_file=abundance_file,
            glycan_file=glycan_file,
            group_file=group_file,
            mode=self.mode,
        )
        formulas = self._formulas
        preprocess(input_data, self.filter_max_na, self.impute_method)
        meta_property_table = build_meta_property_table(
            input_data.glycans, self.mode, self.sia_linkage
        )
        derived_trait_table = calcu_derived_trait(
            input_data.abundance_table, meta_property_table, formulas
        )
        if self.post_filtering:
            formulas, derived_trait_table = filter_invalid(
                formulas, derived_trait_table
            )
            formulas, derived_trait_table = filter_colinearity(
                formulas,
                derived_trait_table,
                self.correlation_threshold,
                method="pearson",
            )
        data_to_export = [
            ("formulas.txt", formulas),
            ("meta_properties.csv", meta_property_table),
            ("derived_traits.csv", derived_trait_table),
            ("glycan_abundance_processed.csv", input_data.abundance_table),
        ]
        export_all(data_to_export, self.output_dir)
