"""The main workflow of GlyTrait.

Functions:
    run_workflow: Run the whole GlyTrait workflow.
"""
from typing import Any, Type, ClassVar, Protocol

import pandas as pd
from attrs import define, field

from glytrait.analysis import differential_analysis, calcu_roc_auc
from glytrait.config import Config
from glytrait.exception import InputError
from glytrait.formula import load_formulas
from glytrait.glycan import load_compositions, load_glycans
from glytrait.io import (
    check_input_file,
    read_structure_file,
    load_default_structures,
    read_group_file,
    write_output,
)

__all__ = ["Workflow"]

from glytrait.meta_property import build_meta_property_table

from glytrait.preprocessing import preprocess_pipeline
from glytrait.trait import calcu_derived_trait, filter_invalid, filter_colinearity


@define
class WorkflowState:
    """The state of workflow.

    Working as a container of all information needed and produced in the workflow.
    This is the way information is passed between workflow steps.
    Information could be some meta-data useful in some steps,
    or some intermediate results produced in some steps.
    Final results will also be stored in this state,
    before being written to the output files in the last step.
    """

    _data: dict[str, Any] = field(factory=dict, repr=False)

    def get(self, key: str) -> Any:
        """Get the value of the key."""
        return self._data.get(key)

    def set(self, key: str, value: Any) -> None:
        """Set the value of the key."""
        self._data[key] = value

    def clear(self) -> None:
        """Clear all data."""
        self._data.clear()


class SupportGetSet(Protocol):
    """A protocol that supports `get` and `set` methods."""

    def get(self, key: str) -> Any:
        ...

    def set(self, key: str, value: Any) -> None:
        ...


@define
class ReadOnlyProxy:
    """A proxy class that prevents modifying the target object.

    This is used to prevent modifying the config and the workflow state
    in the `_execute` method of a workflow step.
    """

    obj: SupportGetSet = field()

    def get(self, key: str) -> Any:
        """Get the value of the key."""
        return self.obj.get(key)

    def set(self, key: str, value: Any, *, verify: bool = False) -> None:
        """Set the value of the key.

        This method is used to prevent modifying the target object accidentally.
        To modify the target object, pass `verify=True` to this method.
        """
        if verify:
            self.obj.set(key, value)
        else:
            raise AttributeError("Cannot set attribute without verification.")

    def __getattr__(self, name: str):
        raise AttributeError(f"Access to method {name} is not allowed")


@define
class WorkflowStep:
    """A workflow step.

    Workflow steps are the building blocks of the workflow.
    The workflow runs by executing these steps one by one.
    The name of a workflow step should ends with "Step",
    for a clear indication of its identity as a workflow step.

    Each step is a class with three methods:
    - '_check_needed': Check if the step is needed.
    - '_execute': Execute the step.
    - 'run': The entry point of the step.

    The `_check_needed` method is used to check if the step is needed.
    Subclasses may not override this method to use the default implementation,
    which always returns True.

    The `_execute` method is used to execute the step.
    Subclasses must override this method to execute the step.
    The method returns a dict of data produced by the step,
    which will be stored in the workflow state in `run`.
    If the step does not produce any data, the method should return None.

    The `run` method is the entry point of the step.
    This method is not supposed to be overridden.

    Extra notes on the `_execute` method:
    1. Don't modify the workflow state manually in the `_execute` method.
    Instead, return the data produced by the step as a dict.
    2. Don't modify the config in the `_execute` method.
    One motivation of doing this is to affect the behavior of the later steps.
    However, this is a bad design.
    Each step should care only about its own behavior.
    Put the logic of deciding whether to run a step in that step's `_check_needed` method.

    In summary, the config and the workflow state should be treated as read-only objects
    in the `_execute` method.
    Also, the `WorkflowStep` base class has a mechanism to prevent modifying them
    in the `_execute` method.


    Examples:
        >>> class SomeStep(WorkflowStep):
        ...     def _execute(self) -> dict[str, Any] | None:
        ...         # Do something
        ...         return {"key": "value"}

        >>> step = SomeStep(...)
        >>> step.run()
    """

    _config: Config = field(repr=False, converter=ReadOnlyProxy)
    _state: WorkflowState = field(repr=False, converter=ReadOnlyProxy)

    def _check_needed(self) -> bool:
        """Check if the step is needed.

        Subclasses may override this method to check if the step is needed.
        """
        return True

    def _execute(self) -> dict[str, Any] | None:
        """Execute the step.

        Subclasses must override this method to execute the step.
        """
        raise NotImplementedError("Subclasses must override this method.")

    def run(self) -> None:
        """The entry point of the step."""
        if self._check_needed():
            data_produced = self._execute()
            if data_produced is not None:
                for key, value in data_produced.items():
                    self._state.set(key, value, verify=True)


@define
class Workflow:
    """The main workflow of GlyTrait.

    Attributes:
        config: Configuration for the workflow.

    Methods:
        run: Run the whole GlyTrait workflow.

    Examples:
        >>> config = Config(...)
        >>> workflow = Workflow(config)
        >>> workflow.run()
    """

    _steps: ClassVar[list[Type[WorkflowStep]]] = []

    config: Config = field(repr=False)
    _state: WorkflowState = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        state_data = {"traits_filtered": False}
        self._state = WorkflowState(state_data)

    @classmethod
    def register_step(cls, step_type: Type[WorkflowStep]) -> Type[WorkflowStep]:
        """Register a workflow step."""
        cls._steps.append(step_type)
        return step_type

    def _reset_state(self) -> None:
        """Reset the workflow state."""
        self._state.clear()

    def run(self) -> None:
        """Run the whole GlyTrait workflow."""
        for step_type in self._steps:
            step = step_type(self.config, self._state)
            step.run()
        self._reset_state()


# ----- Workflow steps -----
@Workflow.register_step
@define
class ReadInputFileStep(WorkflowStep):
    """Read the input file.

    Required meta-data: None

    Produced meta-data:
        input_df: The input file as a pandas.DataFrame.
    """

    def _execute(self) -> dict:
        input_df = pd.read_csv(self._config.get("input_file"))
        return {"input_df": input_df}


@Workflow.register_step
@define
class CheckInputFileStep(WorkflowStep):
    """Check the validity of input file.

    Required meta-data:
        input_df: The input file as a pandas.DataFrame.

    Produced meta-data:
        has_struc_col: Whether the input file has a "Structure" column.

    Errors Raised:
        InputError: If the input file is invalid.
    """

    def _execute(self) -> dict:
        """Execute the step."""
        input_df = self._state.get("input_df")
        check_input_file(input_df)
        # noinspection PyUnreachableCode
        has_struc_col = "Structure" in input_df.columns
        return {"has_struc_col": has_struc_col}


@Workflow.register_step
@define
class CheckHasStructureStep(WorkflowStep):
    """Check if structure information is provided in the structure mode.

    This step is only needed in the "structure" mode (when config.get("mode") == "structure").

    Required meta-data:
        has_struc_col: Whether the input file has a "Structure" column.

    Produced meta-data: None

    Errors Raised:
        InputError: If no structure information is provided.
    """

    def _check_needed(self) -> bool:
        return self._config.get("mode") == "structure"

    def _execute(self) -> None:
        has_struc_col = self._state.get("has_struc_col")
        has_database = self._config.get("database") is not None
        has_structure_file = self._config.get("structure_file") is not None
        if not (has_database or has_structure_file or has_struc_col):
            msg = (
                "Must provide either structure_file or database name when the input file "
                "does not have a 'Structure' column."
            )
            raise InputError(msg)


@Workflow.register_step
@define
class LoadAbundanceTableStep(WorkflowStep):
    """Load the abundance table.

    Required meta-data:
        input_df: The input file as a pandas.DataFrame.

    Produced meta-data:
        abund_df: The abundance table.
    """

    def _execute(self) -> dict[str, Any]:
        """Execute the step."""
        # Get the abundance table
        input_df = self._state.get("input_df")
        input_df = input_df.set_index("Composition")
        if "Structure" in input_df.columns:
            abundance_df = input_df.drop(columns=["Structure"]).T
        else:
            abundance_df = input_df.T
        return {"abund_df": abundance_df}


@Workflow.register_step
@define
class LoadGlycansStep(WorkflowStep):
    """Load the glycan compositions or structures.

    If in the composition mode, the composition strings in the input file
    will be parsed into a list of `Composition` objects.
    If in the structure mode, the structure strings in the input file,
    or the structure strings in a separate structure file,
    or the structures in the build-in database,
    will be parsed into a list of `Structure` objects,
    and the list will be stored in the workflow state.

    Required meta-data:
        input_df: The input file as a pandas.DataFrame.

    Produced meta-data:
        glycans: A list of `Composition` or `NGlycan` objects.

    Errors Raised:
        CompositionParseError: When a composition string in the "composition" mode
            cannot be parsed.
        StructureParseError: When a structure string in the "structure" mode
            cannot be parsed.
    """

    def _execute(self) -> dict[str, Any]:
        input_df = self._state.get("input_df")
        comp_strings = input_df["Composition"].tolist()
        if self._config.get("mode") == "composition":
            glycans = load_compositions(
                comp_strings, sia_linkage=self._config.get("sia_linkage")
            )
        elif self._config.get("mode") == "structure":
            if "Structure" in input_df.columns:
                glycans = load_glycans(comp_strings, input_df["Structure"])
            elif database := self._config.get("database"):
                glycans = load_default_structures(database, comp_strings)
            elif structure_file := self._config.get("structure_file"):
                glycans = read_structure_file(structure_file, comp_strings)
            else:
                raise ValueError("No structure information provided.")
        else:
            raise ValueError(f"Invalid mode: {self._config.get('mode')}")
        return {"glycans": glycans}


@Workflow.register_step
@define
class LoadGroupStep(WorkflowStep):
    """Load the group information.

    Read the group information from the input file,
    and store it in the workflow state.
    This step is optional.

    Required meta-data:
        abund_df: The abundance table.

    Produced meta-data:
        groups: A pandas.Series, with the same index as the abundance table.

    Errors Raised:
        InputError: If the group file has a different set of samples from the input file,
            or there is only one group,
            or there is less than 3 samples in at least one of the groups.
    """

    def _check_needed(self) -> bool:
        return self._config.get("group_file") is not None

    def _execute(self) -> dict[str, Any]:
        groups = read_group_file(self._config.get("group_file"))

        # Check if the group series has the same index as the abundance table (ignore order).
        abund_df = self._state.get("abund_df")
        if not abund_df.index.sort_values().equals(groups.index.sort_values()):
            msg = "The group file must have the same samples as the input file."
            raise InputError(msg)
        # If equal, reindex the group series to the same order as the abundance table.
        groups = groups.reindex(abund_df.index)

        # Check if there are at least 2 groups.
        if groups.nunique() == 1:
            msg = "The group file must have at least 2 groups."
            raise InputError(msg)

        # Check if there are at least 3 samples in each group.
        counts = groups.value_counts()
        counts_less_than_3 = counts[counts < 3]
        if len(counts_less_than_3) > 0:
            msg = (
                f"The following groups have less than 3 samples: "
                f"{', '.join(counts_less_than_3.index)}."
            )
            raise InputError(msg)

        return {"groups": groups}


@Workflow.register_step
@define
class LoadFormulasStep(WorkflowStep):
    """Load the formulas of derived traits.

    Load formulas from the built-in formulas file as well as the user-defined formulas file,
    and store them in the workflow state.
    If "sia_linkage" is False, formulas about sialic acid linkage will be excluded.

    Required meta-data: None

    Produced meta-data:
        formulas: A list of `Formula` objects.

    Errors Raised:
        FormulaError: If a formula string cannot be parsed,
            or the user-provided formula file is in a wrong format.
    """

    def _execute(self) -> dict[str, Any]:
        formulas = load_formulas(
            self._config.get("mode"), self._config.get("formula_file")
        )
        if not self._config.get("sia_linkage"):
            formulas = [f for f in formulas if f.sia_linkage is False]
        return {"formulas": formulas}


@Workflow.register_step
@define
class PreprocessStep(WorkflowStep):
    """Preprocess the abundance table.

    This step will preprocess the abundance table,
    and update "abund_df" and "glycans" in the workflow state.

    Required meta-data:
        abund_df: The abundance table.
        glycans: A list of `Composition` or `NGlycan` objects.

    Produced meta-data:
        abund_df: The preprocessed abundance table.
        glycans: A list of filtered `Composition` or `NGlycan` objects.
    """

    def _execute(self) -> dict[str, Any]:
        abund_df = self._state.get("abund_df")
        glycans = self._state.get("glycans")
        abund_df_processed = preprocess_pipeline(
            abund_df,
            self._config.get("filter_glycan_max_na"),
            self._config.get("impute_method"),
        )
        glycans_filtered = [g for g in glycans if g.name in abund_df_processed.columns]
        return {"abund_df": abund_df_processed, "glycans": glycans_filtered}


@Workflow.register_step
@define
class CalcTraitStep(WorkflowStep):
    """Calculate derived traits.

    This step will calculate the derived traits,
    and store the derived trait table in the workflow state.

    Required meta-data:
        abund_df: The abundance table.
        glycans: A list of `Composition` or `NGlycan` objects.
        formulas: A list of `Formula` objects.

    Produced meta-data:
        meta_property_df: The meta-property table.
        derived_trait_df: The derived trait table.
    """

    def _execute(self) -> dict[str, Any]:
        meta_property_df = build_meta_property_table(
            self._state.get("abund_df").columns.tolist(),
            self._state.get("glycans"),
            self._config.get("mode"),
            self._config.get("sia_linkage"),
        )
        derived_trait_df = calcu_derived_trait(
            self._state.get("abund_df"),
            meta_property_df,
            self._state.get("formulas"),
        )
        return {
            "meta_property_df": meta_property_df,
            "derived_trait_df": derived_trait_df,
        }


@Workflow.register_step
@define
class PostFilteringStep(WorkflowStep):
    """Post-filter the derived traits.

    This step will post-filter the derived traits,
    and update the derived trait table and the formulas in the workflow state.

    Required meta-data:
        abund_df: The abundance table. (For getting the number of samples).
        derived_trait_df: The derived trait table.
        formulas: A list of `Formula` objects.

    Produced meta-data:
        derived_trait_df: The post-filtered derived trait table.
        formulas: The post-filtered formulas.
    """

    def _check_needed(self) -> bool:
        if self._config.get("post_filtering") is False:
            return False
        n_samples = self._state.get("abund_df").shape[0]
        if n_samples < 3:
            msg = "Post-filtering will be skipped when there are less than 3 samples."
            print(msg)
            self._config.set("post_filtering", False, verify=True)
            return False
        return True

    def _execute(self) -> dict[str, Any]:
        formulas = self._state.get("formulas")
        derived_traits = self._state.get("derived_trait_df")
        formulas, derived_traits = filter_invalid(formulas, derived_traits)
        formulas, derived_traits = filter_colinearity(
            formulas,
            derived_traits,
            self._config.get("corr_threshold"),
            self._config.get("corr_method"),
        )
        return {
            "formulas": formulas,
            "derived_trait_df": derived_traits,
            "traits_filtered": True,
        }


@Workflow.register_step
@define
class AnalysisStep(WorkflowStep):
    """Carry out downstream analysis on glycans and derived traits.

    Analysis includes:
    - Univariate analysis
    - ROC analysis

    Required meta-data:
        abund_df: The abundance table.
        derived_trait_df: The derived trait table.
        groups: A series of group labels.
        traits_filtered: Whether the derived traits are filtered.

    Produced meta-data:
        univariate_result: The univariate analysis result.
    """

    def _check_needed(self) -> bool:
        if self._state.get("groups") is None:
            return False
        if self._state.get("traits_filtered") is False:
            print("Downstream analysis will be skipped when post-filtering is off.")
            return False
        return True

    def _execute(self) -> dict[str, Any]:
        abund_df = self._state.get("abund_df")
        derived_trait_df = self._state.get("derived_trait_df")
        combined_df = pd.concat([abund_df, derived_trait_df], axis=1)
        groups = self._state.get("groups")
        univariate_result = differential_analysis(combined_df, groups)
        if groups.nunique() == 2:
            roc_result = calcu_roc_auc(combined_df, self._state.get("groups"))
        else:
            roc_result = None
        return {"univariate_result": univariate_result, "roc_result": roc_result}


@Workflow.register_step
@define
class WriteOutputStep(WorkflowStep):
    """Write the output files.

    Required meta-data:
        derived_trait_df: The derived trait table.
        abund_df: The abundance table.
        meta_property_df: The meta-property table.
        formulas: A list of `Formula` objects.
        groups: A series of group labels. (optional)
        univariate_result: The univariate analysis result. (optional)
        roc_result: The ROC analysis result. (optional)

    Produced files: None
    """

    def _execute(self) -> None:
        write_output(
            self._config.obj,
            self._state.get("derived_trait_df"),
            self._state.get("abund_df"),
            self._state.get("meta_property_df"),
            self._state.get("formulas"),
            self._state.get("groups"),
            self._state.get("univariate_result"),
            self._state.get("roc_result"),
        )


# ---- END of workflow steps ----
