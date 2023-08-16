"""The main page of the GlyTrait Streamlit app."""
from typing import ClassVar, IO

import pandas as pd
import streamlit as st
from attrs import define, field

from glytrait.exception import GlyTraitError
from glytrait.formula import load_formulas
from glytrait.glycan import load_compositions, load_glycans
from glytrait.io import load_default_structures, read_structure_file
from glytrait.meta_property import build_meta_property_table
from glytrait.preprocessing import preprocess_pipeline
from glytrait.trait import calcu_derived_trait, filter_invalid, filter_colinearity


@st.cache_data(ttl=3600)
def load_data(data: IO) -> pd.DataFrame:
    """Load data from a csv file."""
    return pd.read_csv(data)


__BIND_TO_WEIGHT__ = "bind_to_weight"

if "finished" not in st.session_state:
    st.session_state["finished"] = False


@define
class StreamlitStep:
    """A step in the Streamlit app."""

    to_be_initialzed: ClassVar[dict] = {}

    def check_passed(self) -> bool:
        """Check if the step has passed."""
        return True

    def check_needed(self) -> bool:
        """Check if the step is needed."""
        return True

    def run(self) -> None:
        """Run the step."""
        raise NotImplementedError

    def init_state(self):
        """Initialize the state of the step.

        If a value equals to `__BIND_TO_WEIGHT`,
        the key is controlled by a weight,
        and should not be set mannually.
        """
        for key, value in self.to_be_initialzed.items():
            if key not in st.session_state and value != __BIND_TO_WEIGHT__:
                st.session_state[key] = value

    def clear_state(self):
        """Clear the state of the step."""
        for key in self.to_be_initialzed.keys():
            if key in st.session_state:
                del st.session_state[key]


@define
class MainPage:
    """The main page of the Streamlit app."""

    _steps: list[StreamlitStep] = field(factory=list, repr=False)

    def run(self) -> None:
        """Run the main page."""
        it = iter(self._steps)
        while True:
            try:
                step = next(it)
            except StopIteration:
                break
            if step.check_needed():
                step.init_state()
                step.run()
                if not step.check_passed():
                    break
        while True:
            try:
                step = next(it)
            except StopIteration:
                break
            step.clear_state()


@define
class WelcomeStep(StreamlitStep):
    """Step1: show the welcome message."""

    def run(self) -> None:
        st.title("GlyTrait")
        st.markdown(
            "A tool for N-glycan derived trait calculation, and downstream analysis."
        )
        st.markdown(
            "For more information, please visit the [GitHub repository]"
            "(https://github.com/fubin1999/glytrait)."
        )


@define
class ModeStep(StreamlitStep):
    """Step2: ask the user to choose a mode."""

    to_be_initialzed = {"mode": __BIND_TO_WEIGHT__}

    def run(self) -> None:
        st.header("Mode")
        st.markdown(
            """GlyTrait supports two modes: **structure** mode and **composition** mode.
If structure information is available, the structure mode is recommended.
        """
        )
        st.radio("Select a mode: ", ("structure", "composition"), key="mode")


@define
class UploadFileStep(StreamlitStep):
    """Step3: upload the input file."""

    to_be_initialzed = {"input_file": __BIND_TO_WEIGHT__, "input_df": None}

    def run(self) -> None:
        st.header("Upload file")
        self._show_instruction()
        self._inquire_upload_file()

    def check_passed(self) -> bool:
        return st.session_state.input_file is not None

    def _show_instruction(self) -> None:
        if st.session_state.get("mode") == "structure":
            self._show_structure_mode_instruction()
        else:
            self._show_composition_mode_instruction()

    def _show_structure_mode_instruction(self) -> None:
        st.subheader("File format instruction")
        st.markdown(
            """The upload file should be a csv file with the following columns:
- **The first column**: "Composition", the glycan compositions. e.g. H5N4F1S1
- **The second column**: "Structure", the glycan structures in condensed glycoCT format.
- **From the third column**: glycan abundance in each sample, with column names as sample names.

*The "Structure" column is not needed when a structure file is provided, 
or a build-in database is used.
Check out the github page for more information.*"""
        )
        st.markdown("An example input file:")
        example_input_df = pd.DataFrame(
            {
                "Composition": ["H5N4F1S1", "H5N4", "H5N4F1"],
                "Structure": ["RES...", "RES...", "RES..."],
                "Sample1": [0.23, 0.13, 0.23],
                "Sample2": [0.43, 0.12, 0.32],
                "Sample3": [0.12, 0.32, 0.43],
            }
        )
        st.dataframe(example_input_df, hide_index=True)

    def _show_composition_mode_instruction(self) -> None:
        st.subheader("File format instruction")
        st.markdown(
            """The upload file should be a csv file with the following columns:
- The first column: "Composition", the glycan compositions. e.g. H5N4F1S1
- From the second column: glycan abundance in each sample, with column names as sample names."""
        )
        st.markdown("An example input file:")
        example_input_df = pd.DataFrame(
            {
                "Composition": ["H5N4F1S1", "H5N4", "H5N4F1"],
                "Sample1": [0.23, 0.13, 0.23],
                "Sample2": [0.43, 0.12, 0.32],
                "Sample3": [0.12, 0.32, 0.43],
            }
        )
        st.dataframe(example_input_df, hide_index=True)

    def _inquire_upload_file(self) -> None:
        st.subheader("File upload")
        st.file_uploader("Upload a csv file:", type="csv", key="input_file")
        if st.session_state.input_file is not None:
            st.success("File uploaded successfully!")
            st.session_state.input_df = load_data(st.session_state.input_file)
        else:
            st.session_state.input_df = None


@define
class CheckStructureStep(StreamlitStep):
    """Make sure that structure information is provided in the structure mode."""

    to_be_initialzed = {
        "structure_file": __BIND_TO_WEIGHT__,
        "db_name": None,
        "has_struc_col": None,
    }

    def check_needed(self) -> bool:
        return st.session_state.get("mode") == "structure"

    def check_passed(self) -> bool:
        return (
            self._has_struc_col()
            or st.session_state.get("structure_file") is not None
            or st.session_state.get("db_name") is not None
        )

    def run(self) -> None:
        self._update_has_struc_col()
        if st.session_state.has_struc_col is False:
            if not self.check_passed():
                self._show_no_struc_col_warning()
            self._inquire_structures()

    def _update_has_struc_col(self):
        st.session_state.has_struc_col = self._has_struc_col()

    def _has_struc_col(self) -> bool:
        input_df = st.session_state.input_df
        return "Structure" in input_df.columns

    def _show_no_struc_col_warning(self):
        msg = """
        The input file does not have a 'Structure' column.

        You can:
        1. Re-upload the input file with a 'Structure' column.
        2. Upload an additional structure file with two columns: 'Composition' and 'Structure'.
        3. Use the build-in database.
        4. Use the "composition" mode.

        (See the github page for more information.)
            """
        st.warning(msg, icon="⚠️")

    def _inquire_structures(self):
        def clear_state():
            st.session_state.structure_file = None
            st.session_state.db_name = None

        choice = st.selectbox(
            "Your choice:",
            ("Structure file", "Database"),
            index=0,
            on_change=clear_state,
        )
        if choice == "Structure file":
            self._inquire_structure_file()
        else:  # choice == "Database"
            self._inquire_db_name()

    def _inquire_structure_file(self):
        file = st.file_uploader("Upload a structure file:", type="csv")
        if file is not None:
            st.success("File uploaded successfully!")
            valid, reason = self._check_structure_file_validity(file)
            if valid:
                st.session_state.structure_file = file
            else:
                st.error(reason)
                st.stop()
        else:
            st.session_state.structure_file = None

    @staticmethod
    def _check_structure_file_validity(file) -> tuple[bool, str]:
        struc_df = load_data(file)
        input_df = st.session_state.input_df
        if (
            struc_df.shape[1] != 2
            or "Composition" not in struc_df.columns
            or "Structure" not in struc_df.columns
        ):
            msg = (
                "The structure file should only have two columns: "
                "'Composition' and 'Structure'."
            )
            return False, msg
        if missing := set(input_df["Composition"]) - set(struc_df["Composition"]):
            msg = f"These glycans don't have structures: {', '.join(missing)}"
            return False, msg
        return True, ""

    def _inquire_db_name(self):
        choice = st.selectbox(
            "Select a database:",
            ("serum", "IgG"),
            index=0,
        )
        st.session_state.db_name = choice


@define
class LoadAbundanceStep(StreamlitStep):
    """Load the abundance table."""

    to_be_initialzed = {"abund_df": None}

    def run(self) -> None:
        input_df = st.session_state.input_df
        input_df = input_df.set_index("Composition")
        if st.session_state.mode == "structure" and st.session_state.has_struc_col:
            abund_df = input_df.drop("Structure", axis=1).T
        else:
            abund_df = input_df.T
        st.session_state.abund_df = abund_df
        st.subheader("Abundance table confirmation")
        st.markdown(
            f"- Number of samples: **{abund_df.shape[0]}**\n"
            f"- Number of glycans: **{abund_df.shape[1]}**"
        )
        st.dataframe(abund_df)


@define
class UploadGroupFileStep(StreamlitStep):
    """Upload a group file."""

    to_be_initialzed = {"group_file": __BIND_TO_WEIGHT__, "group_series": None}

    def run(self) -> None:
        st.subheader("Group file upload")
        st.info(
            "Group information is essential for downstream analysis. "
            "However, if you just want to get the derived trait table, "
            "you can skip this step."
        )
        st.markdown(
            "1. The group file should be a CSV file with two columns: "
            "'sample' and 'type'.\n"
            "2. The sample names should be the same as the sample names "
            "in the previously uploaded file.\n"
            "3. There must be more than one group, "
            "and each group should have at least three samples."
        )
        st.markdown("An example group file:")
        example_group_df = pd.DataFrame(
            {
                "sample": [
                    "Sample1",
                    "Sample2",
                    "Sample3",
                    "Sample4",
                    "Sample5",
                    "Sample6",
                ],
                "type": ["Control", "Control", "Control", "Case", "Case", "Case"],
            }
        )
        st.dataframe(example_group_df, hide_index=True)
        file = st.file_uploader("Upload a group file:", type="csv")
        if file is not None:
            group_df = load_data(file)
            try:
                group_s = self._get_group_s(group_df)
            except ValueError as e:
                group_s = None
                st.error(e)
                st.stop()
            except Exception as e:
                group_s = None
                st.error("Unknown error occurred: " + str(e))
                st.stop()
            else:
                st.success("File uploaded successfully!")
            finally:
                st.session_state.group_series = group_s
        else:
            st.session_state.group_series = None
        if st.session_state.group_series is not None:
            st.subheader("Group table confirmation")
            st.dataframe(st.session_state.group_series.to_frame().value_counts())

    # @st.cache_data(ttl=600)
    @staticmethod
    def _get_group_s(group_df: pd.DataFrame) -> pd.Series:
        # 1. Check if the group dataframe has two columns: "Sample" and "Group".
        if (
            group_df.shape[1] != 2
            or "sample" not in group_df.columns
            or "type" not in group_df.columns
        ):
            msg = "The group dataframe should have two columns: 'sample' and 'type'."
            raise ValueError(msg)

        # 2. Check if the sample names are the same as the sample names
        # in the abundance dataframe.
        abund_df = st.session_state.abund_df
        if set(abund_df.index) != set(group_df["sample"]):
            msg = ("The sample names in the group dataframe should be the "
                   "same as the sample names in the abundance dataframe.")
            raise ValueError(msg)

        # 3. Check if there are more than one group.
        if len(group_df["type"].unique()) < 2:
            raise ValueError("There must be more than one group.")

        # 4. Check if each group has at least three samples.
        group_sizes = group_df["type"].value_counts()
        if (group_sizes < 3).any():
            raise ValueError("Each group should have at least three samples.")

        # 5. If all checks passed, return the group series.
        group_s = group_df.set_index("sample")
        group_s = group_s.squeeze()
        group_s = group_s.reindex(abund_df.index)
        return group_s


@define
class OptionStep(StreamlitStep):
    """Set options for GlyTrait."""

    to_be_initialzed = {
        "glycan_filter": __BIND_TO_WEIGHT__,
        "impute_method": __BIND_TO_WEIGHT__,
        "post_filtering": __BIND_TO_WEIGHT__,
        "corr_method": __BIND_TO_WEIGHT__,
        "corr_threshold": __BIND_TO_WEIGHT__,
        "sia_linkage": __BIND_TO_WEIGHT__,
    }

    def run(self) -> None:
        st.header("Options")
        col1, col2 = st.columns(2)

        col1.subheader("Preprocessing")
        col1.number_input(
            "Glycan filtering threshold:",
            min_value=0.5,
            max_value=1.0,
            value=1.0,
            step=0.05,
            help="Glycans with missing values in more than this proportion of samples "
            "will be filtered out. Setting to 1.00 means no filtering.",
            key="glycan_filter",
        )
        col1.selectbox(
            "Inputation method:",
            ("zero", "mean", "min", "median", "lod"),
            index=0,
            help="lod: limit of detection, 1/5 of the minimum value in each sample.",
            key="impute_method",
        )
        col1.selectbox(
            "Normalization method:",
            ("Total abundance",),
            index=0,
            help="We only support TA normalization now. "
            "The normalization method will not affect the values of derived traits.",
        )

        col2.subheader("Post-filtering")
        turn_on_post_filtering = st.session_state.get("abund_df").shape[0] >= 3
        if not turn_on_post_filtering:
            col2.info(
                "Post-filtering is disabled because there are less than 3 samples."
            )
        col2.checkbox(
            "Turn on post-filtering",
            value=turn_on_post_filtering,
            help="Filter out uninformative derived traits. "
            "See the github page for more information.",
            key="post_filtering",
            disabled=not turn_on_post_filtering,
        )
        col2.selectbox(
            "Post-filtering correlation method:",
            ("pearson", "spearman"),
            index=0,
            help="The correlation method used in post-filtering.",
            key="corr_method",
            disabled=not st.session_state.post_filtering,
        )
        col2.number_input(
            "Post-filtering correlation threshold:",
            min_value=-1.0,
            max_value=1.0,
            value=1.0,
            step=0.05,
            help="Derived traits with correlation coefficient higher than this threshold"
            " will be reduced (only one will be left)."
            "Setting to 1.00 to filter only traits with perfect correlation."
            "Setting to -1.00 means no filtering.",
            key="corr_threshold",
            disabled=not st.session_state.post_filtering,
        )

        st.subheader("Sialic acid linkage")
        st.checkbox(
            "Calculate sialic acid linkage traits",
            value=False,
            help="Caution: GlyTrait will not check if the input file has sialic acid linkage "
            "information before running. Please make sure the input file "
            "(or the additional structure file) has sialic acid linkage"
            " information. In the structure mode, the linkage of sialic acids should be "
            "specified in glycoCT strings. In the composition mode, all 'S' should be "
            "replaced by 'E' or 'L', for a2,6 and a2,3 linkage, respectively.",
            key="sia_linkage",
        )


@define
class CalculateTraitStep(StreamlitStep):
    """Calculate traits."""

    to_be_initialzed = {
        "formulas": None,
        "abund_df": None,
        "meta_property_df": None,
        "derived_trait_df": None,
    }

    def run(self) -> None:
        st.header("Calculate derived traits")
        if st.button("Run GlyTrait"):
            with st.spinner("Running..."):
                try:
                    self._run()
                except GlyTraitError as e:
                    st.error(e)
                    st.stop()
                else:
                    st.session_state.finished = True
        if st.session_state.finished:
            st.success("Calculation complete!")
            st.info("Please go the 'Results' page to see the results.")

    def _run(self):
        st.session_state.glycans = self._load_glycans()
        st.session_state.formulas = self._load_formulas()
        (
            st.session_state.abund_df,
            st.session_state.glycans,
        ) = self._preprocess()
        st.session_state.meta_property_df = self._calcu_meta_properties()
        st.session_state.derived_trait_df = self._calcu_derived_traits()
        (
            st.session_state.formulas,
            st.session_state.derived_trait_df,
        ) = self._post_filtering()

    def _load_glycans(self) -> list:
        input_df = st.session_state.input_df
        comp_strings = input_df["Composition"].tolist()
        if st.session_state.mode == "composition":
            glycans = load_compositions(
                comp_strings, sia_linkage=st.session_state.sia_linkage
            )
        else:
            if st.session_state.has_struc_col:
                glycans = load_glycans(comp_strings, input_df["Structure"])
            elif database := st.session_state.db_name:
                glycans = load_default_structures(database, comp_strings)
            elif structure_file := st.session_state.structure_file:
                glycans = read_structure_file(structure_file, comp_strings)
            else:
                raise ValueError("No structure information provided.")
        return glycans

    def _load_formulas(self) -> list:
        formulas = load_formulas(st.session_state.mode, None)
        if not st.session_state.sia_linkage:
            formulas = [f for f in formulas if f.sia_linkage is False]
        return formulas

    def _preprocess(self):
        abund_df_processed = preprocess_pipeline(
            st.session_state.abund_df,
            st.session_state.glycan_filter,
            st.session_state.impute_method,
        )
        glycans = st.session_state.glycans
        glycans_filtered = [g for g in glycans if g.name in abund_df_processed.columns]
        return abund_df_processed, glycans_filtered

    def _calcu_meta_properties(self):
        return build_meta_property_table(
            st.session_state.abund_df.columns.tolist(),
            st.session_state.glycans,
            st.session_state.mode,
            st.session_state.sia_linkage,
        )

    def _calcu_derived_traits(self):
        return calcu_derived_trait(
            st.session_state.abund_df,
            st.session_state.meta_property_df,
            st.session_state.formulas,
        )

    def _post_filtering(self):
        if st.session_state.post_filtering:
            formulas = st.session_state.formulas
            derived_traits = st.session_state.derived_trait_df
            formulas, derived_traits = filter_invalid(formulas, derived_traits)
            formulas, derived_traits = filter_colinearity(
                formulas,
                derived_traits,
                st.session_state.corr_threshold,
                st.session_state.corr_method,
            )
            return formulas, derived_traits


if st.session_state.finished is False:
    steps = [
        WelcomeStep(),
        ModeStep(),
        UploadFileStep(),
        CheckStructureStep(),
        LoadAbundanceStep(),
        UploadGroupFileStep(),
        OptionStep(),
        CalculateTraitStep(),
    ]
    MainPage(steps).run()
else:
    st.info("Please go to the 'Results' page to see the results.")
    st.button("Reset", on_click=lambda: st.session_state.clear())
