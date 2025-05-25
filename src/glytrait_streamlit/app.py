# ========== IMPORTS ==========
import tempfile
import zipfile
from io import BytesIO
from pathlib import Path
from contextlib import contextmanager

import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

from glytrait.exception import GlyTraitError
from glytrait import Experiment

# ========== CONSTANTS ==========
BADGES = """[![PyPI - Version](https://img.shields.io/pypi/v/glytrait)](https://pypi.org/project/glytrait/) 
[![GitHub License](https://img.shields.io/github/license/fubin1999/glytrait)](https://github.com/fubin1999/glytrait/blob/main/LICENSE)
[![Static Badge](https://img.shields.io/badge/Github-grey?logo=github)](https://github.com/fubin1999/glytrait)
"""

SIDEBAR_TEXT = """
Q: What is GlyTrait?

A: GlyTrait is a tool to calculate derived traits for N-glycomic data.

Q: Sounds cool! So..., what are derived traits again?

A: Well, derived traits are artificial variables that summarize certain aspects of 
the glycome. 
For example, the proportion of core-fucosylated glycans, 
the average number of sialic acids per glycan, 
or the proportion of bisected glycans within bi-antennary complex glycans, etc. 
Derived traits are more biologically relavant 
and have been used a lot in the glycomics community.

Q: So, GlyTrait does the dirty work for me, 
saving my time and energy for more interesting analysis?

A: You bet!
"""

INPUT_WELCOME = """Upload an abundance file and a structure (or composition) file to start.
See [README](https://github.com/fubin1999/glytrait/blob/main/README.md#input-file-format) 
for details.
Or use our example files below to have a quick start.
"""

HELP_TEXTS = {
    "mode": """If you have structural information of the glycans, 
choose "structure" (recommended).
Otherwise, choose "composition".""",
    
    "abundance_file": """A csv file with samples as rows and glycan names as columns.
The header of the first column should be "Sample", 
and the header of the other columns should be glycan IDs. 
Glycan IDs can be any string, e.g., the composition strings ("H3N4").""",
    
    "glycan_file": """A csv file with two columns: "GlycanID" and "Structure" (or "Composition").
The "GlycanID" column should contain all glycan IDs in the abundance file. 
The "Structure" column should contain the structure strings of the glycans.
For now, only the GlycoCT format is supported.
In the "composition" mode, 
the second column should be "Composition" instead of "Structure", 
and the composition strings should be used instead of the structure strings. 
Condensed format ("H3N4F1S1") is supported.""",
    
    "group_file": """Upload a group file to enable differential analysis after trait calculation.
This is a csv file with two columns: "Sample" and "Group".
The "Sample" column should be the same as the abundance file.
The "Group" column contains the labels of samples, e.g. "Case", "Control".""",
    
    "glycan_filter": """Glycans with missing value proportion larger than 
this threshold will be removed.
Setting to 1 means no glycan will be removed.
Setting to 0 means glycans with only 1 missing value will be removed.""",
    
    "impute_method": """- "min": impute missing values by the minimum value of a glycan 
within all samples.
- "mean": impute missing values by the mean value of a glycan within all samples.
- "median": impute missing values by the median value of a glycan within all samples.
- "zero": impute missing values by 0.
- "lod": impute missing values by the limit of detection (LOD) of the equipment. 
The LOD of a glycan is defined as the minimum value of the glycan 
within all samples divided by 5.""",
    
    "post_filtering": """Post-filtering remove invalid and redundant traits.
This threshold controls the Pearson's correlation coefficient value to 
trigger the collinearity filtering.
Refer to [README](https://github.com/fubin1999/glytrait/blob/main/README.md#post-filtering)
for details,
or just use "1.0" to filter only traits with a perfect collinearity.""",
    
    "sia_linkage": """If True, the linkage of sialic acid will be considered.
Note that in this case all sialic acid residues should have linkage information."""
}

# ========== UTILITY FUNCTIONS ==========
@st.cache_resource
def read_example_file(filename):
    """Read an example file from the data folder."""
    return pd.read_csv(filename)

@contextmanager
def capture_glytrait_error():
    """Capture GlyTraitError, show it in the streamlit app, and stop the app."""
    try:
        yield
    except GlyTraitError as e:
        st.error(str(e))
        st.stop()

def save_df_as_csv(df, filename, dirpath):
    """Save to csv."""
    df.to_csv(Path(dirpath) / filename, index=True)

def make_zipfile(dirpath):
    """将指定目录下的文件压缩成zip文件，并返回对应的BytesIO对象"""
    s = BytesIO()
    with zipfile.ZipFile(s, "w", zipfile.ZIP_DEFLATED) as zipf:
        for file in Path(dirpath).glob("*"):
            zipf.write(file, arcname=file.name)
    s.seek(0)
    return s

def prepare_zip(exp, tmp_dir):
    """Save the results to a temporary directory and return the zip file as bytes."""
    save_df_as_csv(exp.derived_trait_table, "derived_traits.csv", tmp_dir)
    save_df_as_csv(
        exp.filtered_derived_trait_table, "filtered_derived_traits.csv", tmp_dir
    )
    save_df_as_csv(exp.processed_abundance_table, "processed_abundance.csv", tmp_dir)
    save_df_as_csv(exp.meta_property_table, "meta_properties.csv", tmp_dir)
    if exp.groups is not None:
        if len(exp.groups.unique()) == 2:
            t_test_result = exp.diff_results["t_test"]
            save_df_as_csv(t_test_result, "t_test_result.csv", tmp_dir)
        else:
            anova_result = exp.diff_results["anova"]
            post_hoc_result = exp.diff_results["post_hoc"]
            save_df_as_csv(anova_result, "anova_result.csv", tmp_dir)
            save_df_as_csv(post_hoc_result, "post_hoc_result.csv", tmp_dir)
    return make_zipfile(tmp_dir)

# ========== UI COMPONENTS ==========
def render_sidebar():
    """Render the sidebar with logo and about section."""
    with st.sidebar:
        st.image("img/logo.png")
        st.markdown(BADGES)
        st.header("About")
        st.markdown(SIDEBAR_TEXT)

def render_input_section():
    """Render the input section with file uploaders and configuration."""
    st.header("Input")
    st.write(INPUT_WELCOME)
    
    input_c = st.container()
    mode = input_c.selectbox("Mode", ["structure", "composition"], help=HELP_TEXTS["mode"])

    # Upload the abundance file
    abundance_file = input_c.file_uploader("Abundance file", type="csv")
    with input_c.expander("File format instructions"):
        st.markdown(HELP_TEXTS["abundance_file"])
        st.dataframe(abundance_example.head(5), hide_index=True)
        st.download_button(
            "Download Example Abundance File",
            data=abundance_example.to_csv(index=False).encode(),
            file_name="example_abundance.csv",
            mime="text/csv",
        )

    # Upload the structure (composition) file
    glycan_file = input_c.file_uploader("Glycan file", type="csv")
    with input_c.expander("File format instructions"):
        st.markdown(HELP_TEXTS["glycan_file"])
        example_df = structures_example if mode == "structure" else compositions_example
        st.dataframe(example_df.head(5), hide_index=True)
        st.download_button(
            f"Download Example {mode.capitalize()} File",
            data=example_df.to_csv(index=False).encode(),
            file_name=f"example_{mode}.csv",
            mime="text/csv",
        )

    # Upload the group file
    group_file = input_c.file_uploader("Group file", type="csv")
    with input_c.expander("File format instructions"):
        st.markdown(HELP_TEXTS["group_file"])
        st.dataframe(groups_example.head(5), hide_index=True)
        st.download_button(
            "Download Example Group File",
            data=groups_example.to_csv(index=False).encode(),
            file_name="example_groups.csv",
            mime="text/csv",
        )

    return mode, abundance_file, glycan_file, group_file

def render_config_section():
    """Render the configuration section with various settings."""
    st.header("Configuration")
    col1, col2 = st.columns(2)
    
    glycan_filter_threshold = col1.number_input(
        "Glycan Filter Threshold",
        min_value=0.0,
        max_value=1.0,
        step=0.1,
        value=1.0,
        help=HELP_TEXTS["glycan_filter"],
    )
    
    impute_method = col2.selectbox(
        "Imputation Method",
        ["min", "mean", "median", "zero", "lod"],
        help=HELP_TEXTS["impute_method"],
    )
    
    post_filter_threshold = col1.number_input(
        "Post-filtering Threshold",
        min_value=0.0,
        max_value=1.0,
        step=0.1,
        value=1.0,
        help=HELP_TEXTS["post_filtering"],
    )
    
    sia_linkage = col2.selectbox(
        "Consider Sialic Acid Linkage",
        [True, False],
        index=1,
        help=HELP_TEXTS["sia_linkage"],
    )
    
    return glycan_filter_threshold, impute_method, post_filter_threshold, sia_linkage

def render_analysis_section(exp):
    """Render the analysis section for trait visualization."""
    st.markdown("---")
    st.header("Analysis")
    
    # Get all available traits
    all_traits = exp.filtered_derived_trait_table.columns.tolist()
    
    # Let user select a trait to visualize
    selected_trait = st.selectbox(
        "Select a trait to visualize",
        options=all_traits,
        help="Select a trait to see its distribution"
    )
    
    if not selected_trait:
        st.info("Please select a trait to visualize")
        return
    
    # Create visualization based on whether group information is available
    if exp.groups is not None:
        # Create box plot when group information is available
        fig = px.box(
            exp.filtered_derived_trait_table,
            y=selected_trait,
            color=exp.groups,
            title=f"Distribution of {selected_trait} by Group",
            labels={"y": selected_trait, "color": "Group"},
        )
    else:
        # Create histogram when no group information is available
        fig = px.histogram(
            exp.filtered_derived_trait_table,
            x=selected_trait,
            title=f"Distribution of {selected_trait}",
            labels={"x": selected_trait},
            nbins=30,  # 设置直方图的箱数
        )
    
    # 更新布局
    fig.update_layout(
        showlegend=True,
        xaxis_title=selected_trait,
        yaxis_title="Count" if exp.groups is None else selected_trait,
    )
    
    # 显示图表
    st.plotly_chart(fig, use_container_width=True)

def render_results_section(exp):
    """Render the results section with download button."""
    st.markdown("---")
    st.header("Results")
    st.success("Succeed! Click the button below to download the results.")
    
    # Store the experiment object in session state
    st.session_state.exp = exp
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_bytes = prepare_zip(exp, tmp_dir)
        st.download_button(
            "Download Results",
            data=zip_bytes,
            file_name="GlyTrait_results.zip",
            mime="application/zip",
        )
    
    # Add analysis section
    render_analysis_section(exp)

# ========== MAIN PROGRAM ==========
def main():
    """Main program entry point."""
    # Load example files
    global abundance_example, structures_example, compositions_example, groups_example
    abundance_example = read_example_file("example_data/abundance.csv")
    structures_example = read_example_file("example_data/structures.csv")
    compositions_example = read_example_file("example_data/compositions.csv")
    groups_example = read_example_file("example_data/groups.csv")

    # Render UI components
    render_sidebar()
    st.title("GlyTrait")
    
    # Get input files and configuration
    mode, abundance_file, glycan_file, group_file = render_input_section()
    if abundance_file is None or glycan_file is None:
        st.stop()
        
    glycan_filter_threshold, impute_method, post_filter_threshold, sia_linkage = render_config_section()

    # Run analysis when button is clicked
    if st.button("Run GlyTrait"):
        with st.spinner("Running..."):
            exp = Experiment(
                abundance_file=abundance_file,
                glycan_file=glycan_file,
                group_file=group_file,
                mode=mode,
                sia_linkage=sia_linkage,
            )
            with capture_glytrait_error():
                exp.run_workflow(
                    filter_max_na=glycan_filter_threshold,
                    impute_method=impute_method,
                    corr_threshold=post_filter_threshold,
                )
        render_results_section(exp)
    
    # Re-render results section if it exists in session state
    elif "exp" in st.session_state:
        render_results_section(st.session_state.exp)

if __name__ == "__main__":
    main()
