import tempfile
import zipfile
from io import BytesIO
from pathlib import Path

import pandas as pd
import streamlit as st

from glytrait.data_input import AbundanceLoader, GlycanLoader, GroupsLoader, GlyTraitInputData
from glytrait.exception import GlyTraitError
from glytrait import Experiment


BADGES = """[![PyPI - Version](https://img.shields.io/pypi/v/glytrait)](https://pypi.org/project/glytrait/) 
[![GitHub License](https://img.shields.io/github/license/fubin1999/glytrait)](https://github.com/fubin1999/glytrait/blob/main/LICENSE)
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
for details."""
MODE_HELP = """If you have structural information of the glycans, 
choose "structure" (recommended).
Otherwise, choose "composition"."""
ABUNDANCE_FILE_HELP = """A csv file with samples as rows and glycan names as columns.
The header of the first column should be "Sample", 
and the header of the other columns should be glycan IDs. 
Glycan IDs can be any string, e.g., the composition strings ("H3N4")."""
GLYCAN_FILE_HELP = """A csv file with two columns: "GlycanID" and "Structure" (or "Composition").
The "GlycanID" column should contain all glycan IDs in the abundance file. 
The "Structure" column should contain the structure strings of the glycans.
For now, only the GlycoCT format is supported.
In the "composition" mode, 
the second column should be "Composition" instead of "Structure", 
and the composition strings should be used instead of the structure strings. 
Condensed format ("H3N4F1S1") is supported.
"""
GROUP_FILE_HELP = """Upload a group file to enable differential analysis after trait calculation.
This is a csv file with two columns: "Sample" and "Group".
The "Sample" column should be the same as the abundance file.
The "Group" column contains the labels of samples, e.g. "Case", "Control".
"""
GLYCAN_FILTER_HELP = """Glycans with missing value proportion larger than 
this threshold will be removed.
Setting to 1 means no glycan will be removed.
Setting to 0 means glycans with only 1 missing value will be removed."""
IMPUTE_METHOD_HELP = """- "min": impute missing values by the minimum value of a glycan 
within all samples.
- "mean": impute missing values by the mean value of a glycan within all samples.
- "median": impute missing values by the median value of a glycan within all samples.
- "zero": impute missing values by 0.
- "lod": impute missing values by the limit of detection (LOD) of the equipment. 
The LOD of a glycan is defined as the minimum value of the glycan 
within all samples divided by 5."""
POST_FILTERING_HELP = """Post-filtering remove invalid and redundant traits.
This threshold controls the Pearson's correlation coefficient value to 
trigger the collinearity filtering.
Refer to [README](https://github.com/fubin1999/glytrait/blob/main/README.md#post-filtering)
for details,
or just use "1.0" to filter only traits with a perfect collinearity.
"""


# ========== S I D E B A R ==========
with st.sidebar:
    st.image("img/logo.png")
    st.markdown(BADGES)
    st.header("About")
    st.markdown(SIDEBAR_TEXT)


# ========== T I T L E ==========
st.title("GlyTrait")

# ========== I N P U T ==========
st.header("Input")
st.write(INPUT_WELCOME)
input_c = st.container()
mode = input_c.selectbox("Mode", ["structure", "composition"], help=MODE_HELP)

abundance_file = input_c.file_uploader("Abundance File", help=ABUNDANCE_FILE_HELP, type="csv")
if abundance_file:
    try:
        abundance_df = AbundanceLoader().load(pd.read_csv(abundance_file))
    except GlyTraitError as e:
        st.error(str(e))
        st.stop()
else:
    abundance_df = None

glycan_file = input_c.file_uploader(
    "Structure (Composition) File", help=GLYCAN_FILE_HELP, type="csv"
)
if glycan_file:
    try:
        glycans = GlycanLoader(mode=mode).load(pd.read_csv(glycan_file))
    except GlyTraitError as e:
        st.error(str(e))
        st.stop()
else:
    glycans = None

group_file = input_c.file_uploader("Group File (optional)", help=GROUP_FILE_HELP, type="csv")
if group_file:
    try:
        groups = GroupsLoader().load(pd.read_csv(group_file))
    except GlyTraitError as e:
        st.error(str(e))
        st.stop()
else:
    groups = None

if abundance_df is None or glycans is None:
    st.stop()

input_data = GlyTraitInputData(
    abundance_table=abundance_df,
    glycans=glycans,
    groups=groups,
)

# ========== C O N F I G ==========
st.header("Configuration")
col1, col2 = st.columns(2)
glycan_filter_threshold = col1.number_input(
    "Glycan Filter Threshold",
    min_value=0.,
    max_value=1.,
    step=0.1,
    value=1.,
    help=GLYCAN_FILTER_HELP,
)
impute_method = col2.selectbox(
    "Imputation Method",
    ["min", "mean", "median", "zero", "lod"],
    help=IMPUTE_METHOD_HELP,
)
post_filter_threshold = col1.number_input(
    "Post-filtering Threshold",
    min_value=0.,
    max_value=1.,
    step=0.1,
    value=1.,
    help=POST_FILTERING_HELP,
)


# ========== R U N ==========
def save_df_as_csv(df, filename, dirpath):
    """Save to csv."""
    df.to_csv(Path(dirpath) / filename, index=True)


def make_zipfile(dirpath):
    """将指定目录下的文件压缩成zip文件，并返回对应的BytesIO对象"""
    s = BytesIO()
    with zipfile.ZipFile(s, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in Path(dirpath).glob('*'):
            zipf.write(file, arcname=file.name)
    s.seek(0)
    return s


st.markdown("---")
if st.button("Run GlyTrait"):
    exp = Experiment(input_data)
    try:
        exp.run_workflow(
            filter_max_na=glycan_filter_threshold,
            impute_method=impute_method,
            corr_threshold=post_filter_threshold,
        )
    except GlyTraitError as e:
        st.error(str(e))
        st.stop()
    else:
        st.success("Succeed! Click the button below to download the results.")
else:
    exp = None

if exp:
    with tempfile.TemporaryDirectory() as tmp_dir:
        save_df_as_csv(exp.derived_trait_table, "derived_traits.csv", tmp_dir)
        save_df_as_csv(
            exp.filtered_derived_trait_table, "filtered_derived_traits.csv", tmp_dir
        )
        save_df_as_csv(exp.processed_abundance_table, "processed_abundance.csv", tmp_dir)
        save_df_as_csv(exp.meta_property_table, "meta_properties.csv", tmp_dir)
        if groups:
            if len(exp.groups.unique()) == 2:
                t_test_result = exp.diff_results["t_test"]
                save_df_as_csv(t_test_result, "t_test_result.csv", tmp_dir)
            else:
                anova_result = exp.diff_results["anova"]
                post_hoc_result = exp.diff_results["post_hoc"]
                save_df_as_csv(anova_result, "anova_result.csv", tmp_dir)
                save_df_as_csv(post_hoc_result, "post_hoc_result.csv", tmp_dir)
        zip_bytes = make_zipfile(tmp_dir)
        st.download_button(
            "Download Results",
            data=zip_bytes,
            file_name="GlyTrait_results.zip",
            mime="application/zip",
        )
