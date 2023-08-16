import pandas as pd
import plotly.express as px
import streamlit as st

if st.session_state.finished is False:
    st.warning("Please upload a file and click 'Run' to see the results.")
    st.stop()

st.title("Results")

st.markdown(
    "Glycan derived traits are summary statistics of glycans, "
    "calculated based on the shared glycan structural features. "
    "Compared to glycan abundances, "
    "which are sometimes called 'the direct traits', "
    "derived glycosylation traits are biologically more related to "
    "activities of specific enzymes in the glycosylation pathway and "
    "underlying genetic polymorphisms, "
    "making them beneficial for understanding the functional relevance "
    "of obtained results."
)
st.markdown(
    "In a data science perspective, derived traits are just new features or variables, "
    "more relevant to our research question than the original features. "
    "As so, they can be used just as glycan abundances to perform "
    "further analysis, such as clustering, differential analysis, "
    "and machine learning."
)

st.header("Overview")
col1, col2, col3 = st.columns(3)
col1.metric("Number of samples", st.session_state.abund_df.shape[0])
col2.metric("Number of glycans", st.session_state.abund_df.shape[1])
col3.metric("Number of derived traits", st.session_state.derived_trait_df.shape[1])


def get_summary_table(trait_df: pd.DataFrame, formulas: list):
    """Get the summary table of the derived trait table."""
    summary_table = trait_df.describe().T
    summary_table = summary_table.drop("count", axis=1)
    formula_map = {f.name: f.description for f in formulas}
    summary_table["Description"] = summary_table.index.map(formula_map)
    # Reorder the columns, move "Description" to the front
    cols = summary_table.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    summary_table = summary_table[cols]
    return summary_table


st.header("Trait details")
# Get the min, max, mean, and median of each trait from the derived trait table
summary_table = get_summary_table(
    st.session_state.derived_trait_df, st.session_state.formulas
)
trait_selected = st.selectbox("Select a trait", summary_table.index)
st.info("**Definition**: " + summary_table.loc[trait_selected, "Description"])
col1, col2, col3 = st.columns(3)
col1.metric("Mean", round(summary_table.loc[trait_selected, "mean"], 5))
col2.metric("Median", round(summary_table.loc[trait_selected, "50%"], 5))
col3.metric("Std", round(summary_table.loc[trait_selected, "std"], 5))
col1.metric("Min", round(summary_table.loc[trait_selected, "min"], 5))
col2.metric("Max", round(summary_table.loc[trait_selected, "max"], 5))
histgram = px.histogram(st.session_state.derived_trait_df, x=trait_selected)
st.plotly_chart(histgram, use_container_width=True)

st.header("Download results")
st.success(
    "Feel free to download the results GlyTrait generated and start doing "
    "downstream analysis, "
    "or you can click on the 'Analysis' page on the left to perform some "
    "preliminary analysis we provide."
)

st.subheader("Derived trait table")
st.info(
    "The main result produced by GlyTrait. "
    "Feel free to download this table for further analysis,"
    "or use our 'Analysis' page to perform a preliminary analysis."
)
st.dataframe(st.session_state.derived_trait_df, height=210)
st.download_button(
    "Download",
    st.session_state.derived_trait_df.to_csv(),
    file_name="derived_trait_table.csv",
)

st.subheader("Glycan abundance table")
st.info("The glycan abundance table after preprocessing.")
st.dataframe(st.session_state.abund_df, height=210)
st.download_button(
    "Download",
    st.session_state.abund_df.to_csv(),
    file_name="glycan_abundance_table.csv",
)

st.subheader("Derived trait summary table")
st.info("Contains definitions and summary statistics of each trait.")
st.dataframe(summary_table, height=210)
st.download_button(
    "Download",
    summary_table.to_csv(),
    file_name="derived_trait_summary_table.csv",
)

st.subheader("Glycan meta property table")
st.info(
    "The meta properties (intermediate results used to calculate derived traits) "
    "of each glycan."
    "Also provides important information about the glycan structure. "
    "Useful for glycan count based analysis."
)
st.dataframe(st.session_state.meta_property_df, height=210)
st.download_button(
    "Download",
    st.session_state.meta_property_df.to_csv(),
    file_name="glycan_meta_property_table.csv",
)
