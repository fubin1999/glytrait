import dash_bio
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import roc_curve

from glytrait.analysis.diff import differential_analysis as diff_anal
from glytrait.analysis.roc import calcu_roc_auc
from glytrait.streamlit import utils


def differential_analysis():
    st.title("Differential analysis")
    groups = st.session_state.group_series
    if groups.nunique() == 2:
        st.info("Two groups are detected. **Mann-Whitney U test** will be performed.")
        st.info(
            "- **U-val**: U statistic of the Mann-Whitney U test;\n"
            "- **p-val**: p-value of the Mann-Whitney U test;\n"
            "- **p-val-adjusted**: p-value adjusted by Benjamini-Hochberg method;\n"
            "- **CLES**: common language effect size;\n",
            icon="💡",
        )
    st.info("You can click on the column header to sort the table.", icon="💡")
    result = diff_anal(st.session_state.derived_trait_df, groups)
    st.dataframe(
        result,
        column_config={
            "p-val": st.column_config.NumberColumn(
                "p-val", format="%.3e", width="medium"
            ),
            "p-val-adjusted": st.column_config.NumberColumn(
                "p-val-adjust", format="%.3e", width="medium"
            ),
        },
        use_container_width=True,
    )
    st.download_button(
        label="Download result as CSV",
        data=result.to_csv().encode("utf-8"),
        file_name="differential_analysis.csv",
    )

    st.header("Boxplot")
    st.write("Please select the trait you want to plot.")
    trait = st.selectbox("Select trait", st.session_state.derived_trait_df.columns)
    if description := utils.get_description(trait):
        st.info(description)
    col1, col2 = st.columns(2)
    fig1 = px.box(st.session_state.derived_trait_df, x=groups, y=trait)
    fig1.update_xaxes(title_text="")
    col1.plotly_chart(fig1, use_container_width=True)
    fig2 = px.violin(st.session_state.derived_trait_df, x=groups, y=trait)
    fig2.update_xaxes(title_text="")
    col2.plotly_chart(fig2, use_container_width=True)


def heatmap_clustering():
    st.title("Heatmap/clustering")
    with st.form("heatmap"):
        options = ["--All--"] + st.session_state.derived_trait_df.columns.tolist()
        selected = st.multiselect("Select traits", options, default=["--All--"])
        col1, col2 = st.columns(2)
        label_row = col1.checkbox("Label row", value=True)
        label_col = col1.checkbox("Label column", value=True)
        cluster_row = col2.checkbox("Cluster row", value=True)
        cluster_col = col2.checkbox("Cluster column", value=True)
        scale = st.radio("Select scale", ("row", "column", "none"), horizontal=True)
        st.form_submit_button("Update Heatmap")

    if "--All--" in selected:
        selected = st.session_state.derived_trait_df.columns.tolist()
    if len(selected) == 0:
        st.warning("Please select at least one trait.")
        st.stop()
    heatmap_df = st.session_state.derived_trait_df[selected]
    hidden_labels = []
    if not label_row:
        hidden_labels.append("row")
    if not label_col:
        hidden_labels.append("col")
    if cluster_row and cluster_col:
        cluster = "all"
    elif cluster_row and not cluster_col:
        cluster = "row"
    elif not cluster_row and cluster_col:
        cluster = "col"
    else:
        cluster = "none"
    fig = dash_bio.Clustergram(
        heatmap_df.T,
        row_labels=heatmap_df.columns.tolist(),
        column_labels=heatmap_df.index.tolist(),
        standardize=scale,
        cluster=cluster,
        hidden_labels=hidden_labels,
        height=800,
    )
    st.plotly_chart(fig, use_container_width=True)


def correlation_analysis():
    st.title("Correlation analysis")
    with st.form("corr_heatmap"):
        options = ["--All--"] + st.session_state.derived_trait_df.columns.tolist()
        selected = st.multiselect("Select traits", options, default=["--All--"])
        show_corr = st.checkbox("Show Corr. Coef.", value=False)
        corr_method = st.selectbox("Select corr. method", ("pearson", "spearman"))
        label_size = st.slider("Label size", 5, 20, 10)
        st.form_submit_button("Update Heatmap")

    if "--All--" in selected:
        selected = st.session_state.derived_trait_df.columns.tolist()
    if len(selected) == 0:
        st.warning("Please select at least one trait.")
        st.stop()
    corr_df = st.session_state.derived_trait_df[selected]
    corr_matrix = corr_df.corr(method=corr_method).round(2)
    fig = px.imshow(corr_matrix, text_auto=show_corr, color_continuous_scale="RdBu_r")
    fig.update_xaxes(showticklabels=True, ticks="")
    fig.update_yaxes(showticklabels=True, ticks="")
    fig.update_xaxes(side="top")
    fig.update_layout(
        xaxis=dict(
            tickvals=list(range(len(corr_matrix.columns))),
            ticktext=corr_matrix.columns,
            tickfont=dict(size=label_size),
        ),
        yaxis=dict(
            tickvals=list(range(len(corr_matrix.index))),
            ticktext=corr_matrix.index,
            tickfont=dict(size=label_size),
        ),
    )

    # Make the heatmap larger, but not the colorbar
    fig.update_layout(height=800)
    st.plotly_chart(fig, use_container_width=True)


def dimension_reduction():
    st.title("Dimension reduction")
    with st.form("dimension_reduction"):
        options = ["--All--"] + st.session_state.derived_trait_df.columns.tolist()
        selected = st.multiselect("Select traits", options, default=["--All--"])
        method = st.selectbox("Select method", ("PCA", "tSNE"))
        st.form_submit_button("Update Plot")

    if "--All--" in selected:
        selected = st.session_state.derived_trait_df.columns.tolist()
    if len(selected) == 0:
        st.warning("Please select at least one trait.")
        st.stop()

    df_selected = st.session_state.derived_trait_df[selected]
    if method == "PCA":
        pca = PCA(n_components=2)
        components = pca.fit_transform(df_selected)
        fig = px.scatter(components, x=0, y=1, color=st.session_state.group_series)
    elif method == "tSNE":
        tsne = TSNE(n_components=2)
        components = tsne.fit_transform(df_selected)
        fig = px.scatter(components, x=0, y=1, color=st.session_state.group_series)
    else:
        st.warning(f"{method} is not available yet.", icon="⚠️")
        st.stop()
    st.plotly_chart(fig, use_container_width=True)


def roc_analysis():
    st.title("ROC analysis")

    if st.session_state.group_series.nunique() != 2:
        st.warning(
            "ROC analysis with more than two groups is not available yet.", icon="⚠️"
        )
        st.stop()

    st.subheader("ROC Curve")
    with st.form("roc_analysis"):
        options = ["--All--"] + st.session_state.derived_trait_df.columns.tolist()
        selected = st.multiselect("Select traits", options, default=["--All--"])
        st.form_submit_button("Update Plot")

    if "--All--" in selected:
        selected = st.session_state.derived_trait_df.columns.tolist()
    if len(selected) == 0:
        st.warning("Please select at least one trait.")
        st.stop()

    fig = go.Figure()
    fig.add_shape(type="line", line=dict(dash="dash"), x0=0, x1=1, y0=0, y1=1)
    for trait in selected:
        model = LogisticRegression()
        model.fit(
            st.session_state.derived_trait_df[[trait]], st.session_state.group_series
        )
        scores = model.predict_proba(st.session_state.derived_trait_df[[trait]])[:, 1]
        fpr, tpr, _ = roc_curve(st.session_state.group_series, scores)
        fig.add_trace(go.Scatter(x=fpr, y=tpr, name=trait, mode="lines"))
    fig.update_layout(xaxis_title="FPR", yaxis_title="TPR")
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("ROC AUC")
    auc_df = calcu_roc_auc(
        st.session_state.derived_trait_df, st.session_state.group_series
    )
    auc_df["ROC AUC"] = auc_df["ROC AUC"].round(3)
    # Draw a dotplot of the 10 traits with the highest roc auc, using plotly
    fig = px.scatter(
        auc_df.sort_values("ROC AUC", ascending=False).head(10),
        x="ROC AUC",
        y="trait",
        size="ROC AUC",
        size_max=10,
        color="ROC AUC",
        color_continuous_scale="RdBu_r",
        text="trait",
        hover_name="trait",
        hover_data={"ROC AUC": True, "trait": False},
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(
        xaxis_title="ROC AUC",
        yaxis_title="",
        yaxis=dict(autorange="reversed"),
        height=500,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.download_button(
        label="Download result as CSV",
        data=auc_df.to_csv().encode("utf-8"),
        file_name="roc_auc.csv",
    )


@st.cache_data(ttl=600)
def combine_df(df1, df2):
    """Combine two dataframes with the same index."""
    return df1.merge(df2, left_index=True, right_index=True)


if st.session_state.finished is False:
    st.warning("Please upload a file and click 'Run' for further analysis.")
    st.stop()

st.sidebar.subheader("Select analysis")
analysis = st.sidebar.radio(
    "Select the type of analysis you want to perform",
    (
        "Differential analysis",
        "Heatmap/clustering",
        "Correlation analysis",
        "Dimension reduction",
        "ROC analysis",
    ),
)

analysis_dict = {
    "Differential analysis": differential_analysis,
    "Heatmap/clustering": heatmap_clustering,
    "Correlation analysis": correlation_analysis,
    "Dimension reduction": dimension_reduction,
    "ROC analysis": roc_analysis,
}

analysis_dict[analysis]()
