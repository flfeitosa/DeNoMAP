"""Streamlit app for DeNoMAP.

The user uploads a dataset, picks the SMILES column and (optionally) an activity
column, and gets the interactive 3D priority plot with a molecule panel on hover.
"""

import io

import pandas as pd
import streamlit as st

import utils

st.set_page_config(page_title="DeNoMAP", page_icon="🧪", layout="wide")

st.title("🧪 DeNoMAP")
st.caption(
    "Upload a dataset of molecules, choose the SMILES column and an optional "
    "activity column, and explore them in the Novelty / Synthesizability / QED "
    "space. Novelty is computed against ZINC through the SmallWorld API, so it "
    "needs one request per molecule and an internet connection."
)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_dataframe(uploaded_file):
    """Read the uploaded file into a DataFrame, guessing CSV / Excel by suffix."""
    name = uploaded_file.name.lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(uploaded_file)
    if name.endswith((".tsv", ".txt")):
        return pd.read_csv(uploaded_file, sep="\t")
    return pd.read_csv(uploaded_file)


uploaded = st.file_uploader(
    "Dataset (CSV, TSV or Excel)", type=["csv", "tsv", "txt", "xlsx", "xls"]
)

if uploaded is None:
    st.info("Upload a file to get started. It must contain at least a SMILES column.")
    st.stop()

try:
    df = load_dataframe(uploaded)
except Exception as exc:
    st.error(f"Could not read the file: {exc}")
    st.stop()

if df.empty:
    st.error("The uploaded file has no rows.")
    st.stop()

st.success(f"Loaded {len(df)} rows and {len(df.columns)} columns.")

# A set this large is painful to score one molecule at a time, so the warning
# below also pre-checks the clustering box instead of relying on the user to
# remember. They can still turn it off in the sidebar.
large_dataset = len(df) > utils.LARGE_DATASET_SIZE

if large_dataset:
    st.warning(
        f"**WARNING:** this dataset has {len(df)} compounds. Novelty needs one "
        f"SmallWorld API request per molecule, so scoring a set this large can "
        f"take a very long time. **BitBIRCH clustering has been enabled "
        f"automatically** in the sidebar: only one centroid per cluster is sent "
        f"to the API, which keeps the request count (and the runtime) "
        f"manageable. Uncheck it there if you really want every molecule scored."
    )

with st.expander("Preview data", expanded=False):
    st.dataframe(df.head(20), width="stretch")


# ---------------------------------------------------------------------------
# Column selection and options
# ---------------------------------------------------------------------------

columns = list(df.columns)

# Try to preselect a likely SMILES column by name.
def _guess_smiles_col(cols):
    for c in cols:
        if str(c).strip().lower() in ("smiles", "smi", "canonical_smiles"):
            return cols.index(c)
    return 0


col_a, col_b = st.columns(2)
with col_a:
    smiles_col = st.selectbox(
        "SMILES column", columns, index=_guess_smiles_col(columns)
    )
with col_b:
    activity_options = ["(none)"] + [c for c in columns if c != smiles_col]
    activity_choice = st.selectbox(
        "Activity column (optional, colors the points)", activity_options, index=0
    )
activity_col = None if activity_choice == "(none)" else activity_choice

with st.sidebar:
    st.header("Options")
    id_choice = st.selectbox(
        "ID column (optional, shown in tooltip)",
        ["(row index)"] + [c for c in columns if c != smiles_col],
        index=0,
    )
    id_col = None if id_choice == "(row index)" else id_choice

    curate = st.checkbox(
        "Curate SMILES", value=True,
        help="Standardize structures and drop invalid molecules before scoring.",
    )
    drop_duplicates = st.checkbox(
        "Drop duplicates", value=True,
        help="Keep one molecule per InChIKey (only applies when curating).",
    )
    structures = st.checkbox(
        "Draw structures on hover", value=True,
        help="Render a molecule image next to the plot. Turn off for large sets.",
    )

    st.subheader("Clustering")
    cluster = st.checkbox(
        "Cluster with BitBIRCH and keep centroids", value=large_dataset,
        help="Group similar molecules and score only one representative "
             "(the centroid) per cluster, cutting the number of API requests. "
             f"Checked automatically above {utils.LARGE_DATASET_SIZE} molecules.",
    )
    cluster_threshold = st.number_input(
        "Similarity threshold", min_value=0.0,
        max_value=utils.MAX_CLUSTER_THRESHOLD,
        value=utils.DEFAULT_CLUSTER_THRESHOLD, step=0.05, format="%.2f",
        disabled=not cluster,
        help="Tanimoto similarity (ECFP4) above which two molecules fall in the "
             "same cluster. Lower values give fewer, broader clusters - and so "
             "fewer API requests.",
    )

    st.subheader("Novelty / SmallWorld")
    try:
        db_names = list(utils.sw_databases().keys())
        default_idx = db_names.index(utils.DEFAULT_DB) if utils.DEFAULT_DB in db_names else 0
        db = st.selectbox("ZINC collection", db_names, index=default_idx)
    except Exception:
        db = utils.DEFAULT_DB
        st.caption(f"Could not list collections; using {db}.")
    dist = st.text_input("Distance range (dist)", value="0-16")
    timeout = st.number_input("Request timeout (s)", min_value=10, max_value=600, value=300)

    # Only affects the plot and the exported flag, so applying new thresholds
    # reuses the cached scores and never re-queries the API. The form holds the
    # changes back until "Apply" is pressed, so dragging a slider does not
    # redraw the plot on every step.
    st.subheader("Desirability zone")
    with st.form("desirability_form", border=False):
        show_zone = st.checkbox(
            "Show desirability zone", value=True,
            help="Outline a box from the thresholds below up to 1 on "
                 "every axis, and flag the molecules inside it in the "
                 f"{utils.DESIRABILITY_COL} column of the CSV.",
        )
        zone_thresholds = tuple(
            st.slider(f"{axis} ≥", min_value=0.0, max_value=1.0,
                      value=utils.DEFAULT_DESIRABILITY, step=0.05,
                      key=f"zone_{axis}")
            for axis in utils.DEFAULT_AXES
        )
        st.form_submit_button(
            "Apply desirability zone", width="stretch",
            help="Reclassify the molecules with the new thresholds. The scores "
                 "are not recomputed.",
        )


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

# Cache the expensive axis computation so re-coloring the plot does not re-query
# the SmallWorld API. Keyed on the actual data plus the parameters that change it.
@st.cache_data(show_spinner=False)
def compute_axes(df, smiles_col, curate, drop_duplicates, cluster, cluster_threshold,
                 db, dist, timeout):
    progress_bar = st.progress(0.0, text="Starting...")

    def on_curate_progress(done, total):
        progress_bar.progress(done / total, text=f"Curating SMILES... {done}/{total}")

    def on_cluster_progress(done, total):
        progress_bar.progress(done / total,
                              text=f"Clustering with BitBIRCH... step {done}/{total}")

    def on_progress(done, total):
        progress_bar.progress(done / total, text=f"Querying SmallWorld API... {done}/{total}")

    result = utils.calculate_axis(
        df, smiles_col=smiles_col, curate=curate, drop_duplicates=drop_duplicates,
        cluster=cluster, cluster_threshold=float(cluster_threshold),
        db=db, dist=dist, timeout=int(timeout), progress=on_progress,
        curate_progress=on_curate_progress, cluster_progress=on_cluster_progress,
    )
    progress_bar.empty()
    return result


if st.button("Compute axes and plot", type="primary"):
    n = len(df)
    spinner_msg = (
        f"Curating, clustering and then scoring {n} molecule(s) - one SmallWorld "
        f"request per cluster. Watch the progress bar for the current stage."
        if cluster else
        f"Scoring {n} molecule(s)... this queries ZINC once per molecule."
    )
    with st.spinner(spinner_msg):
        try:
            st.session_state.result = compute_axes(
                df, smiles_col, curate, drop_duplicates, cluster, cluster_threshold,
                db, dist, timeout
            )
        except Exception as exc:
            st.error(f"Failed to compute the axes: {exc}")
            st.stop()

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

result = st.session_state.get("result")
if result is None:
    st.stop()

# After curation the SMILES column is renamed to "SMILES".
plot_smiles_col = "SMILES" if "SMILES" in result.columns else smiles_col
# id_col may have survived curation; fall back to the row index otherwise.
plot_id_col = id_col if (id_col and id_col in result.columns) else None
plot_activity_col = activity_col if (activity_col and activity_col in result.columns) else None

n_failed = int(result["Novelty_Score"].isna().sum())
n_ok = int(result["Novelty_Score"].notna().sum())
c1, c2, c3, c4 = st.columns(4)
c1.metric("Molecules plotted", n_ok)
c2.metric("Novelty lookup failed", n_failed)
c3.metric("QED median", f"{result['QED'].median():.3f}")
if "Cluster_Size" in result.columns:
    c4.metric("Clusters", int(result["Cluster_ID"].nunique()),
              help="One centroid per cluster was scored instead of every molecule.")

if show_zone:
    # Work on a copy so the cached result is never mutated between reruns.
    result = result.assign(**{utils.DESIRABILITY_COL: utils.desirability_mask(
        result, zone_thresholds)})
    st.caption(f"**{int(result[utils.DESIRABILITY_COL].sum())}** of {n_ok} "
               "plotted molecule(s) fall inside the desirability zone.")

try:
    fig = utils.denovo_plot(
        result, smiles_col=plot_smiles_col, activity_col=plot_activity_col,
        id_col=plot_id_col, structures=structures,
        desirability=zone_thresholds if show_zone else None,
    )
except Exception as exc:
    st.error(f"Failed to build the plot: {exc}")
    st.stop()

html = utils.plot_with_panel_html(fig, include_plotlyjs="cdn")
st.components.v1.html(html, height=740, scrolling=True)

# ---------------------------------------------------------------------------
# Download the scored table
# ---------------------------------------------------------------------------

csv_buffer = io.StringIO()
result.to_csv(csv_buffer, index=False)
st.download_button(
    "Download scored dataset (CSV)",
    data=csv_buffer.getvalue(),
    file_name="denovo_priority_scored.csv",
    mime="text/csv",
)

with st.expander("Scored data", expanded=False):
    st.dataframe(result, width="stretch")
