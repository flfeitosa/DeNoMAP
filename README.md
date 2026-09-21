# 🧪 De Novo Priority Plot

Interactive web app to prioritize *de novo* generated molecules in a 3D
**Novelty × Synthesizability × QED** space.

Upload a dataset, pick the SMILES column (and, optionally, an activity column to
color the points), and explore your library as an interactive 3D scatter — with
the molecule structure drawn beside the plot as you hover each point.

<p align="center">
  <em>Novelty (vs. ZINC) · Synthesizability (SA score) · QED — all on a 0–1, higher-is-better scale.</em>
</p>

---

## ✨ What it does

- **SMILES curation** — standardizes structures with RDKit's `MolStandardize`
  (cleanup, largest fragment, uncharge, reionize, canonical tautomer), drops
  invalid molecules and, optionally, duplicates (by InChIKey).
- **Novelty Score** — for each molecule, finds the nearest neighbor in a **ZINC**
  collection via the [SmallWorld API](https://sw.docking.org) and computes
  `1 − (Tc_ECFP4 + Tc_path) / 2` (0 = identical, 1 = fully novel).
- **Synthesizability** — Ertl & Schuffenhauer SA score, rescaled to 0–1
  (higher = easier to make).
- **QED** — quantitative estimate of drug-likeness (0–1).
- **Butina clustering** *(optional)* — groups similar molecules by ECFP4
  Tanimoto and scores only the **centroid** of each cluster, cutting the number
  of API requests to the number of clusters.
- **Interactive 3D plot** — Plotly scatter with a live structure panel on hover.
  Color points by any numeric (color bar) or categorical (legend) activity column.
- **Export** — download the fully scored dataset as CSV.

> ⚠️ The novelty score makes **one SmallWorld API request per molecule**, so it
> needs an internet connection and larger datasets take longer. A progress bar
> is shown for both the curation and the API steps, and results are cached so
> re-coloring the plot does not re-query.
>
> Above **500 molecules** the app warns you and enables Butina clustering
> automatically — see [Clustering](#-clustering-large-datasets) below.

---

## 💻 Run it locally

Requires **Python 3.9+**.

```bash
# 1. Clone the repository and 
git clone https://github.com/franciscolucasfo/DeNovo_Priority_Plot.git
cd DeNovo_Priority_Plot/

# 2. Create a environment
conda create -n denovo_plot python=3.11 -y
conda activate denovo_plot

# 3. Install the dependencies
pip install -r requirements.txt

# 4. Launch the app
streamlit run app.py
```

The app opens in your browser at `http://localhost:8501`.

---

## 📥 Input format

Any **CSV**, **TSV** or **Excel** file with at least one column of SMILES.

| SMILES                | pIC50 | id      |
|-----------------------|-------|---------|
| `CCOc1ccccc1C(=O)N`   | 7.2   | mol_001 |
| `Clc1ccc(cc1)C#N`     | 6.4   | mol_002 |

In the app you then choose:

- **SMILES column** *(required)* — used for scoring and depictions.
- **Activity column** *(optional)* — colors the points (numeric → color bar,
  categorical → legend).
- **ID column** *(optional)* — shown in the tooltip and structure panel.

---

## 🧩 Clustering large datasets

Novelty costs one SmallWorld request per molecule, which is what makes large
libraries slow. Clustering trades a little resolution for a large reduction in
requests: molecules are grouped with the **Butina** algorithm on the ECFP4
Tanimoto distance, and only the **centroid** of each cluster is scored and
plotted.

Enable it under **Clustering** in the sidebar:

- **Cluster with Butina and keep centroids** — off by default, but **checked
  automatically when the dataset has more than 500 molecules** (you can still
  uncheck it if you really want every molecule scored).
- **Similarity threshold** — the Tanimoto similarity above which two molecules
  share a cluster. Default **0.30**. Lower values give fewer, broader clusters
  and fewer API requests; higher values give more, tighter clusters.

The scored output gains two columns, `Cluster_ID` and `Cluster_Size`, so you can
tell how many molecules each plotted centroid stands for. Clustering runs after
curation, on the standardized structures; molecules that fail to parse are left
out of every cluster.

From Python, the same thing is available as:

```python
import utils

# One representative per cluster, with Cluster_ID / Cluster_Size appended.
centroids, stats = utils.cluster_centroids(df, smiles_col="SMILES",
                                           similarity_threshold=0.3)

# Or let calculate_axis do it as part of the pipeline.
scored = utils.calculate_axis(df, smiles_col="SMILES",
                              cluster=True, cluster_threshold=0.3)
```

---

## 📁 Project structure

```
De_novoplot_streamlit/
├── app.py               # Streamlit app (upload → select columns → plot)
├── utils.py             # Core logic: curation, clustering, scoring, axes, plotting
├── requirements.txt     # Python dependencies
├── packages.txt         # For streamlit cloud
└── README.md
```

---

## 🧬 How the axes are computed

| Axis              | Source                          | Range | Direction        |
|-------------------|---------------------------------|-------|------------------|
| `Novelty_Score`   | SmallWorld nearest ZINC neighbor | 0–1   | higher = novel   |
| `Synthesizability`| SA score, `(10 − SA) / 9`        | 0–1   | higher = easier  |
| `QED`             | RDKit `QED.qed`                  | 0–1   | higher = drug-like |

When clustering is enabled, two extra columns describe the grouping:

| Column         | Meaning                                                  |
|----------------|----------------------------------------------------------|
| `Cluster_ID`   | Index of the Butina cluster the centroid represents       |
| `Cluster_Size` | How many molecules that cluster contains                  |

