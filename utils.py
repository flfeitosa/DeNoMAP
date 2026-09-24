"""DeNoMAP - core functions.

Curates SMILES, computes the three priority-plot axes (Novelty against ZINC,
Synthesizability from the SA score, and QED), and builds the interactive 3D
plot. Extracted from the De_Novo_Priority_Plot notebook so it can be reused by
the Streamlit app in app.py.
"""

import math
import os
import sys
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

import numpy as np
import pandas as pd
import requests

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import QED, RDConfig, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

import plotly.graph_objects as go

RDLogger.DisableLog("rdApp.*")  # Sometimes there are errors that could flood the output.

# sascorer ships inside RDKit's Contrib directory, which is not on sys.path by default.
_SA_SCORE_DIR = os.path.join(RDConfig.RDContribDir, "SA_Score")
if _SA_SCORE_DIR not in sys.path:
    sys.path.append(_SA_SCORE_DIR)
import sascorer


# ---------------------------------------------------------------------------
# SmallWorld API
# ---------------------------------------------------------------------------

SW = "https://sw.docking.org/search/view"

# Default ZINC collection queried by every novelty calculation.
# Another option: "Zinc-All-25Q2-1.6B". Call sw_databases() to list them all.
DEFAULT_DB = "ZINC20-All-25Q2-1.9B"

# Maximum number of SmallWorld requests in flight at once. Kept low so the
# public server is not flooded.
SW_MAX_CONCURRENT = 5

# A failed lookup is retried this many times when the failure looks transient
# (dropped connection, timeout, 429 or 5xx). The wait grows with each attempt:
# SW_RETRY_DELAY seconds before the first retry, twice that before the second.
SW_RETRIES = 2
SW_RETRY_DELAY = 5.0
SW_RETRY_STATUS = {429, 500, 502, 503, 504}

# One requests.Session per worker thread: Session is not guaranteed to be
# thread-safe, but a per-thread one still reuses its HTTPS connection.
_thread_local = threading.local()


def _session():
    """Return this thread's requests.Session, creating it on first use."""
    if not hasattr(_thread_local, "session"):
        _thread_local.session = requests.Session()
    return _thread_local.session


# ---------------------------------------------------------------------------
# Curation
# ---------------------------------------------------------------------------

# RDKit's MolStandardize is a C++ port of MolVS, and the tautomer step is the
# reason for using it: the MolVS original is pure Python and spends ~40 ms per
# molecule there, about 6x what this costs. The two do not always pick the same
# tautomer, so curated structures differ from MolVS output for a minority of
# molecules.
LARGEST_FRAGMENT = rdMolStandardize.LargestFragmentChooser()
UNCHARGER = rdMolStandardize.Uncharger()
REIONIZER = rdMolStandardize.Reionizer()
TAUTOMER_CANONICALIZER = rdMolStandardize.TautomerEnumerator()

ALLOWED_ELEMENTS = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Se", "Cl", "Br", "I"}

# Rejection codes returned by pretreatment in place of a curated SMILES.
CURATION_ERRORS = ("Error 1", "Error 2", "Error 3")

# How often pretreatment_batch calls its progress callback. Curation runs at
# roughly 40 ms per molecule, so this reports about every half second.
PROGRESS_EVERY = 10


# ---------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------

# ECFP4 is Morgan with radius 2; the generator is stateless, so one is enough.
MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

# The Ertl & Schuffenhauer score is defined on a 1-10 scale.
SA_MIN, SA_MAX = 1.0, 10.0


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

# Tanimoto similarity above which two molecules fall in the same cluster.
# Set below the 0.3 the BitBIRCH authors recommend for ECFP4, which gave
# better-balanced clusters in practice. Lowering it gives broader clusters,
# and since one SmallWorld request is spent per cluster, it is the knob that
# sets how long a large run takes.
DEFAULT_CLUSTER_THRESHOLD = 0.2

# Upper bound of the similarity threshold offered in the app; above it the
# clusters get so tight that clustering barely cuts the request count.
MAX_CLUSTER_THRESHOLD = 0.5

# Parameters for the BitBIRCH clustering algorithm.
BITBIRCH_BRANCHING_FACTOR = 100
BITBIRCH_MERGE_CRITERION = "diameter"

# Above this many molecules the SmallWorld API calls dominate the runtime, so
# clustering first is strongly recommended. Used by the app to warn the user.
LARGE_DATASET_SIZE = 500

# ---------------------------------------------------------------------------
# Priority plot
# ---------------------------------------------------------------------------

# Fallback color when the user does not map an activity column.
DEFAULT_POINT_COLOR = "#2b7bba"

# The three axes produced by calculate_axis, all on a 0-1 higher-is-better scale.
DEFAULT_AXES = ("Novelty_Score", "Synthesizability", "QED")

# Continuous scale for a numeric activity column.
DEFAULT_COLORSCALE = "Viridis"

# Qualitative colors for a categorical activity column; they wrap around when
# there are more categories than colors.
CATEGORY_PALETTE = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Default lower bound of the desirability zone on each axis; the zone runs
# from this value up to 1 on all three axes.
DEFAULT_DESIRABILITY = 0.7
DESIRABILITY_COL = "In_Desirability_Zone"
DESIRABILITY_COLOR = "#d62728"

# Passed to fig.show(): the camera button then downloads a vector SVG.
PLOT_CONFIG = {
    "displaylogo": False,
    "toImageButtonOptions": {"format": "svg", "filename": "denomap_plot"},
}


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------

def _as_smiles_list(smiles):
    """Normalize any supported input into (list_of_smiles, was_scalar, index).

    A single SMILES is a str, and a str is iterable: without this check,
    list("CCO") would silently score three separate atoms instead of one
    molecule. Series keep their index so the result can be assigned back to
    the DataFrame it came from.
    """
    if isinstance(smiles, str) or smiles is None:
        return [smiles], True, None
    if isinstance(smiles, pd.Series):
        return smiles.tolist(), False, smiles.index
    if hasattr(smiles, "__iter__"):
        return list(smiles), False, None
    return [smiles], True, None  # NaN and other non-iterable scalars


# ---------------------------------------------------------------------------
# Curation
# ---------------------------------------------------------------------------

def pretreatment(smi):
    """Standardize one SMILES, or return the code that rejects it.

    Error 1: unparsable. Error 2: contains an element outside ALLOWED_ELEMENTS.
    Error 3: the standardization pipeline itself raised.
    """
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return "Error 1"
        mol = rdMolStandardize.Cleanup(mol)
        mol = LARGEST_FRAGMENT.choose(mol)
        if {atom.GetSymbol() for atom in mol.GetAtoms()} - ALLOWED_ELEMENTS:
            return "Error 2"
        mol = UNCHARGER.uncharge(mol)
        mol = REIONIZER.reionize(mol)
        mol = TAUTOMER_CANONICALIZER.Canonicalize(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return "Error 3"


def pretreatment_batch(smiles, progress=None):
    """Run pretreatment over many SMILES, reporting how far along it is.

    Curation is the slow local step - tautomer canonicalization alone costs
    about 40 ms per molecule - so a few thousand rows take minutes with no
    output of their own. The callback is what keeps the caller able to show
    that it is still working.

    progress : optional callable(done, total).
    """
    smiles_list, _, _ = _as_smiles_list(smiles)
    total = len(smiles_list)

    results = []
    for done, smi in enumerate(smiles_list, start=1):
        results.append(pretreatment(smi))
        if progress is not None and (done % PROGRESS_EVERY == 0 or done == total):
            progress(done, total)
    return results


def _inchikey(smiles):
    """InChIKey of an already curated SMILES, used only to spot duplicates."""
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


def curate_smiles(df, smiles_col="SMILES", drop_duplicates=True, progress=None):
    """Standardize a DataFrame of molecules and report validity / uniqueness.

    Returns (curated_df, stats). The curated structures replace the input
    column under the name SMILES, and an InChIKey column is added. Rejected
    molecules are dropped, so the result can be shorter than df.
    """
    df = df.copy()
    df["SMILES_Curated"] = pretreatment_batch(df[smiles_col], progress=progress)
    n_total = len(df)

    valid = df[~df["SMILES_Curated"].isin(CURATION_ERRORS)].reset_index(drop=True)
    n_valid = len(valid)
    valid["InChIKey"] = valid["SMILES_Curated"].apply(_inchikey)
    n_unique = valid["InChIKey"].nunique()

    validity = n_valid / n_total if n_total else 0
    uniqueness = n_unique / n_valid if n_valid else 0

    if drop_duplicates:
        valid = valid.drop_duplicates(subset="InChIKey", keep="first").reset_index(drop=True)

    valid = valid.drop(columns=[smiles_col]).rename(columns={"SMILES_Curated": "SMILES"})
    return valid, {"validity": validity, "uniqueness": uniqueness, "n_total": n_total,
                   "n_valid": n_valid, "n_unique": n_unique}


# ---------------------------------------------------------------------------
# BitBIRCH clustering
# ---------------------------------------------------------------------------

def bitbirch_clusters(smiles, similarity_threshold=DEFAULT_CLUSTER_THRESHOLD,
                      progress=None):
    """Group molecules with BitBIRCH on the ECFP4 Tanimoto similarity.

    Returns a list of clusters, each a tuple of positional indices into the
    input whose first element is the cluster representative. Molecules that
    cannot be parsed are left out of every cluster.

    BitBIRCH inserts each molecule into a BIRCH tree once rather than comparing
    every pair, so it stays linear where an all-pairs method does not - the
    reason it is used here instead of RDKit's Butina.

    The representative is each cluster's medoid: the real input molecule
    closest to the cluster's centroid, so it can be scored like any other row.

    Requires the bblean package.
    """
    if not 0.0 <= similarity_threshold <= 1.0:
        raise ValueError("similarity_threshold must be between 0 and 1.")

    try:
        import bblean
    except ImportError as exc:  # optional dependency, so say what is missing
        raise ImportError(
            "BitBIRCH clustering needs the 'bblean' package: pip install bblean"
        ) from exc

    smiles_list, _, _ = _as_smiles_list(smiles)

    # fps_from_smiles can skip unparsable entries, but then the indices it
    # returns no longer address the caller's rows. Filter here instead, so
    # the returned clusters still point at the input positions.
    positions, valid = [], []
    for position, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) and smi.strip() else None
        if mol is None:
            continue
        positions.append(position)
        valid.append(smi)

    if not valid:
        return []

    if progress is not None:
        progress(0, 2)
    fingerprints = bblean.fps_from_smiles(valid, kind="ecfp4", n_features=2048,
                                          pack=True)
    if progress is not None:
        progress(1, 2)

    tree = bblean.BitBirch(threshold=similarity_threshold,
                           branching_factor=BITBIRCH_BRANCHING_FACTOR,
                           merge_criterion=BITBIRCH_MERGE_CRITERION)
    tree.fit(fingerprints)

    members = tree.get_cluster_mol_ids()
    # medoid_idxs index into each cluster's own member list, not into the
    # fingerprint array, so they have to be resolved through it.
    medoids = tree.get_medoids_mol_ids(fingerprints)["medoid_idxs"]
    if progress is not None:
        progress(2, 2)

    clusters = []
    for cluster, medoid in zip(members, medoids):
        cluster = [int(i) for i in cluster]
        seat = int(medoid)
        cluster.insert(0, cluster.pop(seat))  # representative goes first
        clusters.append(tuple(positions[i] for i in cluster))
    return clusters


def cluster_centroids(df, smiles_col="SMILES",
                      similarity_threshold=DEFAULT_CLUSTER_THRESHOLD,
                      progress=None):
    """Reduce a DataFrame to one representative molecule per cluster.

    Returns (centroids_df, stats). Only the centroid row of each cluster is
    kept, with Cluster_ID and Cluster_Size columns appended, so the result can
    be scored in place of the full set - one SmallWorld request per cluster
    instead of one per molecule.
    """
    clusters = bitbirch_clusters(df[smiles_col],
                                 similarity_threshold=similarity_threshold,
                                 progress=progress)

    centroid_positions = [cluster[0] for cluster in clusters]
    out = df.iloc[centroid_positions].copy()
    out["Cluster_ID"] = range(len(clusters))
    out["Cluster_Size"] = [len(cluster) for cluster in clusters]
    out = out.reset_index(drop=True)

    n_clustered = sum(len(cluster) for cluster in clusters)
    return out, {"n_total": len(df), "n_clustered": n_clustered,
                 "n_clusters": len(clusters),
                 "similarity_threshold": similarity_threshold}


# ---------------------------------------------------------------------------
# SmallWorld API - nearest neighbor lookup in ZINC
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def sw_databases():
    """List the collections the SmallWorld server currently exposes."""
    r = requests.get("https://sw.docking.org/search/maps", timeout=60)
    r.raise_for_status()
    return {v["name"]: v for v in r.json().values() if v.get("enabled")}


def _is_transient(exc):
    """True for failures worth retrying: network blips and server overload."""
    if isinstance(exc, (requests.ConnectionError, requests.Timeout)):
        return True
    response = getattr(exc, "response", None)
    return response is not None and response.status_code in SW_RETRY_STATUS


def sw_nearest_neighbor(smiles, db=None, dist="0-16", timeout=300):
    """Return (zinc_id, smiles) of the closest ZINC molecule, or None.

    Transient failures are retried SW_RETRIES times before giving up.
    """
    params = {"smi": smiles, "db": db or DEFAULT_DB, "fmt": "tsv",
              "length": 1, "top": 1, "dist": dist}
    for attempt in range(SW_RETRIES + 1):
        try:
            r = _session().get(SW, params=params, timeout=timeout)
            r.raise_for_status()
            break
        except requests.RequestException as exc:
            if attempt == SW_RETRIES or not _is_transient(exc):
                raise
            time.sleep(SW_RETRY_DELAY * (attempt + 1))
    lines = r.text.strip().split("\n")[1:]
    if not lines:
        return None
    parts = lines[0].split("\t")[0].split()
    return (parts[1], parts[0]) if len(parts) >= 2 else None


def _fetch_neighbor(smiles, db, dist, timeout):
    """Look one SMILES up in ZINC, turning any failure into an error record.

    Network-bound, so it runs in a small thread pool (SW_MAX_CONCURRENT
    requests at a time); threads release the GIL while waiting on the network.
    """
    try:
        neighbor = sw_nearest_neighbor(smiles, db=db, dist=dist, timeout=timeout)
    except Exception as exc:
        return {"smiles": smiles, "neighbor_zinc_id": None, "neighbor_smiles": None, "error": str(exc)}
    if neighbor is None:
        return {"smiles": smiles, "neighbor_zinc_id": None, "neighbor_smiles": None, "error": "No neighbor found"}
    zinc_id, ref_smiles = neighbor
    return {"smiles": smiles, "neighbor_zinc_id": zinc_id, "neighbor_smiles": ref_smiles, "error": None}


# ---------------------------------------------------------------------------
# Fingerprints and ZINC Novelty Score
# ---------------------------------------------------------------------------

def _tanimoto(smi_query, smi_ref, fingerprint_fn):
    """Compute the Tanimoto coefficient between two SMILES using the given fingerprint function."""
    mol_query = Chem.MolFromSmiles(smi_query)
    mol_ref = Chem.MolFromSmiles(smi_ref)
    if mol_query is None or mol_ref is None:
        raise ValueError("Could not parse one of the SMILES strings.")
    return DataStructs.TanimotoSimilarity(fingerprint_fn(mol_query), fingerprint_fn(mol_ref))


def tc_ecfp4(smi_query, smi_ref):
    """Tanimoto coefficient using ECFP4 (Morgan, radius=2) fingerprints."""
    return _tanimoto(smi_query, smi_ref, MORGAN_GENERATOR.GetFingerprint)


def tc_path(smi_query, smi_ref):
    """Tanimoto coefficient using RDKit's topological (path-based) fingerprint."""
    return _tanimoto(smi_query, smi_ref, Chem.RDKFingerprint)


def novelty_score(smi_query, smi_ref):
    """ZINC Novelty Score (ZNS), from 0 (identical) to 1 (fully novel).

    ZNS = 1.0 - (Tc_ecfp4 + Tc_path) / 2. Compares the query molecule against
    its nearest ZINC reference molecule using both ECFP4 and RDKit path-based
    fingerprints; higher values mean the query is more novel relative to the
    ZINC reference. Both Tanimoto coefficients are already in 0-1, so the
    score needs no further rescaling.
    """
    return 1.0 - (tc_ecfp4(smi_query, smi_ref) + tc_path(smi_query, smi_ref)) / 2


def _score_pair(record):
    """Compute the fingerprints and ZNS metric for one query/reference pair.

    Runs in a worker thread. RDKit's C++ fingerprint/Tanimoto routines release
    the GIL for most of their work, so threads still parallelize this CPU step.
    """
    if record["error"] is not None:
        return {**record, "novelty_score": None}
    try:
        return {**record, "novelty_score": novelty_score(record["smiles"], record["neighbor_smiles"])}
    except Exception as exc:
        return {**record, "novelty_score": None, "error": str(exc)}


def novelty_score_batch(smiles, db=None, dist="0-16", timeout=300, max_workers=None,
                        progress=None):
    """Score a batch of SMILES against their nearest ZINC neighbors.

    Accepts a single SMILES, a list, or a pandas Series / DataFrame column, and
    returns a DataFrame with one row per input SMILES, in input order.

    Runs in two stages:
      1. Concurrent SmallWorld API requests (at most SW_MAX_CONCURRENT at a
         time), one nearest-neighbor lookup per SMILES.
      2. Once all lookups are done, the descriptor (ECFP4 / path fingerprint) and
         novelty score calculations run in parallel across threads.

    progress : optional callable(done, total) invoked after each API lookup, so
        a caller (e.g. the Streamlit app) can render a progress bar.
    """
    smiles_list, _, _ = _as_smiles_list(smiles)

    # Stage 1: concurrent API requests. Results are slotted back by index so
    # the rows keep input order; progress is reported from this (the calling)
    # thread, which Streamlit requires.
    total = len(smiles_list)
    records = [None] * total
    with ThreadPoolExecutor(max_workers=SW_MAX_CONCURRENT) as executor:
        futures = {executor.submit(_fetch_neighbor, smi, db, dist, timeout): i
                   for i, smi in enumerate(smiles_list)}
        for done, future in enumerate(as_completed(futures), start=1):
            records[futures[future]] = future.result()
            if progress is not None:
                progress(done, total)

    # Stage 2: parallel descriptor + ZNS calculation. executor.map yields in
    # input order, so the rows still line up with the input SMILES.
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_score_pair, records))

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# QED and synthetic accessibility
# ---------------------------------------------------------------------------

def _score_smiles(smiles, score_fn):
    """Apply a per-molecule scoring function over a scalar or a batch of SMILES.

    Unparsable or missing SMILES score as None instead of raising, so one bad
    row cannot abort a whole library. Repeated SMILES are computed once.
    """
    smiles_list, was_scalar, index = _as_smiles_list(smiles)

    cache = {}
    scores = []
    for smi in smiles_list:
        if smi not in cache:
            mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) and smi.strip() else None
            try:
                cache[smi] = None if mol is None else float(score_fn(mol))
            except Exception:
                cache[smi] = None
        scores.append(cache[smi])

    if was_scalar:
        return scores[0]
    if index is not None:
        return pd.Series(scores, index=index, dtype="float64")
    return scores


def qed_score(smiles):
    """Quantitative Estimate of Drug-likeness, from 0 (poor) to 1 (drug-like).

    Accepts a single SMILES, a list, or a pandas Series / DataFrame column,
    and mirrors the input shape: a float for a single SMILES, a Series for a
    Series, a list otherwise.
    """
    return _score_smiles(smiles, QED.qed)


def sa_score(smiles):
    """Ertl & Schuffenhauer synthetic accessibility, from 1 (easy) to 10 (hard).

    See synthesizability for the 0-1 version used as a priority-plot axis.
    Same input/output contract as qed_score.
    """
    return _score_smiles(smiles, sascorer.calculateScore)


def _synthesizability(mol):
    """SA score of one molecule, rescaled to 0-1 and inverted."""
    scaled = (SA_MAX - sascorer.calculateScore(mol)) / (SA_MAX - SA_MIN)
    return min(1.0, max(0.0, scaled))  # sascorer can drift slightly outside 1-10


def synthesizability(smiles):
    """SA score rescaled to 0-1 and inverted: higher means easier to synthesize.

    Uses (10 - SA) / 9, which puts synthetic accessibility on the same 0-1,
    higher-is-better scale as QED, so both axes of the priority plot read in the
    same direction. Same input/output contract as qed_score.
    """
    return _score_smiles(smiles, _synthesizability)


# ---------------------------------------------------------------------------
# Priority plot axes
# ---------------------------------------------------------------------------

def calculate_axis(df, smiles_col="SMILES", curate=True, drop_duplicates=True,
                   cluster=False, cluster_threshold=DEFAULT_CLUSTER_THRESHOLD,
                   db=None, dist="0-16", timeout=300,
                   max_workers=None, keep_neighbor=True, progress=None,
                   curate_progress=None, cluster_progress=None):
    """Compute the three priority-plot axes for a DataFrame of molecules.

    Returns a copy of df with the Novelty_Score, Synthesizability and QED
    columns appended, all three on a 0-1 scale where higher is better. The
    novelty score is the slow step, since it needs one SmallWorld API call per
    molecule; QED and synthesizability are local RDKit calculations.

    With curate=True (the default) the molecules go through curate_smiles
    first, so all three axes describe the standardized structure rather than
    the raw input.

    With cluster=True the set is reduced to one centroid per cluster
    before scoring, which cuts the number of SmallWorld requests down to the
    number of clusters. cluster_threshold is the Tanimoto similarity above
    which two molecules share a cluster.

    progress : optional callable(done, total) forwarded to novelty_score_batch.
    curate_progress : optional callable(done, total) for the curation step,
        which is the slow local one - roughly 40 ms per molecule.
    cluster_progress : optional callable(done, total) for the clustering
        pass, which is quadratic and otherwise reports nothing at all.
    """
    out = df.copy()

    if curate:
        out, _ = curate_smiles(out, smiles_col=smiles_col, drop_duplicates=drop_duplicates,
                               progress=curate_progress)
        smiles_col = "SMILES"

    if cluster:
        out, _ = cluster_centroids(out, smiles_col=smiles_col,
                                   similarity_threshold=cluster_threshold,
                                   progress=cluster_progress)

    novelty = novelty_score_batch(out[smiles_col], db=db, dist=dist,
                                  timeout=timeout, max_workers=max_workers,
                                  progress=progress)

    # novelty_score_batch returns a fresh RangeIndex, so align by position
    # (row order is preserved) rather than by index label.
    if keep_neighbor:
        out["Neighbor_ZINC_ID"] = novelty["neighbor_zinc_id"].values
        out["Neighbor_SMILES"] = novelty["neighbor_smiles"].values
    out["Novelty_Score"] = novelty["novelty_score"].values

    out["Synthesizability"] = synthesizability(out[smiles_col])
    out["QED"] = qed_score(out[smiles_col])

    return out


# ---------------------------------------------------------------------------
# Molecule depictions
# ---------------------------------------------------------------------------

def _wrap(text, width=38):
    """Break a long SMILES so the tooltip does not grow past the plot."""
    text = str(text)
    return "<br>".join(text[i:i + width] for i in range(0, len(text), width)) or "-"


def _camera_eye(azimuth, elevation, distance=1.9):
    """Turn azimuth / elevation in degrees into a Plotly scene camera position."""
    a, e = math.radians(azimuth), math.radians(elevation)
    return dict(x=distance * math.cos(e) * math.cos(a),
                y=distance * math.cos(e) * math.sin(a),
                z=distance * math.sin(e))


def _axis_style(title):
    """Shared styling for the three scene axes."""
    return dict(title=dict(text=title, font=dict(size=12)),
                backgroundcolor="#fbfcfd", gridcolor="#dfe4ea",
                zerolinecolor="#c8d0d9", showbackground=True,
                tickfont=dict(size=10, color="#52616f"))


# ---------------------------------------------------------------------------
# Desirability zone
# ---------------------------------------------------------------------------

def _as_thresholds(thresholds):
    """Accept one number for all axes or one per axis; return three floats."""
    if np.isscalar(thresholds):
        return (float(thresholds),) * 3
    thresholds = tuple(float(t) for t in thresholds)
    if len(thresholds) != 3:
        raise ValueError("Pass one desirability threshold or exactly three.")
    return thresholds


def desirability_mask(df, thresholds=DEFAULT_DESIRABILITY, axes=DEFAULT_AXES):
    """Boolean Series: True where every axis is at or above its threshold.

    Molecules missing any axis value (e.g. a failed novelty lookup) are False.
    """
    mask = pd.Series(True, index=df.index)
    for col, t in zip(axes, _as_thresholds(thresholds)):
        mask &= df[col].ge(t).fillna(False)
    return mask


def _desirability_traces(thresholds):
    """Wireframe box from the thresholds up to 1 on every axis.

    Only the edges are drawn: in a 3D scene every surface, even one with
    hoverinfo="skip", is written to the WebGL pick buffer and blocks the hover
    of the points behind it, so a shaded Mesh3d would hide the molecules
    inside the zone from the tooltip and the structure panel.
    """
    (x0, y0, z0), (x1, y1, z1) = _as_thresholds(thresholds), (1.0, 1.0, 1.0)
    xs = [x0, x1, x1, x0, x0, x1, x1, x0]
    ys = [y0, y0, y1, y1, y0, y0, y1, y1]
    zs = [z0, z0, z0, z0, z1, z1, z1, z1]
    # Trace the 12 edges as one polyline; None breaks the line between them.
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
    ex, ey, ez = [], [], []
    for a, b in edges:
        ex += [xs[a], xs[b], None]
        ey += [ys[a], ys[b], None]
        ez += [zs[a], zs[b], None]
    return [go.Scatter3d(x=ex, y=ey, z=ez, mode="lines", hoverinfo="skip",
                         line=dict(color=DESIRABILITY_COLOR, width=4),
                         name="Desirability zone", showlegend=True,
                         legend="legend2")]


# ---------------------------------------------------------------------------
# De novo priority plot
# ---------------------------------------------------------------------------

def denovo_plot(df, smiles_col, activity_col=None, axes=DEFAULT_AXES,
                id_col=None, palette=None, point_size=5, width=820, height=680,
                azimuth=35.0, elevation=22.0,
                title="DeNoMAP", structures=True, axis_kwargs=None,
                desirability=None):
    """Interactive 3D scatter of the three priority-plot axes.

    Every molecule is one point in the Novelty / Synthesizability / QED space.
    Pass the result to plot_with_panel_html() to render it next to a panel that
    draws whichever molecule is hovered.

    If the default axis columns are not in df yet, they are computed on the fly
    by calling calculate_axis(df, smiles_col=smiles_col). That step is the slow
    one - one SmallWorld API request per molecule.

    activity_col : optional column driving the point color. Numeric columns get
        a continuous scale plus a color bar; anything else is treated as
        categories, one trace each, and gets a legend.
    desirability : optional threshold (one for all axes, or one per axis). Draws
        a wireframe box from the thresholds up to 1 on every axis and fixes
        the axis ranges to 0-1 so the box keeps its true proportions.
    """
    x_col, y_col, z_col = axes
    if smiles_col not in df.columns:
        raise ValueError(f"Column not found in the DataFrame: {smiles_col!r}")

    # Only the standard axes can be reconstructed; a custom set of axis names
    # says the values come from somewhere calculate_axis knows nothing about.
    missing_axes = [c for c in axes if c not in df.columns]
    if missing_axes and tuple(axes) == DEFAULT_AXES:
        df = calculate_axis(df, smiles_col=smiles_col, **(axis_kwargs or {}))
        if smiles_col not in df.columns and "SMILES" in df.columns:
            smiles_col = "SMILES"

    needed = list(axes) + ([activity_col] if activity_col else []) + ([id_col] if id_col else [])
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Columns not found in the DataFrame: {missing}")

    # Molecules without coordinates cannot be placed in the scene at all; the
    # usual cause is a novelty lookup that failed against the SmallWorld API.
    data = df.dropna(subset=list(axes)).reset_index(drop=True)
    if data.empty:
        raise ValueError("No molecule has all three axis values.")

    labels = data[id_col].astype(str) if id_col else data.index.astype(str)
    smiles = data[smiles_col].astype(str)

    # The panel draws the hovered molecule in the browser, so the point only
    # has to carry its SMILES.
    depiction_smiles = smiles if structures else pd.Series([""] * len(data))

    # customdata carries what the tooltip prints, plus the SMILES that
    # the panel draws on hover.
    activity_text = (data[activity_col].astype(str) if activity_col
                     else pd.Series([""] * len(data)))
    customdata = np.column_stack([labels, smiles.map(_wrap), activity_text,
                                  depiction_smiles])

    activity_row = f"{activity_col}: %{{customdata[2]}}<br>" if activity_col else ""
    hovertemplate = (
        "<b>%{customdata[0]}</b><br>"
        f"{x_col}: %{{x:.3f}}<br>"
        f"{y_col}: %{{y:.3f}}<br>"
        f"{z_col}: %{{z:.3f}}<br>"
        f"{activity_row}"
        "<span style='color:#8a94a0'>%{customdata[1]}</span>"
        "<extra></extra>"
    )

    marker = dict(size=point_size, opacity=0.85,
                  line=dict(width=0.5, color="#33404d"))
    scatter_kwargs = dict(mode="markers", hovertemplate=hovertemplate,
                          hoverlabel=dict(align="left", bgcolor="white",
                                          bordercolor="#c8d0d9"))

    traces = []
    if activity_col is None:
        traces.append(go.Scatter3d(
            x=data[x_col], y=data[y_col], z=data[z_col], customdata=customdata,
            marker={**marker, "color": DEFAULT_POINT_COLOR},
            showlegend=False, **scatter_kwargs))
    elif pd.api.types.is_numeric_dtype(data[activity_col]):
        traces.append(go.Scatter3d(
            x=data[x_col], y=data[y_col], z=data[z_col], customdata=customdata,
            marker={**marker,
                    "color": data[activity_col].astype(float),
                    "colorscale": palette or DEFAULT_COLORSCALE,
                    "colorbar": dict(title=dict(text=activity_col, side="right"),
                                     thickness=14, len=0.6, x=1.02)},
            showlegend=False, **scatter_kwargs))
    else:
        # One trace per category, which is what gives Plotly a clickable legend.
        categories = data[activity_col].astype(str)
        factors = sorted(categories.unique())
        colors = list(palette) if palette else list(CATEGORY_PALETTE)
        for i, factor in enumerate(factors):
            mask = (categories == factor).values
            traces.append(go.Scatter3d(
                x=data.loc[mask, x_col], y=data.loc[mask, y_col],
                z=data.loc[mask, z_col], customdata=customdata[mask],
                name=factor, marker={**marker, "color": colors[i % len(colors)]},
                **scatter_kwargs))

    scene_axes = dict(xaxis=_axis_style(x_col), yaxis=_axis_style(y_col),
                      zaxis=_axis_style(z_col))
    if desirability is not None:
        traces.extend(_desirability_traces(desirability))
        # Pin every axis to the same 0-1 span; with autoscaled axes each one
        # stretches differently and the box no longer looks like a box.
        for ax in scene_axes.values():
            ax["range"] = [0, 1]

    fig = go.Figure(traces)
    fig.update_layout(
        title=dict(text=title, x=0.02, xanchor="left", font=dict(size=15)),
        width=width, height=height, template="plotly_white",
        margin=dict(l=0, r=0, t=48, b=36 if desirability is not None else 0),
        legend=dict(title=dict(text=activity_col or ""), itemsizing="constant",
                    yanchor="top", y=0.95, xanchor="left", x=0.98),
        # Separate legend for the zone, centered under the scene, so it stays
        # apart from the activity categories at the top right.
        legend2=dict(orientation="h", xanchor="center", x=0.5,
                     yanchor="top", y=0),
        scene=dict(**scene_axes, aspectmode="cube",
                   camera=dict(eye=_camera_eye(azimuth, elevation))),
    )
    return fig


# JS twin of the old Bokeh hover: on plotly_hover, copy the depiction that
# denovo_plot stored in customdata[3] into the panel next to the plot.
# Draws the hovered molecule client-side, so the figure only ever carries
# SMILES strings. Kept as a classic script tag (not a module) and placed ahead
# of the plot so SmilesDrawer is defined by the time this runs.
SMILES_DRAWER_CDN = "https://cdn.jsdelivr.net/npm/smiles-drawer@2.1.7/dist/smiles-drawer.min.js"

HOVER_JS = """
(function() {
    var plot = document.getElementById("__PLOT_ID__");
    if (!plot || !plot.on) { return; }

    // SmiDrawer draws into an <svg>; the plain Drawer class targets a canvas
    // through SvgWrapper and throws on one. One drawer serves every hover.
    var drawer = null;
    function getDrawer() {
        if (!drawer && window.SmilesDrawer && window.SmilesDrawer.SmiDrawer) {
            drawer = new SmilesDrawer.SmiDrawer({width: __IMG_W__, height: __IMG_H__});
        }
        return drawer;
    }

    function note(panel, text) {
        var el = document.createElement("div");
        el.style.cssText = "color:#b04a4a; font-size:11px;";
        el.textContent = text;
        panel.appendChild(el);
    }

    plot.on("plotly_hover", function(event) {
        var panel = document.getElementById("__PANEL_ID__");
        if (!panel) { return; }
        var point = event.points[0] || {};
        var row = point.customdata;
        if (!row && point.data && point.data.customdata) {
            row = point.data.customdata[point.pointNumber];
        }
        panel.textContent = "";
        panel.style.display = "block";

        var caption = document.createElement("div");
        caption.style.cssText = "font-weight:600; margin-bottom:4px; color:#3e4c59;";
        caption.textContent = row ? row[0] : "no data on this point";
        panel.appendChild(caption);
        if (!row || !row[3]) { return; }

        var d = getDrawer();
        if (!d) { note(panel, "Structure renderer did not load."); return; }

        // createElementNS, not createElement: an HTML <svg> element gets no SVG
        // behaviour and SmiDrawer silently draws nothing into it.
        var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("width", __IMG_W__);
        svg.setAttribute("height", __IMG_H__);
        svg.style.cssText = "display:block; background:#ffffff; border:1px solid #dfe4ea; border-radius:6px;";
        panel.appendChild(svg);

        try {
            d.draw(row[3], svg, "light", null, function(err) {
                svg.remove();
                note(panel, "Could not draw this structure.");
            });
        } catch (err) {
            svg.remove();
            note(panel, "Could not draw this structure.");
        }
    });
})();
"""


def plot_with_panel_html(fig, image_size=(230, 190), include_plotlyjs="cdn"):
    """Build the HTML of the figure and the structure panel side by side.

    Returns a self-contained HTML string suitable for embedding with
    streamlit.components.v1.html. include_plotlyjs="cdn" keeps it small but
    needs a connection; "inline" embeds plotly.js so it also works offline.
    """
    plot_id = "denovo-plot-" + uuid.uuid4().hex[:8]
    panel_id = plot_id + "-panel"
    img_width, img_height = image_size

    plot_html = fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs,
                            div_id=plot_id, config=PLOT_CONFIG,
                            post_script=HOVER_JS.replace("__PLOT_ID__", plot_id)
                                                .replace("__PANEL_ID__", panel_id)
                                                .replace("__IMG_W__", str(img_width))
                                                .replace("__IMG_H__", str(img_height)))
    panel_html = (
        f"<div id='{panel_id}' style='width:{img_width}px; min-height:{img_height}px;"
        " background:#ffffff;"
        " display:flex; align-items:center; justify-content:center; color:#8a94a0;"
        " font:12px Helvetica, Arial, sans-serif; border:1px dashed #dfe4ea;"
        " border-radius:6px; padding:8px;'>Hover a point</div>"
    )
    return (f"<script src='{SMILES_DRAWER_CDN}'></script>"
            "<div style='display:flex; align-items:center; gap:16px; flex-wrap:wrap;'>"
            f"<div>{plot_html}</div>{panel_html}</div>")
