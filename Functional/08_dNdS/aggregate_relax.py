#!/usr/bin/env python3
"""
aggregate_relax.py

Walk a directory of HyPhy RELAX output JSONs, pull out the key test
statistics for each gene/sweep, and apply Benjamini-Hochberg FDR
correction for multiple testing.

Usage:
    python3 aggregate_relax.py --results-dir relax_results --out relax_summary.csv

Expected input layout (adjust parsing in parse_relax_json() if yours differs):
    relax_results/
        sweep_7/
            OHFJIJCL_00001.RELAX.json
            OHFJIJCL_00002.RELAX.json
            ...
        sweep_12/
            ...

Only needs numpy + pandas (BH-FDR is implemented directly, no statsmodels
dependency).
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

# HyPhy writes bare `inf` / `-inf` / `nan` tokens (e.g. in "Evidence Ratios"
# for branches with an infinite likelihood ratio). These are not valid JSON --
# Python's json module accepts "Infinity"/"-Infinity"/"NaN" as an extension,
# but not the bare C-style tokens HyPhy emits -- so json.load chokes on them
# with a generic "Expecting value" error. Patch them to the JSON-extension
# spelling before parsing. The colon-anchored pattern only matches numeric
# value position, so it won't touch "inf"/"nan" occurring inside strings.
_INF_NAN_RE = [
    (re.compile(r'(:\s*)-inf\b'), r'\1-Infinity'),
    (re.compile(r'(:\s*)inf\b'), r'\1Infinity'),
    (re.compile(r'(:\s*)nan\b', re.IGNORECASE), r'\1NaN'),
]


def load_hyphy_json(path: Path) -> dict:
    with open(path) as f:
        text = f.read()
    for pattern, repl in _INF_NAN_RE:
        text = pattern.sub(repl, text)
    return json.loads(text)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. NaNs are ignored and returned as NaN."""
    pvals = np.asarray(pvals, dtype=float)
    q = np.full_like(pvals, np.nan)
    valid = ~np.isnan(pvals)
    p = pvals[valid]
    n = len(p)
    if n == 0:
        return q
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q_raw = p * n / ranks
    q_sorted = q_raw[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]  # enforce monotonicity
    q_valid = np.empty(n)
    q_valid[order] = np.clip(q_sorted, 0, 1)
    q[valid] = q_valid
    return q


def parse_relax_json(path: Path) -> dict:
    """Extract the fields we care about from one RELAX.json. Returns a dict
    with NaNs for anything missing, and a 'status' column flagging problems
    so failed/incomplete runs don't silently vanish from the table."""
    sweep = path.parent.name
    gene = path.name.replace(".RELAX.json", "")

    row = {
        "sweep": sweep,
        "gene": gene,
        "file": str(path),
        "status": "ok",
        "p_value": np.nan,
        "K": np.nan,
        "LRT": np.nan,
        "n_sequences": np.nan,
        "n_sites": np.nan,
        "ll_alt": np.nan,
        "ll_null": np.nan,
        "aicc_alt": np.nan,
        "aicc_null": np.nan,
    }

    try:
        data = load_hyphy_json(path)
    except (json.JSONDecodeError, OSError) as e:
        row["status"] = f"failed_to_parse: {e}"
        return row

    test = data.get("test results")
    if test is None:
        row["status"] = "missing_test_results"
    else:
        row["p_value"] = test.get("p-value", np.nan)
        row["K"] = test.get("relaxation or intensification parameter", np.nan)
        row["LRT"] = test.get("LRT", np.nan)

    inp = data.get("input", {})
    row["n_sequences"] = inp.get("number of sequences", np.nan)
    row["n_sites"] = inp.get("number of sites", np.nan)

    fits = data.get("fits", {})
    alt = fits.get("RELAX alternative", {})
    null = fits.get("RELAX null", {})
    row["ll_alt"] = alt.get("Log Likelihood", np.nan)
    row["ll_null"] = null.get("Log Likelihood", np.nan)
    row["aicc_alt"] = alt.get("AIC-c", np.nan)
    row["aicc_null"] = null.get("AIC-c", np.nan)

    # Flag suspicious boundary estimates of K -- a common symptom of a poor fit
    # (e.g. too little data / too few substitutions on the test branches).
    if row["status"] == "ok" and not pd.isna(row["K"]):
        if row["K"] < 1e-3 or row["K"] > 100:
            row["status"] = "boundary_K_check_convergence"

    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, type=Path,
                     help="Top-level directory containing per-sweep subfolders of *.RELAX.json")
    ap.add_argument("--out", required=True, type=Path, help="Output CSV path")
    ap.add_argument("--fdr-scope", choices=["per-sweep", "global"], default="per-sweep",
                     help="Apply BH-FDR within each sweep separately (default) or pooled across all genes/sweeps")
    ap.add_argument("--alpha", type=float, default=0.05, help="Significance threshold on q-value")
    args = ap.parse_args()

    files = sorted(args.results_dir.rglob("*.RELAX.json"))
    if not files:
        raise SystemExit(f"No *.RELAX.json files found under {args.results_dir}")

    rows = [parse_relax_json(f) for f in files]
    df = pd.DataFrame(rows)

    n_flagged = (df["status"] != "ok").sum()
    print(f"Parsed {len(df)} files ({n_flagged} flagged: parse failures, missing results, or boundary K)")

    # --- Multiple testing correction ---
    if args.fdr_scope == "per-sweep":
        df["q_value"] = df.groupby("sweep")["p_value"].transform(bh_fdr)
    else:
        df["q_value"] = bh_fdr(df["p_value"].to_numpy())

    df["significant"] = df["q_value"] < args.alpha

    def direction(row):
        if not row["significant"] or pd.isna(row["K"]):
            return "not significant"
        return "intensified (K>1)" if row["K"] > 1 else "relaxed (K<1)"

    df["direction"] = df.apply(direction, axis=1)

    df = df.sort_values(["sweep", "q_value"])
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")

    print("\nSummary by sweep:")
    print(df.groupby(["sweep", "direction"]).size().unstack(fill_value=0))

    if n_flagged:
        print(f"\n{n_flagged} genes flagged as non-'ok' status -- check the 'status' "
              f"column before trusting their q-values.")


if __name__ == "__main__":
    main()