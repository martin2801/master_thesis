#!/usr/bin/env python3
"""
diagnose_relax_power.py

For each sweep's genome-wide RELAX JSON, quantify how much real evolutionary
signal actually sits on the Test (Foreground) branch set, on top of the
plain K / LRT / p-value result. This is the same check done by hand for
sweep_1: sum the estimated branch lengths (under the "RELAX alternative"
fit) across all Test-labeled branches, convert to an expected substitution
count over the whole alignment, and check how concentrated that signal is
on a single branch (a "stem-only" pattern -- one normal-length branch
leading into the clade, with all descendant tips at ~0 -- means the clade's
own genomes are essentially identical to each other, which caps how much
RELAX can ever say about selection on that lineage, independent of
alignment length).

Handles HyPhy's bare inf/-inf/nan tokens the same way aggregate_relax.py
does, and the large file sizes these whole-genome RELAX runs produce.

Usage:
    python3 diagnose_relax_power.py --results-dir relax_results_concat --out power_diagnostics.csv

Expected layout:
    relax_results_concat/
        sweep_1/concatenated_core.RELAX.json
        sweep_2/concatenated_core.RELAX.json
        ...
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd

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


def diagnose_one(path: Path) -> dict:
    sweep = path.parent.name
    row = {
        "sweep": sweep, "file": str(path), "status": "ok",
        "K": None, "LRT": None, "p_value": None,
        "n_sites": None, "n_test_branches": None, "n_reference_branches": None,
        "sum_test_branch_len": None, "mean_test_branch_len": None,
        "sum_reference_branch_len": None, "mean_reference_branch_len": None,
        "expected_subs_on_test_set": None,
        "top_branch_frac_of_test_signal": None,
        "flag": "",
    }

    try:
        data = load_hyphy_json(path)
    except (json.JSONDecodeError, OSError) as e:
        row["status"] = f"failed_to_parse: {e}"
        return row

    test = data.get("test results", {})
    row["K"] = test.get("relaxation or intensification parameter")
    row["LRT"] = test.get("LRT")
    row["p_value"] = test.get("p-value")
    row["n_sites"] = data.get("input", {}).get("number of sites")

    tested = data.get("tested", {}).get("0", {})
    battr = data.get("branch attributes", {}).get("0", {})

    test_branches = [b for b, v in tested.items() if v == "Test"]
    ref_branches = [b for b, v in tested.items() if v == "Reference"]

    test_lens = [battr.get(b, {}).get("RELAX alternative", 0) for b in test_branches]
    ref_lens = [battr.get(b, {}).get("RELAX alternative", 0) for b in ref_branches]

    row["n_test_branches"] = len(test_branches)
    row["n_reference_branches"] = len(ref_branches)
    row["sum_test_branch_len"] = sum(test_lens) if test_lens else 0.0
    row["mean_test_branch_len"] = (sum(test_lens) / len(test_lens)) if test_lens else 0.0
    row["sum_reference_branch_len"] = sum(ref_lens) if ref_lens else 0.0
    row["mean_reference_branch_len"] = (sum(ref_lens) / len(ref_lens)) if ref_lens else 0.0

    if row["n_sites"] and row["sum_test_branch_len"] is not None:
        row["expected_subs_on_test_set"] = row["sum_test_branch_len"] * row["n_sites"]

    if test_lens and sum(test_lens) > 0:
        row["top_branch_frac_of_test_signal"] = max(test_lens) / sum(test_lens)

    # Heuristic flags to draw your eye to sweeps like sweep_1
    flags = []
    if row["expected_subs_on_test_set"] is not None and row["expected_subs_on_test_set"] < 10:
        flags.append("LOW_SUBSTITUTION_COUNT(<10 expected on whole Test set)")
    if row["top_branch_frac_of_test_signal"] is not None and row["top_branch_frac_of_test_signal"] > 0.8:
        flags.append("SIGNAL_CONCENTRATED_ON_ONE_BRANCH(>80%)")
    if row["n_test_branches"] is not None and row["n_test_branches"] <= 5:
        flags.append("VERY_FEW_TEST_BRANCHES(<=5)")
    row["flag"] = "; ".join(flags)

    return row


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--results-dir", type=Path,
                      help="Directory containing sweep_*/concatenated_core.RELAX.json (batch mode)")
    src.add_argument("--file", type=Path,
                      help="A single RELAX.json to diagnose on its own")
    ap.add_argument("--pattern", default="*.RELAX.json", help="Only used with --results-dir")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if args.file:
        files = [args.file]
        if not files[0].exists():
            raise SystemExit(f"{args.file} not found")
    else:
        files = sorted(args.results_dir.rglob(args.pattern))
        if not files:
            raise SystemExit(f"No files matching {args.pattern} under {args.results_dir}")

    rows = [diagnose_one(f) for f in files]
    df = pd.DataFrame(rows)
    df = df.sort_values("sweep")
    df.to_csv(args.out, index=False)

    print(f"Diagnosed {len(df)} sweep result(s). Wrote {args.out}\n")
    for _, r in df.iterrows():
        print(f"{r['sweep']}: K={r['K']}, p={r['p_value']}, "
              f"n_test_branches={r['n_test_branches']}, "
              f"expected_subs_on_test_set={r['expected_subs_on_test_set']:.2f}" if r['expected_subs_on_test_set'] is not None
              else f"{r['sweep']}: {r['status']}")
        if r["flag"]:
            print(f"    FLAGS: {r['flag']}")


if __name__ == "__main__":
    main()