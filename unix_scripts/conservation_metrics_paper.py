#!/usr/bin/env python3
"""
conservation_metrics.py
=======================
Extract per-interval evolutionary conservation metrics from BigWig tracks for
CRE foreground and matched background BED files.

Purpose
-------
Following covariate-matched background generation (cre_matcher.py), this script
quantifies evolutionary conservation for each interval in both the CRE foreground
and the matched background set. The resulting per-interval metrics can then be
compared statistically to test whether CREs are more conserved than expected given
their genomic covariates.

Conservation tracks
------------------------------
phastCons
    Per-base probability of evolutionary conservation under a phylogenetic HMM.
    Values are in [0, 1]. Summarised as:
        mean_phastcons_{label}              : mean score across all bases in the interval
        frac_phastcons_{label}_ge_{thr}     : fraction of bases >= threshold (default 0.5)

phyloP
    Per-base log-likelihood ratio score under a neutral substitution model.
    Positive values indicate conservation; negative values indicate acceleration.
    Values span roughly [-20, +30] depending on the clade. Summarised as:
        mean_phyloP_{label}                 : mean score across the interval
        p95_phyloP_{label}                  : 95th-percentile score (captures peaks)
        frac_phyloP_{label}_ge_{thr}        : fraction of bases >= threshold (default 2.0)

Outputs (written to --outdir)
-----------------------------
CREs_metrics.tsv        : Per-CRE conservation metrics (one row per CRE)
BG_metrics.tsv          : Per-background-peak conservation metrics
metrics_meta.json       : Run parameters (tracks, thresholds, thread count)
[optional log file]     : Timestamped progress log if --log-file is specified

Usage example
-------------
python conservation_metrics.py \\
    --cre-bed  results/matched_background/CREs.bed \\
    --bg-bed   results/matched_background/ATAC_matched_1to5.bed \\
    --phastcons  hg38.phastCons100way.bw|verts100 \\
    --phastcons  hg38.phastCons30way.bw|mammals30 \\
    --phylop     hg38.phyloP100way.bw|verts100 \\
    --extra-bw   ENCODE_H3K27ac.bw|H3K27ac \\
    --phastcons-threshold 0.5 \\
    --phylop-threshold    2.0 \\
    --threads 16 \\
    --outdir  results/conservation \\
    --log-file results/conservation/run.log
"""

import argparse
import os
import math
import json
import sys
import traceback
import time
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np


# =============================================================================
# Logging
# =============================================================================

LOG_FH  = None   # Optional file handle for --log-file
VERBOSE = True   # Print progress messages to stdout
QUIET   = False  # Suppress all stdout output (overrides VERBOSE)


def log(msg):
    """
    Write a timestamped log message to stdout and optionally to a log file.

    Parameters
    ----------
    msg : str
        Message text. A timestamp is prepended automatically.

    Notes
    -----
    Output is suppressed when QUIET=True (e.g. when running as a batch job
    whose stdout is discarded). The log file (if open) always receives the
    message regardless of QUIET.
    """
    ts   = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    if not QUIET:
        print(line, flush=True)
    if LOG_FH is not None:
        LOG_FH.write(line + "\n")
        LOG_FH.flush()


# =============================================================================
# I/O helpers
# =============================================================================

def parse_labeled_paths(items):
    """
    Parse a list of "path|label" strings into (path, label) tuples.

    Parameters
    ----------
    items : list of str or None
        Each element is either "path|label" (label taken from after the pipe)
        or just "path" (label defaults to the file basename without extension).

    Returns
    -------
    list of (str, str)
        Pairs of (file_path, track_label). Track labels are used as column
        name suffixes in the output TSVs (e.g., mean_phastcons_verts100).

    Example
    -------
    parse_labeled_paths(["hg38.phastCons100way.bw|verts100"])
    → [("hg38.phastCons100way.bw", "verts100")]
    """
    pairs = []
    for it in items or []:
        if "|" in it:
            path, name = it.split("|", 1)
        else:
            path, name = it, os.path.splitext(os.path.basename(it))[0]
        pairs.append((path, name))
    return pairs


def read_bed(path, max_preview=3):
    """
    Parse a BED file into a pandas DataFrame with columns chrom, start, end, id.

    Parameters
    ----------
    path : str
        Path to the BED file. Lines beginning with '#' are treated as comments.
        Accepts 3-column or 4-column BED; if 4 columns are present, column 4
        is used as the interval ID. If only 3 columns, IDs are auto-generated.
    max_preview : int
        Unused; retained for API compatibility.

    Returns
    -------
    pandas.DataFrame
        Columns: chrom (str), start (int), end (int), id (str).

    Raises
    ------
    SystemExit
        If the file does not exist.
    """
    if not os.path.exists(path):
        raise SystemExit(f"BED file not found: {path}")
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] >= 4:
        df = df.iloc[:, :4]
        df.columns = ["chrom", "start", "end", "id"]
    else:
        df = df.iloc[:, :3]
        df.columns = ["chrom", "start", "end"]
        df["id"] = [f"ID_{i+1}" for i in range(len(df))]
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    log(f"Loaded {path}: {len(df)} intervals.")
    return df


# =============================================================================
# Worker process globals
# =============================================================================
# BigWig file handles (pyBigWig objects) are opened once per worker process in
# _init_worker() and reused across all tasks assigned to that worker. Opening a
# BigWig handle is expensive (reads the index into memory); storing handles as
# process-level globals amortises this cost across many intervals.
# These objects are NOT picklable, so they cannot be passed as task arguments.

G_PCS     = None   # list of (pyBigWig handle, label) for phastCons tracks
G_PLP     = None   # list of (pyBigWig handle, label) for phyloP tracks
G_EXS     = None   # list of (pyBigWig handle, label) for extra BigWig tracks
G_thr_pc  = 0.5    # phastCons threshold for fractional coverage metric
G_thr_pl  = 2.0    # phyloP threshold for fractional coverage metric


def _check_bigwig_openable(path):
    """
    Verify that a BigWig file can be opened and its chromosome table read.

    Called during preflight checks before spawning worker processes, so that
    missing or malformed BigWig files are detected early rather than causing
    cryptic worker crashes mid-run.

    Parameters
    ----------
    path : str
        Path to the BigWig file.

    Returns
    -------
    True on success.

    Raises
    ------
    SystemExit
        If pyBigWig is not installed, or if the file cannot be opened.
    """
    try:
        import pyBigWig
    except ImportError:
        raise SystemExit("pyBigWig is required. Install with: pip install pyBigWig")
    try:
        bw = pyBigWig.open(path)
        _  = bw.chroms()  # Force read of chromosome table to catch truncated files
        bw.close()
        return True
    except Exception as e:
        raise SystemExit(f"Failed to open BigWig: {path}\n{e}")


def _init_worker(pcs_specs, plp_specs, ex_specs, thr_pc, thr_pl, _verbose):
    """
    Initialise worker-process globals: open BigWig handles and set thresholds.

    Called once per worker process via Pool(initializer=_init_worker, ...).
    Each worker opens its own set of file handles independently; this avoids
    concurrent read conflicts that can occur when multiple processes share a
    single file handle.

    Parameters
    ----------
    pcs_specs : list of (str, str)
        (path, label) pairs for phastCons BigWig files.
    plp_specs : list of (str, str)
        (path, label) pairs for phyloP BigWig files.
    ex_specs : list of (str, str)
        (path, label) pairs for extra BigWig files.
    thr_pc : float
        phastCons threshold: bases with score >= thr_pc contribute to
        frac_phastcons_{label}_ge_{thr} columns.
    thr_pl : float
        phyloP threshold: bases with score >= thr_pl contribute to
        frac_phyloP_{label}_ge_{thr} columns.
    _verbose : bool
        Controls per-interval progress logging within the worker.
    """
    global G_PCS, G_PLP, G_EXS, G_thr_pc, G_thr_pl, VERBOSE
    VERBOSE  = _verbose
    import pyBigWig
    G_PCS    = [(pyBigWig.open(p), name) for p, name in pcs_specs]
    G_PLP    = [(pyBigWig.open(p), name) for p, name in plp_specs]
    G_EXS    = [(pyBigWig.open(p), name) for p, name in ex_specs]
    G_thr_pc = thr_pc
    G_thr_pl = thr_pl


# =============================================================================
# Per-interval metric computation
# =============================================================================

def _bw_values(bw, chrom, start, end):
    """
    Retrieve per-base BigWig values for a genomic interval as a dense array.

    Parameters
    ----------
    bw : pyBigWig.open object
        Open BigWig file handle.
    chrom : str
        Chromosome name (must match BigWig chromosome table exactly).
    start : int
        0-based start coordinate (BED convention).
    end : int
        0-based end coordinate, exclusive.

    Returns
    -------
    numpy.ndarray
        1-D float array of per-base scores with NaN positions removed.
        Returns an empty array if the chromosome is absent from the BigWig,
        if the interval is out of range, or if any other error occurs.

    Notes
    -----
    pyBigWig returns NaN for positions not covered by the BigWig track.
    Removing NaNs before computing summary statistics avoids biasing means
    downward in sparsely covered regions; it also means that metrics reflect
    only the *covered* fraction of the interval. The fraction of bases covered
    can be inferred from len(vals) / (end - start) if needed.
    """
    try:
        vals = np.array(bw.values(chrom, int(start), int(end)))
        vals = vals[~np.isnan(vals)]
    except Exception:
        vals = np.array([])
    return vals


def _compute_metrics_for_row(row):
    """
    Compute all conservation metrics for a single genomic interval.

    Uses the worker-process BigWig handles and thresholds stored in module
    globals (G_PCS, G_PLP, G_EXS, G_thr_pc, G_thr_pl). Returns NaN for
    metrics where the BigWig has no coverage over the interval.

    Parameters
    ----------
    row : dict
        Must contain keys: chrom, start, end, id.

    Returns
    -------
    dict
        Contains the interval coordinates and ID, plus one or more of the
        following metric columns (depending on which tracks were supplied):

        For each phastCons track (label = track label from --phastcons):
            mean_phastcons_{label}            : mean per-base score [0, 1]
            frac_phastcons_{label}_ge_{thr}   : fraction of bases >= thr_pc

        For each phyloP track (label = track label from --phylop):
            mean_phyloP_{label}               : mean per-base score
            p95_phyloP_{label}                : 95th percentile score
            frac_phyloP_{label}_ge_{thr}      : fraction of bases >= thr_pl

        For each extra BigWig track (label = track label from --extra-bw):
            mean_{label}                      : mean per-base signal

        On exception: an "ERROR" key with the exception message is added
        instead of metric values, allowing the DataFrame to be written and
        inspected rather than crashing the entire run.

    Notes
    -----
    The 95th-percentile phyloP score (p95) is included in addition to the
    mean because CRE conservation often manifests as short constrained
    sub-elements embedded in less-constrained flanks. The p95 captures
    these peaks more sensitively than a mean that is diluted by flanking
    unconstrained bases.
    """
    try:
        out = {
            "chrom": row["chrom"],
            "start": int(row["start"]),
            "end":   int(row["end"]),
            "id":    row["id"],
        }

        # phastCons metrics
        for bw, name in G_PCS or []:
            vals = _bw_values(bw, row["chrom"], row["start"], row["end"])
            out[f"mean_phastcons_{name}"] = float(
                np.nan if len(vals) == 0 else np.mean(vals)
            )
            out[f"frac_phastcons_{name}_ge_{G_thr_pc}"] = float(
                np.nan if len(vals) == 0 else np.mean(vals >= G_thr_pc)
            )

        # phyloP metrics
        for bw, name in G_PLP or []:
            vals = _bw_values(bw, row["chrom"], row["start"], row["end"])
            out[f"mean_phyloP_{name}"] = float(
                np.nan if len(vals) == 0 else np.mean(vals)
            )
            out[f"p95_phyloP_{name}"] = float(
                np.nan if len(vals) == 0 else np.percentile(vals, 95)
            )
            out[f"frac_phyloP_{name}_ge_{G_thr_pl}"] = float(
                np.nan if len(vals) == 0 else np.mean(vals >= G_thr_pl)
            )

        # Extra BigWig tracks
        for bw, name in G_EXS or []:
            vals = _bw_values(bw, row["chrom"], row["start"], row["end"])
            out[f"mean_{name}"] = float(
                np.nan if len(vals) == 0 else np.mean(vals)
            )

        return out

    except Exception as e:
        # Return a partial dict with an ERROR field rather than letting the
        # worker crash. Errors are visible in the output TSV for inspection.
        return {
            "chrom": row.get("chrom", ""),
            "start": int(row.get("start", 0)),
            "end":   int(row.get("end", 0)),
            "id":    row.get("id", ""),
            "ERROR": f"{type(e).__name__}: {e}",
        }


# =============================================================================
# Main entry point
# =============================================================================

def main():
    """
    Parse command-line arguments and orchestrate conservation metric extraction.

    Pipeline steps
    --------------
    1.  Preflight: load BED files, validate BigWig accessibility.
    2.  For each input BED (CREs and background), compute per-interval metrics
        across all supplied BigWig tracks, parallelised across intervals.
    3.  Write per-interval metric TSVs and a JSON metadata file.
    """
    ap = argparse.ArgumentParser(
        description=(
            "Extract per-interval evolutionary conservation metrics from BigWig "
            "tracks for CRE foreground and matched background BED files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required BED inputs
    ap.add_argument("--cre-bed", required=True,
                    help="BED file of foreground CREs.")
    ap.add_argument("--bg-bed",  required=True,
                    help="BED file of matched background peaks (e.g., ATAC_matched_1to5.bed).")

    # Conservation track inputs (repeatable; each as "path|label")
    ap.add_argument("--phastcons", action="append", default=[],
                    help='phastCons BigWig: "path/to/file.bw|label". Repeatable.')
    ap.add_argument("--phylop",    action="append", default=[],
                    help='phyloP BigWig:    "path/to/file.bw|label". Repeatable.')
    ap.add_argument("--extra-bw",  action="append", default=[],
                    help='Any additional BigWig: "path/to/file.bw|label". Repeatable.')

    # Threshold parameters for fractional coverage metrics
    ap.add_argument("--phastcons-threshold", type=float, default=0.5,
                    help="phastCons score threshold for frac_ge metric. Default: 0.5.")
    ap.add_argument("--phylop-threshold",    type=float, default=2.0,
                    help="phyloP score threshold for frac_ge metric. Default: 2.0.")

    # Parallelisation
    ap.add_argument("--threads", type=int, default=max(1, cpu_count() // 2),
                    help="Worker processes. Default: half of available CPUs.")

    # Output
    ap.add_argument("--outdir",    default="results",
                    help="Output directory (created if absent). Default: 'results'.")

    # Logging and diagnostic options
    ap.add_argument("--verbose",    action="store_true",
                    help="Print per-interval progress every 1000/2000 intervals.")
    ap.add_argument("--quiet",      action="store_true",
                    help="Suppress all stdout output.")
    ap.add_argument("--log-file",   default="",
                    help="Optional log file path. Receives all log() output.")
    ap.add_argument("--check-only", action="store_true",
                    help="Run preflight checks only (BigWig openability, BED existence) "
                         "then exit. Useful before submitting long cluster jobs.")

    args = ap.parse_args()

    # ── Logging setup ──────────────────────────────────────────────────────────
    global LOG_FH, VERBOSE, QUIET
    VERBOSE = args.verbose or True  # Default verbose=True for transparency
    QUIET   = args.quiet
    if args.log_file:
        LOG_FH = open(args.log_file, "w")

    os.makedirs(args.outdir, exist_ok=True)

    # ── Step 1: Preflight checks ──────────────────────────────────────────────
    log("=== Preflight checks ===")
    cres = read_bed(args.cre_bed)
    bg   = read_bed(args.bg_bed)

    pcs_specs = parse_labeled_paths(args.phastcons)
    plp_specs = parse_labeled_paths(args.phylop)
    ex_specs  = parse_labeled_paths(args.extra_bw)

    if not (pcs_specs or plp_specs or ex_specs):
        raise SystemExit(
            "No BigWig tracks specified. Use --phastcons, --phylop, or --extra-bw."
        )
    log(f"Tracks: phastCons={len(pcs_specs)}, phyloP={len(plp_specs)}, "
        f"extra={len(ex_specs)}")

    # Validate all BigWig files before spawning workers
    for p, _ in pcs_specs + plp_specs + ex_specs:
        log(f"Checking BigWig: {p}")
        _check_bigwig_openable(p)
    log("All BigWig checks passed.")

    log(f"CREs={len(cres)}  BG={len(bg)}  threads={args.threads}")

    if args.check_only:
        log("--check-only mode: exiting after preflight.")
        return

    # ── Step 2: Extract metrics ───────────────────────────────────────────────

    def run_set(name, df):
        """
        Compute conservation metrics for all intervals in a BED DataFrame.

        Parameters
        ----------
        name : str
            Label for logging (e.g. "CREs" or "BG").
        df : pandas.DataFrame
            BED DataFrame with columns: chrom, start, end, id.

        Returns
        -------
        pandas.DataFrame
            One row per interval with all metric columns.

        Notes
        -----
        In parallel mode, P.imap_unordered() is used rather than P.map() to
        return results as they complete, enabling progress logging at
        fine-grained intervals without waiting for the full chunk. chunksize=200
        balances task-dispatch overhead against load-balancing granularity.
        """
        log(f"=== {name} === starting ({len(df)} intervals)")
        rows = df.to_dict(orient="records")

        if args.threads <= 1:
            # Sequential mode: initialise globals in the main process
            _init_worker(
                pcs_specs, plp_specs, ex_specs,
                args.phastcons_threshold, args.phylop_threshold, VERBOSE
            )
            out = []
            for i, r in enumerate(rows, 1):
                if VERBOSE and i % 1000 == 0:
                    log(f"{name}: processed {i}/{len(rows)}")
                out.append(_compute_metrics_for_row(r))
        else:
            with Pool(
                processes=args.threads,
                initializer=_init_worker,
                initargs=(
                    pcs_specs, plp_specs, ex_specs,
                    args.phastcons_threshold, args.phylop_threshold, VERBOSE
                ),
            ) as P:
                out = []
                for i, res in enumerate(
                    P.imap_unordered(_compute_metrics_for_row, rows, chunksize=200), 1
                ):
                    if VERBOSE and i % 2000 == 0:
                        log(f"{name}: processed {i}/{len(rows)}")
                    out.append(res)

        log(f"=== {name} === done")
        return pd.DataFrame(out)

    df_cre = run_set("CREs", cres)
    df_bg  = run_set("BG",   bg)

    # ── Step 3: Write outputs ─────────────────────────────────────────────────
    cre_out = os.path.join(args.outdir, "CREs_metrics.tsv")
    bg_out  = os.path.join(args.outdir, "BG_metrics.tsv")
    df_cre.to_csv(cre_out, sep="\t", index=False)
    df_bg.to_csv(bg_out,  sep="\t", index=False)

    # Save run parameters for reproducibility and audit trail
    meta = {
        "phastcons":            pcs_specs,
        "phylop":               plp_specs,
        "extra_bw":             ex_specs,
        "phastcons_threshold":  args.phastcons_threshold,
        "phylop_threshold":     args.phylop_threshold,
        "threads":              args.threads,
    }
    with open(os.path.join(args.outdir, "metrics_meta.json"), "w") as f:
        f.write(json.dumps(meta, indent=2))

    log(f"Wrote: {cre_out}")
    log(f"Wrote: {bg_out}")
    log(f"Wrote: {os.path.join(args.outdir, 'metrics_meta.json')}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        if str(e):
            print(f"[exit] {e}", file=sys.stderr)
        raise
    except Exception as e:
        print("[fatal]", repr(e), file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
