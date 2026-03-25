#!/usr/bin/env python3
"""
cre_matcher.py
==============
Covariate-matched background sampler for cis-regulatory element (CRE) conservation analyses.

Purpose
-------
When testing whether a set of CREs is enriched for evolutionary conservation (or any
other genomic signal), a naïve comparison against all open-chromatin peaks is confounded
by systematic differences in interval length, GC content, distance to the nearest TSS,
and chromatin accessibility (ATAC signal). This script constructs a background set of
ATAC-seq peaks that is matched to the CRE foreground on all four covariates simultaneously,
ensuring that downstream enrichment signals reflect genuine conservation rather than
covariate bias.

Matching strategy
-----------------
1.  Non-overlap filtering — ATAC peaks that overlap any CRE are removed from the
    background pool, preventing the same locus appearing in both foreground and background.

2.  Covariate binning — Both CREs and remaining ATAC peaks are jointly assigned to a
    quantile-bin grid defined by:
        • interval length           (default 10 quantile bins)
        • GC content                (default 10 quantile bins)
        • log10(TSS distance + 1)   (default 10 quantile bins)
        • mean ATAC signal          (default  4 quantile bins; optional)
    Bin edges are derived from the *pooled* CRE + ATAC distribution so that both sets
    share the same reference frame.

3.  Matched sampling — For each CRE, a random background peak is drawn from the bin
    cell containing that CRE (1-to-1 match) or up to N peaks are drawn (1-to-N, default
    N=5). If a CRE falls in an empty bin cell, a kNN fallback in z-scored covariate space
    selects the nearest available peak(s).

All steps are parallelised across chromosomes (non-overlap) or across CREs (matching)
using Python's multiprocessing.Pool.

Outputs (written to --outdir)
-----------------------------
CREs_covariates.tsv      Per-CRE covariate values and bin assignments
ATAC_covariates.tsv      Per-ATAC-peak covariate values and bin assignments (post-filter)
bins_summary.tsv         Number of CREs and background peaks per bin cell
matches_1to1.tsv         1-to-1 CRE → ATAC matched pairs
matches_1to5.tsv         1-to-N CRE → ATAC matched pairs (N = --n-per-cre)
ATAC_matched_1to1.bed    BED file of background peaks selected in 1-to-1 matching
ATAC_matched_1to5.bed    BED file of background peaks selected in 1-to-N matching
bin_edges.json           Quantile bin edges for each covariate (for reproducibility)

Usage example
-------------
python cre_matcher.py \\
    --cre-bed   heart_cres.bed \\
    --atac-bed  all_atac_peaks.bed \\
    --cre-gc    heart_cres_gc.bed \\
    --atac-gc   all_atac_gc.bed \\
    --cre-tssd  heart_cres_tssd.bed \\
    --atac-tssd all_atac_tssd.bed \\
    --cre-atac  heart_cres_atac.bed \\
    --atac-atac all_atac_signal.bed \\
    --len-bins 10 --gc-bins 10 --tssd-bins 10 --atac-bins 4 \\
    --n-per-cre 5 --threads 8 --seed 42 \\
    --outdir results/matched_background

Input format notes
------------------
--cre-bed / --atac-bed
    Standard 4-column BED (chrom start end id). A unique ID in column 4 is strongly
    recommended; if absent, IDs are auto-generated as CRE_1, CRE_2, ... etc.

--cre-gc / --atac-gc
    BED-like file with GC fraction in column 4.
    Can be generated with: bedtools nuc -fi genome.fa -bed input.bed

--cre-tssd / --atac-tssd
    BED-like file with distance to the nearest annotated TSS in column 4 (bp).
    Can be generated with:
        bedtools closest -a input.bed -b tss.bed -d | awk '{OFS="\\t"; print $1,$2,$3,$NF}'

--cre-atac / --atac-atac (optional)
    BED-like file with mean ATAC-seq signal per interval in column 4.
    Can be generated with: bigWigAverageOverBed atac.bw input.bed /dev/stdout

"""

import argparse
import os
import math
import random
import json
from collections import defaultdict, Counter
import csv
from multiprocessing import Pool, cpu_count


# =============================================================================
# I/O helpers
# =============================================================================

def read_bed3(path, with_id=True, prefix="ID_"):
    """
    Parse a BED file and return a list of (chrom, start, end, id) tuples.

    Parameters
    ----------
    path : str
        Path to the BED file. Lines beginning with '#' and blank lines are skipped.
        At least 3 columns (chrom, start, end) are required.
    with_id : bool
        If True and a 4th column is present, use it as the region identifier.
        If False or the 4th column is absent, IDs are auto-generated as
        f"{prefix}{i+1}" (1-based).
    prefix : str
        Prefix for auto-generated IDs.

    Returns
    -------
    list of (str, int, int, str)
        Each tuple is (chrom, start, end, region_id).

    Notes
    -----
    Using unique, stable IDs in column 4 is strongly recommended so that the output
    mapping tables (matches_1to1.tsv, matches_1to5.tsv) are unambiguously
    interpretable across runs and downstream tools.
    """
    rows = []
    with open(path) as f:
        i = 0
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[:3]
            start = int(start)
            end   = int(end)
            rid   = parts[3] if (with_id and len(parts) >= 4) else f"{prefix}{i+1}"
            rows.append((chrom, start, end, rid))
            i += 1
    return rows


def read_bed_value(path):
    """
    Parse a BED-like file where column 4 holds a numeric value per interval.

    Parameters
    ----------
    path : str
        Path to a tab-separated file: chrom, start, end, value [, extra_cols...].
        If column 4 is non-numeric (e.g. a header token), the *last* column is
        tried — this makes the function robust to extra annotation columns.

    Returns
    -------
    dict
        Keys are (chrom, start, end) tuples; values are floats.
        Rows that cannot be parsed to a float in any tried column are silently
        skipped and will appear as NaN when joined via add_covariates().

    Typical use
    -----------
    Covariates such as GC content, TSS distance, or mean ATAC signal are stored
    as one-line-per-interval files and loaded with this function. They are then
    joined to the main interval table via (chrom, start, end) keys in
    add_covariates().
    """
    d = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 4:
                continue
            chrom, start, end = parts[:3]
            key    = (chrom, int(start), int(end))
            valtok = parts[3]
            try:
                v = float(valtok)
            except ValueError:
                try:
                    v = float(parts[-1])
                except ValueError:
                    continue
            d[key] = v
    return d


def add_covariates(bed_rows, gc_map, tssd_map, atac_map=None):
    """
    Attach covariate values to each interval by joining on (chrom, start, end) keys.

    Parameters
    ----------
    bed_rows : list of (chrom, start, end, id)
        Parsed BED rows from read_bed3().
    gc_map : dict
        (chrom, start, end) → GC fraction, from read_bed_value().
    tssd_map : dict
        (chrom, start, end) → TSS distance in bp, from read_bed_value().
    atac_map : dict or None
        Optional (chrom, start, end) → mean ATAC signal.
        If None, ATAC is stored as float('nan') for all intervals.

    Returns
    -------
    list of dict
        Each dict has keys: chrom, start, end, id, length, gc, tssd, atac.
        Unmatched covariate values are stored as float('nan') so they are
        handled gracefully during binning (assign_bin returns None for NaN,
        triggering the kNN fallback at matching time).

    Notes
    -----
    Interval length is computed directly from (end - start) and is always
    available regardless of whether a length covariate file is provided.
    """
    out = []
    for chrom, start, end, rid in bed_rows:
        key    = (chrom, start, end)
        length = end - start
        gc     = gc_map.get(key, float('nan'))
        tssd   = tssd_map.get(key, float('nan'))
        atac   = atac_map.get(key, float('nan')) if atac_map else float('nan')
        out.append({
            "chrom":  chrom,
            "start":  start,
            "end":    end,
            "id":     rid,
            "length": float(length),
            "gc":     gc,
            "tssd":   float(tssd),
            "atac":   atac,
        })
    return out


def write_tsv(path, rows, header):
    """
    Write a list of dicts to a tab-separated file.

    Parameters
    ----------
    path : str
        Output file path.
    rows : list of dict
        Each dict should contain keys corresponding to `header`. Missing keys
        are written as empty strings.
    header : list of str
        Column names — defines both the header row and the column order.
    """
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(h, "") for h in header])


# =============================================================================
# Covariate binning
# =============================================================================

def quantile_bins(values, nbins):
    """
    Compute quantile-based bin edges for a 1-D numeric array.

    Quantile binning is preferred over equal-width binning because covariate
    distributions are often highly skewed (e.g., TSS distance spans several
    orders of magnitude). Quantile bins ensure approximately equal numbers of
    intervals per bin, maximising the density of the background pool in each
    cell and reducing the fraction of CREs that require the kNN fallback.

    Parameters
    ----------
    values : list of float
        Input values; NaN values are excluded before computing quantiles.
    nbins : int
        Number of bins requested. If <= 1 or input is empty, returns a single
        catch-all bin (-inf, +inf) so the covariate is effectively ignored.

    Returns
    -------
    list of float
        Sorted list of (nbins + 1) edge values, always spanning (-inf, +inf).
        Duplicate edges are collapsed, so the actual number of usable bins may
        be less than nbins for low-cardinality data.

    Notes
    -----
    Uses linear interpolation for non-integer quantile positions, consistent
    with numpy.percentile and pandas.quantile default behaviour.
    """
    xs = sorted([v for v in values if not math.isnan(v)])
    if nbins <= 1 or len(xs) == 0:
        return [-float('inf'), float('inf')]
    edges = [-float('inf')]
    for i in range(1, nbins):
        qpos  = (len(xs) - 1) * i / nbins
        lower = int(math.floor(qpos))
        upper = int(math.ceil(qpos))
        if lower == upper:
            qv = xs[lower]
        else:
            frac = qpos - lower
            qv   = xs[lower] * (1 - frac) + xs[upper] * frac
        # Only append if strictly greater than the previous edge to prevent
        # zero-width bins in heavily tied distributions.
        if qv > edges[-1]:
            edges.append(qv)
    edges.append(float('inf'))
    return edges


def assign_bin(value, edges):
    """
    Return the 0-based bin index for a single value given a set of bin edges.

    Uses binary search — O(log k) for k bins — to locate the correct bin.

    Parameters
    ----------
    value : float
        Value to bin. If NaN, returns None; intervals with None bin assignments
        are excluded from bin-cell matching and fall back to kNN.
    edges : list of float
        Sorted bin edges as returned by quantile_bins().

    Returns
    -------
    int or None
        0-based bin index, or None for NaN input.
    """
    if math.isnan(value):
        return None
    lo = 0
    hi = len(edges) - 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if value < edges[mid]:
            hi = mid
        else:
            lo = mid
    return lo


def build_bins_table(rows, len_bins, gc_bins, tssd_bins, atac_bins):
    """
    Assign all intervals to their respective multi-dimensional bin cells and
    return the bin-edge metadata.

    TSS distance is log10-transformed before binning because the distribution
    spans several orders of magnitude (tens of bp near promoters to >1 Mb for
    distal enhancers). The transform compresses the long tail so quantile edges
    are placed more evenly across the biologically relevant range.

    Parameters
    ----------
    rows : list of dict
        Interval dicts as returned by add_covariates(). Each dict is modified
        *in place* to add keys: len_bin, gc_bin, tssd_bin, atac_bin.
    len_bins, gc_bins, tssd_bins, atac_bins : int
        Number of quantile bins for each covariate.
        Set atac_bins <= 1 to disable ATAC binning (all intervals map to one bin,
        effectively ignoring the ATAC dimension).

    Returns
    -------
    rows : list of dict
        Input rows with bin assignments added.
    meta : dict
        Bin edges under keys: len_edges, gc_edges, tssd_log10_edges, atac_edges.
        Written to bin_edges.json for downstream reproducibility and to allow
        re-binning of new datasets against the same reference frame.

    Notes
    -----
    In main(), this function is called on the *pooled* CRE + ATAC rows so that
    both sets share identical bin boundaries, which is required for the bin-cell
    matching strategy to be valid.
    """
    len_edges  = quantile_bins([r["length"] for r in rows], len_bins)
    gc_edges   = quantile_bins([r["gc"] for r in rows], gc_bins)
    tssd_edges = quantile_bins([math.log10(r["tssd"] + 1.0) for r in rows], tssd_bins)

    # ATAC is optional; disable if no valid values or atac_bins <= 1.
    atac_vals  = [r["atac"] for r in rows if not math.isnan(r["atac"])]
    atac_edges = (
        quantile_bins(atac_vals, atac_bins)
        if (atac_bins and atac_bins > 1 and atac_vals)
        else [-float('inf'), float('inf')]
    )

    for r in rows:
        r["len_bin"]  = assign_bin(r["length"],                   len_edges)
        r["gc_bin"]   = assign_bin(r["gc"],                        gc_edges)
        r["tssd_bin"] = assign_bin(math.log10(r["tssd"] + 1.0),   tssd_edges)
        r["atac_bin"] = assign_bin(r["atac"],                      atac_edges)

    return rows, {
        "len_edges":        len_edges,
        "gc_edges":         gc_edges,
        "tssd_log10_edges": tssd_edges,
        "atac_edges":       atac_edges,
    }


# =============================================================================
# Non-overlap filtering (parallelised per chromosome)
# =============================================================================

def _nonoverlap_one_chrom(args):
    """
    Worker: remove ATAC peaks that overlap any CRE on a single chromosome.

    Overlap is defined as any base-pair overlap:
        atac_start < cre_end  AND  cre_start < atac_end

    Parameters
    ----------
    args : tuple
        (chrom, atac_list, cre_intervals)
        chrom         — chromosome name (str; used for return identification only)
        atac_list     — list of ATAC peak dicts for this chromosome
        cre_intervals — sorted list of (start, end) tuples for CREs on this chrom

    Returns
    -------
    (chrom, kept)
        kept is the subset of atac_list with no CRE overlap.

    Algorithm
    ---------
    For each ATAC peak, binary search locates the first CRE whose end > atac_start.
    All CREs before this index end before the ATAC peak starts and cannot overlap.
    A short linear scan then checks subsequent CREs until their start >= atac_end.
    Complexity: O(N log M + total_overlap_count) per chromosome vs. O(N*M) naïve.
    """
    chrom, atac_list, cre_intervals = args
    kept = []
    arr  = cre_intervals  # sorted list of (start, end)

    for a in atac_list:
        s, e = a["start"], a["end"]
        overlaps = False

        # Binary search: find first CRE index whose end > atac_start.
        lo, hi = 0, len(arr)
        while lo < hi:
            mid = (lo + hi) // 2
            if arr[mid][1] <= s:
                lo = mid + 1
            else:
                hi = mid
        idx = lo

        # Linear scan from first candidate.
        while idx < len(arr) and arr[idx][0] < e:
            cs, ce = arr[idx]
            if cs < e and ce > s:
                overlaps = True
                break
            idx += 1

        if not overlaps:
            kept.append(a)

    return chrom, kept


def enforce_nonoverlap_parallel(atac_rows, cre_rows, threads):
    """
    Remove ATAC peaks overlapping any CRE, parallelised across chromosomes.

    Parallelisation is by chromosome because intervals on different chromosomes
    are fully independent, enabling embarrassingly parallel processing. For a
    typical mammalian genome (25–30 chromosomes), this provides near-linear
    speedup up to ~25 threads.

    Parameters
    ----------
    atac_rows : list of dict
        All ATAC peak dicts from add_covariates().
    cre_rows  : list of dict
        All CRE dicts. Only (chrom, start, end) are used.
    threads : int
        Number of worker processes. If <= 1, runs sequentially (useful for
        debugging or memory-constrained environments).

    Returns
    -------
    list of dict
        ATAC peaks with no CRE overlap, across all chromosomes.
    """
    bychr_cre  = defaultdict(list)
    bychr_atac = defaultdict(list)

    for r in cre_rows:
        bychr_cre[r["chrom"]].append((r["start"], r["end"]))
    for a in atac_rows:
        bychr_atac[a["chrom"]].append(a)

    # Sort CRE intervals per chromosome to enable the binary search in worker.
    for c in bychr_cre:
        bychr_cre[c].sort()

    tasks = [
        (c, bychr_atac.get(c, []), bychr_cre.get(c, []))
        for c in bychr_atac.keys()
    ]

    if threads <= 1:
        outs = [_nonoverlap_one_chrom(t) for t in tasks]
    else:
        with Pool(processes=threads) as P:
            outs = P.map(_nonoverlap_one_chrom, tasks)

    kept = []
    for chrom, lst in outs:
        kept.extend(lst)
    return kept


# =============================================================================
# Covariate matching (parallelised per CRE)
# =============================================================================

def scale(vals):
    """
    Compute mean and standard deviation of a list of values, ignoring NaNs.

    Used to z-score covariates for the kNN fallback distance calculation,
    placing all four dimensions on a common scale before computing Euclidean
    distances.

    Parameters
    ----------
    vals : list of float
        Input values; NaN entries are excluded.

    Returns
    -------
    (mean, sd) : (float, float)
        Returns (0.0, 1.0) if no valid values are present (safe no-op z-score).
        Returns (mean, 1.0) if variance is zero (constant-valued covariate).
    """
    xs = [v for v in vals if not math.isnan(v)]
    if not xs:
        return (0.0, 1.0)
    m   = sum(xs) / len(xs)
    var = sum((v - m) * (v - m) for v in xs) / max(1, len(xs) - 1)
    sd  = math.sqrt(var) if var > 0 else 1.0
    return (m, sd)


def prepare_knn_scaling(atac_rows):
    """
    Pre-compute z-score parameters (mean, sd) for each covariate from the ATAC pool.

    These parameters are used in dist_point() to place all four covariates on a
    common scale before computing Euclidean distances in the kNN fallback.

    Parameters
    ----------
    atac_rows : list of dict
        The non-overlapping ATAC pool after enforce_nonoverlap_parallel().

    Returns
    -------
    dict with keys mL, sL, mG, sG, mT, sT, mA, sA
        Mean (m) and standard deviation (s) for Length (L), GC (G),
        log10-TSS-distance (T), and ATAC signal (A).

    Notes
    -----
    Scaling is derived from the ATAC pool rather than the CREs to avoid
    contaminating the background reference frame with foreground statistics.
    """
    L = [p["length"] for p in atac_rows]
    G = [p["gc"]     for p in atac_rows]
    T = [math.log10(p["tssd"] + 1.0) for p in atac_rows]
    A = [p["atac"]   for p in atac_rows]
    return {
        "mL": scale(L)[0], "sL": scale(L)[1],
        "mG": scale(G)[0], "sG": scale(G)[1],
        "mT": scale(T)[0], "sT": scale(T)[1],
        "mA": scale(A)[0], "sA": scale(A)[1],
    }


def z(v, m, s):
    """
    Z-score a single value. Returns 0.0 for NaN inputs (safe default that
    contributes zero to the Euclidean distance for missing covariates).

    Parameters
    ----------
    v : float  — raw value
    m : float  — population mean
    s : float  — population standard deviation (returns 0.0 if s == 0)
    """
    if math.isnan(v):
        return 0.0
    return (v - m) / s if s > 0 else 0.0


def dist_point(cre, p, S):
    """
    Squared Euclidean distance between a CRE and a candidate ATAC peak in
    z-scored 4-D covariate space (length, GC, log10-TSS-distance, ATAC signal).

    Used only in the kNN fallback when a CRE's bin cell contains no ATAC peaks.
    Squared distance is sufficient for ranking (avoids sqrt computation).

    Parameters
    ----------
    cre : dict  — CRE covariate dict (keys: length, gc, tssd, atac)
    p   : dict  — ATAC peak covariate dict
    S   : dict  — scaling parameters from prepare_knn_scaling()

    Returns
    -------
    float : squared Euclidean distance in z-scored 4-D space.
    """
    c_vec = (
        z(cre["length"],                   S["mL"], S["sL"]),
        z(cre["gc"],                        S["mG"], S["sG"]),
        z(math.log10(cre["tssd"] + 1.0),   S["mT"], S["sT"]),
        z(cre["atac"],                      S["mA"], S["sA"]),
    )
    p_vec = (
        z(p["length"],                     S["mL"], S["sL"]),
        z(p["gc"],                          S["mG"], S["sG"]),
        z(math.log10(p["tssd"] + 1.0),     S["mT"], S["sT"]),
        z(p["atac"],                        S["mA"], S["sA"]),
    )
    return sum((a - b) * (a - b) for a, b in zip(c_vec, p_vec))


# ---------------------------------------------------------------------------
# Worker process globals
# ---------------------------------------------------------------------------
# These module-level variables are populated in each worker process by
# _init_worker() via Pool(initializer=...). They cannot be passed as regular
# task arguments to _match_one() because multiprocessing serialises arguments
# via pickle for every task, which would be prohibitively slow for large
# pool_by_bin and atac_rows structures. Storing them as process-local globals
# avoids repeated pickling overhead — each worker receives a single copy via
# the initializer and retains it for the lifetime of the process.

G_pool_by_bin = None   # dict: (len_bin, gc_bin, tssd_bin, atac_bin) → list of ATAC dicts
G_atac_rows   = None   # flat list of all ATAC rows (used for kNN fallback)
G_scaling     = None   # z-score parameters from prepare_knn_scaling()
G_seed        = 1      # base random seed; per-CRE seed = G_seed + idx * 9973
G_n_per       = 5      # number of background peaks to sample per CRE (1-to-N)


def _init_worker(pool_by_bin, atac_rows, scaling, seed, n_per):
    """
    Initialise worker-process globals before processing CRE matching tasks.

    Called once per worker process by Pool(initializer=_init_worker, ...).
    Populates the module-level globals consumed by _match_one().

    Parameters
    ----------
    pool_by_bin : dict
        Pre-built lookup: bin-cell 4-tuple → list of ATAC peak dicts.
    atac_rows : list of dict
        Full ATAC pool (for kNN fallback when a bin cell is empty).
    scaling : dict
        Z-score parameters from prepare_knn_scaling().
    seed : int
        Base random seed for deterministic per-CRE sampling.
    n_per : int
        Number of background peaks to draw per CRE in 1-to-N matching.
    """
    global G_pool_by_bin, G_atac_rows, G_scaling, G_seed, G_n_per
    G_pool_by_bin = pool_by_bin
    G_atac_rows   = atac_rows
    G_scaling     = scaling
    G_seed        = seed
    G_n_per       = n_per


def _match_one(args):
    """
    Match a single CRE to background ATAC peaks (both 1:1 and 1:N).

    Called in parallel by Pool.map(), one invocation per CRE.

    Matching proceeds in two stages:
    1.  Bin-cell match (primary): Draw at random from ATAC peaks in the same
        4-D bin cell as the CRE. This ensures covariate-matched backgrounds
        in O(1) lookup time.
    2.  kNN fallback: If the CRE's bin cell is empty (unusual covariate
        combination), the nearest neighbours in z-scored 4-D covariate space
        are selected instead. This prevents unmatched CREs at the cost of
        slightly less precise covariate matching for edge cases.

    Reproducibility
    ---------------
    Each CRE uses an independent deterministic seed derived from the base seed
    and the CRE's index: seed = G_seed + idx * 9973. The large prime multiplier
    spreads seeds widely to avoid correlated random sequences across adjacent
    CREs while keeping the scheme reproducible given the same --seed argument.

    Parameters
    ----------
    args : (int, dict)
        (cre_index, cre_dict), where cre_dict contains covariate values and
        bin assignments as produced by assign_all() in main().

    Returns
    -------
    (one, many) : tuple
        one  — (cre_id, atac_id) for the 1-to-1 match, or (cre_id, None) if
               no ATAC peaks exist anywhere (degenerate edge case).
        many — list of (cre_id, atac_id) for the 1-to-N matches.
    """
    idx, cre = args
    rnd  = random.Random(G_seed + idx * 9973)  # Deterministic per-CRE RNG

    key  = (cre["len_bin"], cre["gc_bin"], cre["tssd_bin"], cre["atac_bin"])
    pool = G_pool_by_bin.get(key, [])

    # ── 1-to-1 match ──────────────────────────────────────────────────────────
    if pool:
        pick1 = rnd.choice(pool)
    else:
        # kNN fallback: sort entire ATAC pool by covariate distance, take closest
        nn    = sorted(G_atac_rows, key=lambda p: dist_point(cre, p, G_scaling))
        pick1 = nn[0] if nn else None
    one = (cre["id"], pick1["id"] if pick1 else None)

    # ── 1-to-N match ──────────────────────────────────────────────────────────
    if pool:
        k         = min(G_n_per, len(pool))
        picks     = rnd.sample(pool, k)      # sample without replacement within bin cell
        picks_ids = [p["id"] for p in picks]
    else:
        # kNN fallback: take the G_n_per nearest neighbours
        nn        = sorted(G_atac_rows, key=lambda p: dist_point(cre, p, G_scaling))[:G_n_per]
        picks_ids = [p["id"] for p in nn]

    many = [(cre["id"], aid) for aid in picks_ids]
    return one, many


# =============================================================================
# Main entry point
# =============================================================================

def main():
    """
    Parse command-line arguments and orchestrate the full matching pipeline.

    Pipeline steps
    --------------
    1.  Read CRE and ATAC BED files and their covariate annotation files.
    2.  Remove ATAC peaks overlapping CREs (parallelised per chromosome).
    3.  Compute quantile bin edges from the pooled CRE + ATAC distribution.
    4.  Assign all intervals to bin cells.
    5.  Match each CRE to 1 (1:1) and N (1:N) background ATAC peaks.
    6.  Write all output tables, BED files, and bin-edge JSON.
    """
    ap = argparse.ArgumentParser(
        description=(
            "Parallel CRE–ATAC covariate binning and matched background sampler. "
            "Produces covariate-matched background ATAC peaks for conservation "
            "enrichment analyses, controlling for length, GC content, TSS distance, "
            "and (optionally) ATAC signal."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required inputs
    ap.add_argument("--cre-bed",   required=True,
                    help="BED file of foreground CREs (>=3 cols; col 4 = unique ID recommended).")
    ap.add_argument("--atac-bed",  required=True,
                    help="BED file of background ATAC peaks (>=3 cols; col 4 = unique ID recommended).")
    ap.add_argument("--cre-gc",    required=True,
                    help="BED-like file; GC fraction in col 4 for each CRE.")
    ap.add_argument("--atac-gc",   required=True,
                    help="BED-like file; GC fraction in col 4 for each ATAC peak.")
    ap.add_argument("--cre-tssd",  required=True,
                    help="BED-like file; TSS distance (bp) in col 4 for each CRE.")
    ap.add_argument("--atac-tssd", required=True,
                    help="BED-like file; TSS distance (bp) in col 4 for each ATAC peak.")

    # Optional ATAC signal covariates
    ap.add_argument("--cre-atac",  default=None,
                    help="Optional BED-like; mean ATAC signal in col 4 for each CRE.")
    ap.add_argument("--atac-atac", default=None,
                    help="Optional BED-like; mean ATAC signal in col 4 for each ATAC peak.")

    # Binning parameters
    ap.add_argument("--len-bins",  type=int, default=10,
                    help="Quantile bins for interval length. Default: 10.")
    ap.add_argument("--gc-bins",   type=int, default=10,
                    help="Quantile bins for GC content. Default: 10.")
    ap.add_argument("--tssd-bins", type=int, default=10,
                    help="Quantile bins for log10(TSS distance + 1). Default: 10.")
    ap.add_argument("--atac-bins", type=int, default=4,
                    help="Quantile bins for ATAC signal. Set to 1 to disable. Default: 4.")

    # Matching parameters
    ap.add_argument("--n-per-cre", type=int, default=5,
                    help="Background peaks per CRE in 1:N matching. Default: 5.")
    ap.add_argument("--threads",   type=int, default=max(1, cpu_count() // 2),
                    help="Parallel worker processes. Default: half of available CPUs.")
    ap.add_argument("--seed",      type=int, default=1,
                    help="Base random seed for reproducible matching. Default: 1.")

    # Output
    ap.add_argument("--outdir",    default="results",
                    help="Output directory (created if absent). Default: 'results'.")

    args = ap.parse_args()

    random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # ── Step 1: Read all inputs ───────────────────────────────────────────────
    print("Reading input files...")
    cre_bed  = read_bed3(args.cre_bed,  with_id=True, prefix="CRE_")
    atac_bed = read_bed3(args.atac_bed, with_id=True, prefix="ATAC_")

    cre_gc    = read_bed_value(args.cre_gc)
    atac_gc   = read_bed_value(args.atac_gc)
    cre_tssd  = read_bed_value(args.cre_tssd)
    atac_tssd = read_bed_value(args.atac_tssd)
    cre_atac  = read_bed_value(args.cre_atac)  if args.cre_atac  else None
    atac_atac = read_bed_value(args.atac_atac) if args.atac_atac else None

    cre_rows  = add_covariates(cre_bed,  cre_gc,  cre_tssd,  cre_atac)
    atac_rows = add_covariates(atac_bed, atac_gc, atac_tssd, atac_atac)
    print(f"  CREs: {len(cre_rows)}  |  ATAC peaks (before filter): {len(atac_rows)}")

    # ── Step 2: Non-overlap filtering ─────────────────────────────────────────
    print("Filtering ATAC peaks overlapping CREs (parallel, per chromosome)...")
    atac_rows = enforce_nonoverlap_parallel(atac_rows, cre_rows, threads=args.threads)
    print(f"  ATAC peaks after non-overlap filter: {len(atac_rows)}")

    # ── Steps 3 & 4: Compute bin edges and assign intervals to bins ───────────
    # Edges are derived from the *pooled* distribution so that CREs and ATAC
    # peaks share the same bin boundaries — a prerequisite for valid bin-cell
    # matching.
    print("Computing quantile bin edges (pooled CRE + ATAC distribution)...")
    pooled = cre_rows + atac_rows
    _, meta = build_bins_table(
        pooled, args.len_bins, args.gc_bins, args.tssd_bins, args.atac_bins
    )
    len_edges  = meta["len_edges"]
    gc_edges   = meta["gc_edges"]
    tssd_edges = meta["tssd_log10_edges"]
    atac_edges = meta["atac_edges"]

    def assign_all(rows):
        """Re-assign bins to a subset of rows using the shared pooled edges."""
        out = []
        for r in rows:
            r = dict(r)  # shallow copy to avoid mutating pooled rows
            r["len_bin"]  = assign_bin(r["length"],                 len_edges)
            r["gc_bin"]   = assign_bin(r["gc"],                      gc_edges)
            r["tssd_bin"] = assign_bin(math.log10(r["tssd"] + 1.0), tssd_edges)
            r["atac_bin"] = assign_bin(r["atac"],                    atac_edges)
            out.append(r)
        return out

    cre_rows  = assign_all(cre_rows)
    atac_rows = assign_all(atac_rows)

    # Save covariate tables (useful for QC / verifying matching quality downstream)
    hdr = ["chrom", "start", "end", "id", "length", "gc", "tssd", "atac",
           "len_bin", "gc_bin", "tssd_bin", "atac_bin"]
    write_tsv(os.path.join(args.outdir, "CREs_covariates.tsv"),  cre_rows,  hdr)
    write_tsv(os.path.join(args.outdir, "ATAC_covariates.tsv"),  atac_rows, hdr)

    # Save bin summary — inspect for empty bins (CREs with no matching ATAC peers)
    def key(r):
        return (r["len_bin"], r["gc_bin"], r["tssd_bin"], r["atac_bin"])

    c    = Counter([key(r) for r in cre_rows])
    a    = Counter([key(r) for r in atac_rows])
    keys = sorted(set(list(c.keys()) + list(a.keys())))
    bin_rows = [
        {"len_bin": k[0], "gc_bin": k[1], "tssd_bin": k[2], "atac_bin": k[3],
         "n_cre": c.get(k, 0), "n_atac": a.get(k, 0)}
        for k in keys
    ]
    write_tsv(
        os.path.join(args.outdir, "bins_summary.tsv"), bin_rows,
        ["len_bin", "gc_bin", "tssd_bin", "atac_bin", "n_cre", "n_atac"]
    )

    # ── Step 5: Matching ──────────────────────────────────────────────────────
    print(f"Matching {len(cre_rows)} CREs to background (threads={args.threads}, "
          f"seed={args.seed}, n_per_cre={args.n_per_cre})...")

    # Pre-build bin-cell lookup for O(1) pool retrieval per CRE
    pool_by_bin = defaultdict(list)
    for arow in atac_rows:
        pool_by_bin[(arow["len_bin"], arow["gc_bin"],
                     arow["tssd_bin"], arow["atac_bin"])].append(arow)

    scaling = prepare_knn_scaling(atac_rows)  # for kNN fallback

    tasks = list(enumerate(cre_rows))

    if args.threads <= 1:
        # Sequential mode: initialise globals directly (no subprocess overhead)
        _init_worker(pool_by_bin, atac_rows, scaling, args.seed, args.n_per_cre)
        results = [_match_one(t) for t in tasks]
    else:
        # Parallel mode: each worker process receives globals via initializer
        with Pool(
            processes=args.threads,
            initializer=_init_worker,
            initargs=(pool_by_bin, atac_rows, scaling, args.seed, args.n_per_cre),
        ) as P:
            results = P.map(_match_one, tasks)

    # ── Collect results ───────────────────────────────────────────────────────
    matches_1to1 = []
    matches_1to5 = []
    for one, many in results:
        if one[1] is not None:
            matches_1to1.append(one)
        matches_1to5.extend(many)

    n_unmatched = len(cre_rows) - len(matches_1to1)
    print(f"  Matched 1:1 = {len(matches_1to1)}  |  Unmatched = {n_unmatched}  "
          f"|  Matched 1:N pairs = {len(matches_1to5)}")

    # ── Step 6: Write outputs ─────────────────────────────────────────────────
    with open(os.path.join(args.outdir, "matches_1to1.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["CRE_ID", "ATAC_ID"])
        w.writerows(matches_1to1)

    with open(os.path.join(args.outdir, "matches_1to5.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["CRE_ID", "ATAC_ID"])
        w.writerows(matches_1to5)

    # Export matched ATAC intervals as BED files for use with bedtools / bigWigAverageOverBed
    atac_idx = {a["id"]: a for a in atac_rows}

    ids1 = sorted({b for _, b in matches_1to1})
    ids5 = sorted({b for _, b in matches_1to5})

    def write_bed(ids, outp):
        """Write a set of ATAC peak IDs to a 4-column BED file."""
        with open(outp, "w") as f:
            for aid in ids:
                a = atac_idx.get(aid)
                if a is None:
                    continue
                f.write(f"{a['chrom']}\t{a['start']}\t{a['end']}\t{a['id']}\n")

    write_bed(ids1, os.path.join(args.outdir, "ATAC_matched_1to1.bed"))
    write_bed(ids5, os.path.join(args.outdir, "ATAC_matched_1to5.bed"))

    # Save bin edges as JSON for exact reproducibility across runs and datasets
    with open(os.path.join(args.outdir, "bin_edges.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print("Done. Outputs written to:", args.outdir)


if __name__ == "__main__":
    main()
