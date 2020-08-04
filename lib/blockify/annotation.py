import numpy as np
import pandas as pd
from pybedtools import BedTool
from . import segmentation
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import sys
from . import utilities
import warnings

# Suppress certain warnings
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=ResourceWarning)

# Sub-list collapsing method from https://stackoverflow.com/a/952952
CORRECTION_METHODS = [method for sublist in multitest._alias_list for method in sublist]


def validateAnnotationArguments(
    input_file,
    regions_bed,
    background_file,
    measure,
    alpha,
    correction,
    p_value,
    distance,
    min_size,
    max_size,
    pseudocount,
):
    # Check that the input file is sorted
    assert utilities.isSortedBEDObject(input_file), "input file must be sorted"
    # Check that regions_bed is sorted
    assert utilities.isSortedBEDObject(regions_bed), "regions BED file must be sorted"
    # Check that background_file is sorted
    assert utilities.isSortedBEDObject(
        background_file
    ), "background file must be sorted"
    # Check that measure is a valid parameter
    assert measure in ["enrichment", "depletion"], "measurement must be either 'enrichment' or 'depletion'"
    # If regions has been supplied, check that they are also sorted
    if regions_bed:
        assert utilities.isSortedBEDFile(regions_bed), "regions BED file must be sorted"
    # Check that alpha or pValueCutoff are valid
    if alpha:
        assert 0 <= alpha <= 1, "--alpha should be between 0 and 1, inclusive"
        # If alpha was provided, correction should have been provided
        # correction should be any valid method in statsmodels.stats.multitest
        assert (
            correction in CORRECTION_METHODS
        ), "invalid multiple hypothesis correction method"
    if p_value:
        assert 0 <= p_value <= 1, "--pValueCutoff should be non-negative"
    # Check that the other integers, if provided, are valid
    if distance:
        assert distance >= 0, "--distance should be a non-negative integer"
    if min_size:
        assert min_size >= 0, "--min should be a non-negative integer"
    if max_size:
        assert max_size >= 0, "--max should be a non-negative integer"
    # Check pseudocount is a positive real number
    assert pseudocount >= 0, "--pseudocount should be non-negative"


def tighten(data):
    # data is a BedTool; return value is also BedTool
    # Calculate new boundaries based on first and last event locations in the region
    df = data.to_dataframe()
    df = df.iloc[:, :6]
    df = df.rename(
        index=str,
        columns={"name": "TTAA_chrom", "score": "TTAA_start", "strand": "TTAA_end"},
    )
    df = df.sort_values(
        ["chrom", "start", "end", "TTAA_chrom", "TTAA_start", "TTAA_end"]
    )
    groups = df.groupby(["chrom", "start", "end"])
    first = groups.nth(0)["TTAA_start"]
    last = groups.nth(-1)["TTAA_end"]
    joined = pd.concat([first, last], axis=1).reset_index()
    refined = joined[["chrom", "TTAA_start", "TTAA_end"]]
    return BedTool.from_dataframe(refined)


def parcelConsecutiveBlocks(df):
    # df is a pandas DataFrame; return value is a list of
    # DataFrames, each of which contains consecutive blocks
    outlist = []
    chroms = np.unique(df.chrom)
    for chrom in chroms:
        chrom_rows = df[df.chrom == chrom].reset_index()
        for i, row in chrom_rows.iterrows():
            i = int(i)
            if i == 0:
                # Start of a new chromosome
                outlist.append([row])
            else:
                # Check if this block is adjacent to the previous block
                if chrom_rows.iloc[i].start == chrom_rows.iloc[i - 1].end:
                    # If so, append to the previous list
                    outlist[-1].append(row)
                else:
                    # If not, create a new list
                    outlist.append([row])
    return [pd.DataFrame(_) for _ in outlist]


def getPeakSummits(df, metric="pValue"):
    # df is a pandas DataFrame; return value is also a DataFrame
    assert metric in ["pValue", "density"]
    if metric == "pValue":
        column = "negLog10corrected"
    elif metric == "density":
        column = "Net_density"
    else:
        # Something unexpected happened
        print("Unexpected error:", sys.exc_info()[0])
        raise
    return pd.concat([_[_[column] == _[column].max()] for _ in parcelConsecutiveBlocks(df)]).reset_index().drop(columns=["level_0", "index"])


def sizeFilter(bed, min_size, max_size):
    # For size filter, convert to DataFrame, filter, and go back to BED
    df = bed.to_dataframe()
    df["size"] = df["end"] - df["start"]
    df = df[df["size"] <= max_size]
    df = df[df["size"] >= min_size]
    # Also set columns for name, score, and strand
    df["name"] = "peak_" + df.index.astype(str)
    df["score"] = 1
    df["strand"] = "."
    return BedTool.from_dataframe(df[["chrom", "start", "end", "name", "score", "strand"]])


def annotate(
    input_file,
    regions_bed,
    background_file,
    measure="enrichment",
    intermediate=None,
    alpha=None,
    p_value=None,
    correction=None,
    distance=None,
    min_size=None,
    max_size=None,
    pseudocount=1,
    tight=False,
    summit=False,
):
    # input_file, regions, and background_file are BedTool objects
    # Validate annotation arguments
    validateAnnotationArguments(
        input_file,
        regions_bed,
        background_file,
        measure,
        alpha,
        correction,
        p_value,
        distance,
        min_size,
        max_size,
        pseudocount,
    )
    # Calculate scaling factor
    scalingFactor = len(input_file.to_dataframe()) / len(background_file.to_dataframe())

    # Pull region edges to the nearest event in input, if specified
    if tight:
        data = regions_bed.intersect(input_file, wa=True, wb=True, sorted=True)
        regions_bed = tighten(data)

    # Intersect regions with the input file
    data = regions_bed.intersect(input_file, c=True, sorted=True).intersect(background_file, c=True, sorted=True)
    # Convert to DataFrame
    df = data.to_dataframe()
    df = df.rename(index=str, columns={"name": "Input", "score": "Background"})

    # Calculate the normalized background (Normed_bg) number of events
    # by multiplying Background by scalingFactor. Then add the pseudocount,
    # in the case of Normed_bg; and the floor of the pseudocount to Input.
    # This preserves log-fold change if Background is 0, and keeps the value
    # added to Input an integer (if pseudocount is a float).
    df["Input"] += np.floor(pseudocount)
    df["Normed_bg"] = df["Background"] * scalingFactor + pseudocount

    # Calculate density of insertions in each block
    df["Net_density"] = (df["Input"] - df["Normed_bg"]) / (df["end"] - df["start"])

    if measure == "enrichment":
        # Calculate the one-tail Poisson p-value of observing Input or more number of events
        # given a lambda of Normed_bg. Use the survival function (sf), which is 1 - CDF.
        # Need to specify Input - 1 because for discrete distributions, 1 - cdf(x) is p(X ≥ x + 1);
        # sf(x - 1) is thus p(X ≥ x).
        df["pValue"] = stats.poisson.sf(df["Input"] - 1, df["Normed_bg"])
    elif measure == "depletion":
        # Calculate the one-tail Poisson p-value of observing Input or fewer number of events
        # given a lambda of Normed_bg. Use the cumulative distribution function (cdf).
        # cdf(x) is  p(X ≤ x).
        df["pValue"] = stats.poisson.cdf(df["Input"], df["Normed_bg"])
    else:
        # Something unexpected happened
        print("Unexpected error:", sys.exc_info()[0])
        raise
    # Replace 0's with 1/FLOAT_MAX to obtain finite -log10(pValue)
    df.replace(to_replace=0, value=1 / utilities.FLOAT_MAX, inplace=True)
    df["negLog10pValue"] = -np.log10(df["pValue"])

    # If p_value has been provided, filter by it; else, perform multiple hypothesis correction
    if p_value:
        df["rejected"] = df["pValue"] <= p_value
        out_df = df[df["pValue"] <= p_value]
    else:
        df["rejected"], df[
            "corrected_pValue"
        ], alphacSidak, alphacBonf = multitest.multipletests(
            df["pValue"],
            alpha=alpha,
            method=correction,
            is_sorted=False,
            returnsorted=False,
        )
        df["negLog10corrected"] = -np.log10(df["corrected_pValue"])
        out_df = df[df["rejected"]]

    # Return peak summits, if specified
    if summit:
        out_df = getPeakSummits(out_df)

    # Convert out_df to out_bed
    out_bed = BedTool.from_dataframe(out_df)

    # Merge peaks within distance, if specified
    if distance is not None:
        out_bed = out_bed.merge(d=distance)

    # Filter by minimum and/or maximum size, if specified
    if min_size:
        minSize = min_size
    else:
        minSize = 0

    if max_size:
        maxSize = max_size
    else:
        maxSize = np.inf
    out_bed = sizeFilter(out_bed, minSize, maxSize)

    # Return out_bed and intermediate file, if any
    if intermediate:
        return out_bed, df
    else:
        return out_bed, None


def annotate_from_command_line(args):
    input_file = BedTool(args.input)
    # If regions has been supplied, use it;
    if args.regions:
        regions_bed = BedTool(args.regions)
    # otherwise, segment the file
    else:
        region_segmentation = segmentation.segment(
            input_file, args.method, p0=args.p0, prior=args.prior
        )
        regions_bed = BedTool.from_dataframe(region_segmentation.df)
    background_file = BedTool(args.background)

    return annotate(
        input_file,
        regions_bed,
        background_file,
        measure=args.measure,
        intermediate=args.intermediate,
        alpha=args.alpha,
        p_value=args.pValueCutoff,
        correction=args.correction,
        distance=args.distance,
        min_size=args.min,
        max_size=args.max,
        pseudocount=args.pseudocount,
        tight=args.tight,
    )
