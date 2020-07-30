from pybedtools import BedTool
from . import segmentation
from . import utilities
import warnings

# Generates a normalized bedgraph of event rate for BED-formatted genomic data.
# The output file can be visualized on genome browsers.

# Suppress certain warnings
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=ResourceWarning)


def validateNormalizationArguments(input_file, regions_bed, libraryFactor, lengthFactor):
    # Check that the input file is sorted
    assert utilities.isSortedBEDObject(input_file), "input file must be sorted"
    # Check that regions is also sorted
    assert utilities.isSortedBEDObject(regions_bed), "regions BED file must be sorted"
    # Check libraryFactor
    assert libraryFactor > 0, "--libraryFactor should be a positive number"
    # Check lengthFactor
    if lengthFactor is not None:
        assert lengthFactor > 0, "--lengthFactor should be a positive number"


def normalize(input_file, regions_bed, libraryFactor, lengthFactor):
    # input_file and regions_bed are BedTool objects
    # Calculate library scaling constant, which is the total number
    # Validate normalization arguments
    validateNormalizationArguments(input_file, regions_bed, libraryFactor, lengthFactor)
    # of events in input BED divided by the library factor
    library_scaling_constant = len(input_file.to_dataframe()) / libraryFactor

    # For each interval in regions, count the number of events;
    # normalize the count by library_scaling_constant.
    # The last call to .iloc should be able to accomodate region BED files with arbitrary numbers of fields
    intersect_df = regions_bed.intersect(input_file, c=True, sorted=True).to_dataframe().iloc[:, [0, 1, 2, -1]]
    intersect_df.columns = ["chrom", "start", "end", "rawCount"]
    intersect_df["normCount"] = intersect_df["rawCount"] / library_scaling_constant

    # If lengthFactor has been provided, calculate normalized rates of events
    if lengthFactor:
        intersect_df["normRate"] = intersect_df["normCount"] / (
            (intersect_df["end"] - intersect_df["start"]) / lengthFactor
        )

    # Return a BedTool object
    if lengthFactor:
        return BedTool.from_dataframe(
            intersect_df[["chrom", "start", "end", "normRate"]]
        )
    else:
        return BedTool.from_dataframe(
            intersect_df[["chrom", "start", "end", "normCount"]]
        )


# Normalize the rate of events in a file over a set of regions
# Thin wrapper for normalize()
def normalize_from_command_line(args):
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
    return normalize(input_file, regions_bed, args.libraryFactor, args.lengthFactor)
