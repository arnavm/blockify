from . import algorithms
from io import StringIO
import numpy as np
import pandas as pd
from pybedtools import BedTool
import sys
from . import utilities
import warnings

# Suppress certain warnings
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=ResourceWarning)


# A class to store Bayesian block segmentation, to facilitate
# passing segmentations between modules.
class SegmentationRecord(object):
    def __init__(self):
        # # File to be segmented
        # self.filename = None
        # p0 value that was used for segmentation
        self.p0 = None
        # dict to store priors, keyed by chromosome
        self.priors = {}
        # dict to store counts of how many blocks each chromosome was segmented into
        self.nblocks = {}
        # dict to store fitness values from segmentation of each chromosome
        self.fitness = {}
        # The actual segmentation itself, stored as a pandas DataFrame
        self.df = pd.DataFrame()
        # The segmentation as a BedTool object
        self.blocks = None
        # total_priors is sum of values in priors
        self.total_priors = 0
        # total_blocks is sum of values in nblocks
        self.total_blocks = 0
        # total_fitness is sum of values in fitness
        self.total_fitness = 0

    def finalize(self):
        self.total_priors = np.sum(list(self.priors.values()))
        self.total_blocks = np.sum(list(self.nblocks.values()))
        self.total_fitness = np.sum(list(self.fitness.values()))
        self.blocks = BedTool.from_dataframe(self.df)


def validateSegmentationArguments(input_file, p0, prior):
    # Check that input_file is sorted
    assert utilities.isSortedBEDObject(input_file), "input file must be sorted"
    # If prior has been provided, check that it is positive
    if prior:
        assert prior >= 0, "--prior should be non-negative"
    # If p0 has been provided, check that it is between 0 and 1
    if p0:
        assert 0 <= p0 <= 1, "--p0 should be between 0 and 1, inclusive"


# Convert a set of continuous Bayesian blocks to DataFrame format
def blocksToDF(chrom, ranges):
    output = ""
    # Chromosomes need at least two events at different positions
    # to be able to report a block. If output is empty,
    # return an emtpy DataFrame.
    for i in range(len(ranges) - 1):
        interval = "{}\t{}\t{}\n".format(chrom, ranges[i], ranges[i + 1])
        output += interval
    if output:
        return pd.read_csv(StringIO(output), header=None, sep="\t")
    else:
        return pd.DataFrame()


# Returns a SegmentationRecord object
def segment(input_file, method, p0=None, prior=None):
    # input_file is a BedTool object
    # Validate segmentation arguments
    validateSegmentationArguments(input_file, p0, prior)
    # Get the algorithm class, then instantiate with specified parameters
    algorithm = algorithms.ALGORITHM_DICT.get(method, method)
    if prior:
        alg = algorithm(ncp_prior=prior)
    elif p0 is not None:
        alg = algorithm(p0=p0)
    else:
        raise ValueError("--p0 or --prior argument are invalid")

    # Open input file as DataFrame
    input_df = input_file.to_dataframe()
    # Pre-process coordinates by taking the floor of the mean of the start and end values
    input_df["coordinate"] = (input_df["start"] + input_df["end"]) // 2
    # Get list of chromosomes specified in input_df
    chroms = utilities.getChromosomesInDF(input_df)
    n_chroms = len(chroms)
    # Now we are ready to segment. Instantiate a SegmentationRecord object
    segmentation = SegmentationRecord()
    # segmentation.filename = input_df

    if p0:
        segmentation.p0 = p0

    # Segment by chromosome
    for i, chrom in enumerate(chroms):
        # Print progress to stdout
        print("[{}/{}] Processing {}".format(i + 1, n_chroms, chrom), file=sys.stderr)
        # Get the block boundaries of the segmentation
        block_boundaries = alg.segment(
            input_df[input_df["chrom"] == chrom]["coordinate"]
        )
        # Conditional on if there are ≥ 1 block (≥ 2 boundaries)
        if len(block_boundaries) > 1:
            # Store the prior in the SegmentationRecord
            segmentation.priors[chrom] = alg.prior
            # Store the fitness of each segmentation
            segmentation.fitness[chrom] = alg.best_fitness
            # Convert block boundaries to DataFrame and append to the SegmentationRecord
            df = blocksToDF(chrom, block_boundaries.astype(int))
            segmentation.df = segmentation.df.append(df)
            # Store the number of blocks in the SegmentationRecord
            segmentation.nblocks[chrom] = len(df)
            print("--Found {} blocks".format(segmentation.nblocks[chrom]), file=sys.stderr)
        else:
            print("--Skipped, no blocks found", file=sys.stderr)
    segmentation.finalize()

    return segmentation


# Segment a genomic BED/qBED file from the command line
# Thin wrapper for segment()
def segment_from_command_line(args):
    input_file = BedTool(args.input)
    # Segment the input file
    return segment(input_file, args.method, p0=args.p0, prior=args.prior)
