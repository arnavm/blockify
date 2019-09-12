import argparse
import sys

# Initialize defaults
DEFAULT_SEGMENTATION_P0 = 0.05
DEFAULT_NORMALIZATION_LIBRARY_FACTOR = 10 ** 6
DEFAULT_NORMALIZATION_LENGTH_FACTOR = None
DEFAULT_MULTIPLE_HYPOTHESIS_CORRECTION = "bonferroni"
DEFAULT_PSEUDOCOUNT = 1

# Top-most parser
blockify_parser = argparse.ArgumentParser(
    description="Genomic peak caller for one-dimensional data"
)
# Sub-parsers
subcommands = blockify_parser.add_subparsers(help="Subcommands", dest="command")
# Parent parser (shared by sub-parsers)
input_parser = argparse.ArgumentParser(add_help=False)
input_parser.add_argument("-i", "--input", required=True, help="Input .ccf file")
regions_parser = argparse.ArgumentParser(add_help=False)
regions_parser.add_argument(
    "-r",
    "--regions",
    required=False,
    help="Regions over which to normalize event counts; should be supplied as a BED file. If not provided, the input file will be segmented using Bayesian blocks.",
)

# Segmentation sub-command
segment = subcommands.add_parser(
    name="segment",
    description="Segment a .ccf file using Bayesian blocks",
    help="Segment a .ccf file using Bayesian blocks",
    parents=[input_parser],
)
segment.add_argument("-o", "--output", required=True, help="Output file (BED format)")
prior_group = segment.add_mutually_exclusive_group(required=False)
prior_group.add_argument(
    "--prior", type=float, help="Explicit prior on the number of blocks"
)
prior_group.add_argument(
    "--p0",
    type=float,
    default=DEFAULT_SEGMENTATION_P0,
    help="Empirical prior based on a specified false-positive rate (between 0 and 1)",
)
# prior_group.add_argument("--gamma",
#                          type=float,
#                          help="False positive rate of default prior rate")
segment.add_argument(
    "--method",
    choices=["OP", "PELT"],
    default="PELT",
    help="Segment using the optimal partitioning (OP) or pruned exact linear time (PELT) algorithm",
)
# segment.add_argument("-t",
#                      "--time",
#                      action="store_true",
#                      default=False,
#                      help="Time the execution of the algorithm")

# Normalize sub-command
normalize = subcommands.add_parser(
    name="normalize",
    description="Calculate normalized rates of events in a .ccf file",
    help="Calculate normalized rates of events in a .ccf file",
    parents=[segment, regions_parser],
    conflict_handler="resolve",
)
normalize.add_argument(
    "-o", "--output", required=True, help="Output file (bedGraph format)"
)
normalize.add_argument(
    "-k",
    "--libraryFactor",
    type=float,
    default=DEFAULT_NORMALIZATION_LIBRARY_FACTOR,
    help="Normalization factor for library size",
)
normalize.add_argument(
    "-l",
    "--lengthFactor",
    type=float,
    default=DEFAULT_NORMALIZATION_LENGTH_FACTOR,
    help="Normalization factor for the length of regions; used to calculate scaled rates of events per interval",
)

# Peak calling (annotation) sub-command
annotate = subcommands.add_parser(
    name="call",
    description="Call peaks in a .ccf file",
    help="Call peaks in a .ccf file",
    parents=[segment, regions_parser],
    conflict_handler="resolve",
)
annotate.add_argument(
    "-bg", "--background", type=str, help="Background .ccf file", required=True
)
annotate.add_argument(
    "--intermediate",
    type=str,
    required=False,
    help="Intermediate file to write verbose output (CSV format)",
)
alpha_group = annotate.add_mutually_exclusive_group(required=True)
alpha_group.add_argument(
    "-a",
    "--alpha",
    type=float,
    help="Alpha for multiple hypothesis correction (must be between 0 and 1)",
)
alpha_group.add_argument(
    "-p",
    "--pValueCutoff",
    type=float,
    help="p-value cutoff (NOTE: This is a straight cutoff and will not take into account multiple hypothesis correction!)",
)
annotate.add_argument(
    "--correction",
    type=str,
    required="-a" in sys.argv or "--alpha" in sys.argv,
    default=DEFAULT_MULTIPLE_HYPOTHESIS_CORRECTION,
    help="If alpha provided, need to specificity method of multiple hypothesis correction. See statsmodels.stats.multitest for a complete list of choices",
)
annotate.add_argument(
    "-d",
    "--distance",
    type=int,
    required=False,
    help="Merge features closer than this distance (bp)",
)
annotate.add_argument(
    "--min", type=int, required=False, help="Report peaks larger than this cutoff (bp)"
)
annotate.add_argument(
    "--max", type=int, required=False, help="Report peaks smaller than this cutoff (bp)"
)
annotate.add_argument("-t", "--tight", action="store_true", default=False)
annotate.add_argument(
    "-c",
    "--pseudocount",
    type=float,
    default=DEFAULT_PSEUDOCOUNT,
    help="Pseudocount for background regions (default: %(default)s)",
)

# Downsample sub-command
downsample = subcommands.add_parser(
    name="downsample",
    description="Downsample a .ccf file in proportion to the value column",
    help="Downsample a .ccf file in proportion to the value column",
    parents=[input_parser],
)
downsample.add_argument(
    "-n",
    "--number",
    type=int,
    help="Number of events to downsample (cannot exceed length of input file)",
    required=True,
)
downsample.add_argument("-s", "--seed", type=int, help="Random seed")
downsample.add_argument(
    "--naive", action="store_true", help="Sample every row with equal likelihood"
)
downsample.add_argument(
    "-o", "--output", required=True, help="Output file (CCF format)"
)
