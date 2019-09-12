from blockify import annotation
from blockify import segmentation
from blockify import normalization
from blockify import downsampling
from blockify.parsers import blockify_parser
import sys

# Disable warnings for command line use (https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def run():
    args = blockify_parser.parse_args()
    if args.command == "segment":
        result = segmentation.segment_from_command_line(args)
        # result is a SegmentationResult object
        result.df.to_csv(args.output, sep="\t", index=None, header=None)
    elif args.command == "normalize":
        result = normalization.normalize_from_command_line(args)
        # result is a BedTool object
        result.saveas(args.output)
    elif args.command == "call":
        result, intermediate = annotation.annotate_from_command_line(args)
        # result is a BedTool object
        result.saveas(args.output)
        if intermediate:
            intermediate.to_csv(args.intermediate)
    elif args.command == "downsample":
        result = downsampling.downsample_from_command_line(args)
        # result is a BedTool object
        result.saveas(args.output)
    else:
        raise Exception


if __name__ == "__main__":
    run()
