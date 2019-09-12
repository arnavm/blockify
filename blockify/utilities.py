from functools import reduce
import numpy as np
from pybedtools import BedTool

# numpy log10 float max
FLOAT_MAX = np.finfo(np.float64).max
LOG10_FLOAT_MAX = np.log10(FLOAT_MAX)


# Fast method for getting number of lines in a file
# For BED files, much faster than calling len() on file
# From https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def getChromosomesInDF(df):
    return reduce(lambda l, x: l if x in l else l + [x], df["chrom"], [])


def isSortedBEDObject(bed_object):
    # Convert BedTool object to pandas DataFrame
    df = bed_object.to_dataframe()
    # First, check that chrom is in sorted order
    if df["chrom"].is_monotonic:
        # If so, check that the start coordinates are in order
        chroms = getChromosomesInDF(df)
        for c in chroms:
            if not df[df["chrom"] == c]["start"].is_monotonic:
                return False
        return True
    return False


def isSortedBEDFile(bed_file_path):
    return isSortedBEDObject(BedTool(bed_file_path))
