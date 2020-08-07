import numpy as np
import pandas as pd
from pybedtools import BedTool


def downsample(input_file, n, seed=None, naive=False):
    """Core downsampling method

    Parameters
    ----------
    input_file: ``pandas`` DataFrame
        Input data (e.g. BED, qBED, CCF) as a ``pandas`` DataFrame
    n: int
        Number of entries to sample
    seed: int
        Seed for random number generator
    naive: bool
        Choose whether to sample each entry with equal probability (True) or weighted by the value in the fourth column (if supplied)

    Returns
    -------
    downsampled_file: BedTool object
        Input file after downsampling
    """

    # input_file is a pandas DataFrame
    total = np.sum(input_file[3])  # total
    if naive:
        # Select rows with equal likelihood
        p = None
    else:
        # Select rows with proportional weights
        p = input_file[3] / total
    if seed is not None:  # set random seed if provided
        np.random.seed(seed)
    # Sample rows
    indexes = np.random.choice(np.arange(len(input_file)), size=n, replace=False, p=p)
    indexes.sort()
    downsampled_file = input_file.iloc[indexes]
    return BedTool.from_dataframe(downsampled_file)


# Downsample a qBED file from the command line
# Thin wrapper for downsample()
def downsample_from_command_line(args):
    """Wrapper function for the command line function ``blockify downsample``

    Parameters
    ----------
    args: ``argparse.Namespace`` object
        Input from command line

    Returns
    -------
    downsampled_file: BedTool
        Downsampled command line data
    """

    input_file = pd.read_table(args.input, header=None)
    # Segment the input file
    return downsample(input_file, args.number, seed=args.seed, naive=args.naive)
