import numpy as np
import pandas as pd
from pybedtools import BedTool


def downsample(input_ccf, n, seed=None, naive=False):
    # input_ccf is a pandas DataFrame
    total = np.sum(input_ccf[3])  # total
    if naive:
        # Select rows with equal likelihood
        p = None
    else:
        # Select rows with proportional weights
        p = input_ccf[3] / total
    if seed is not None:  # set random seed if provided
        np.random.seed(seed)
    # Sample rows
    indexes = np.random.choice(np.arange(len(input_ccf)), size=n, replace=False, p=p)
    indexes.sort()
    downsampled_ccf = input_ccf.iloc[indexes]
    return BedTool.from_dataframe(downsampled_ccf)


# Downsample a CCF file from the command line
# Thin wrapper for downsample()
def downsample_from_command_line(args):
    input_ccf = pd.read_table(args.input, header=None)
    # Segment the input file
    return downsample(input_ccf, args.number, seed=args.seed, naive=args.naive)
