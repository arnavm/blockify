import numpy as np
import pandas as pd
from pybedtools import BedTool


def downsample(input_file, n, seed=None, naive=False):
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
    input_file = pd.read_table(args.input, header=None)
    # Segment the input file
    return downsample(input_file, args.number, seed=args.seed, naive=args.naive)
