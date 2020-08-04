from blockify.parsers import blockify_parser
import blockify.downsampling as downsampling
import os
import pandas as pd
import sys
import unittest


# Enable warnings (see https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("default")  # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default"  # Also affect subprocesses


class TestDownsampling(unittest.TestCase):
    def test_yeast_normal(self):
        args = blockify_parser.parse_args(
            [
                "downsample",
                "--input",
                "tests/data/S288C_dSIR4.qbed",
                "-n",
                "10000",
                "--seed",
                "0",
            ]
        )
        result = downsampling.downsample_from_command_line(args)
        self.assertEqual(result.to_dataframe()["start"][0], 2091)

    def test_yeast_naive(self):
        args = blockify_parser.parse_args(
            [
                "downsample",
                "--input",
                "tests/data/S288C_dSIR4.qbed",
                "-n",
                "10000",
                "--seed",
                "0",
                "--naive",
            ]
        )
        result = downsampling.downsample_from_command_line(args)
        self.assertEqual(result.to_dataframe()["start"][0], 2118)


class TestDownsamplingAPI(unittest.TestCase):
    def test_downsampling_API_normal(self):
        result = downsampling.downsample(
            pd.read_table("tests/data/S288C_dSIR4.qbed", header=None),
            10000,
            seed=0,
            naive=False
        )
        self.assertEqual(result.to_dataframe()["start"][0], 2091)

    def test_downsampling_API_naive(self):
        result = downsampling.downsample(
            pd.read_table("tests/data/S288C_dSIR4.qbed", header=None),
            10000,
            seed=0,
            naive=True
        )
        self.assertEqual(result.to_dataframe()["start"][0], 2118)


if __name__ == "__main__":
    unittest.main()
