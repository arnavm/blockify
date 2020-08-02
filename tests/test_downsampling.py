from blockify.parsers import blockify_parser
import blockify.downsampling as downsampling
import gzip
import os
import pandas as pd
import sys
import unittest
import urllib.request


# Enable warnings (see https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("default")  # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default"  # Also affect subprocesses


class TestDownsampling(unittest.TestCase):
    def test_PBase_normal(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "downsample",
                "--input",
                "HCT-116_PBase.ccf",
                "-n",
                "10000",
                "--seed",
                "0",
            ]
        )
        result = downsampling.downsample_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")
        self.assertEqual(result.to_dataframe()["start"][0], 900198)

    def test_PBase_naive(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "downsample",
                "--input",
                "HCT-116_PBase.ccf",
                "-n",
                "10000",
                "--seed",
                "0",
                "--naive",
            ]
        )
        result = downsampling.downsample_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")
        self.assertEqual(result.to_dataframe()["start"][0], 845648)


class TestDownsamplingAPI(unittest.TestCase):
    def test_downsampling_API_normal(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        result = downsampling.downsample(
            pd.read_table("HCT-116_PBase.ccf", header=None),
            10000,
            seed=0,
            naive=False
        )
        self.assertEqual(result.to_dataframe()["start"][0], 900198)
        os.remove("HCT-116_PBase.ccf")

    def test_downsampling_API_naive(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        result = downsampling.downsample(
            pd.read_table("HCT-116_PBase.ccf", header=None),
            10000,
            seed=0,
            naive=True
        )
        self.assertEqual(result.to_dataframe()["start"][0], 845648)
        os.remove("HCT-116_PBase.ccf")


if __name__ == "__main__":
    unittest.main()
