from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_NORMALIZATION_LIBRARY_FACTOR, DEFAULT_NORMALIZATION_LENGTH_FACTOR
import blockify.normalization as normalization
import blockify.segmentation as segmentation
import gzip
from pybedtools import BedTool
import sys
import unittest
import urllib.request


# Enable warnings
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestNormalization(unittest.TestCase):
    def test_PBase(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_PBase.blocks", sep='\t', index=None, header=None)
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "HCT-116_PBase.ccf",
                "--regions",
                "HCT-116_PBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")
        os.remove("HCT-116_PBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 44.04857703372938)

    def test_SP1_PBase(self):
        SP1_PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471637&format=file&file=GSM4471637%5FHCT%2D116%5FSP1%2DPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(SP1_PBase_url) as response, open("HCT-116_SP1-PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_SP1-PBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_SP1-PBase.blocks", sep='\t', index=None, header=None)
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "HCT-116_SP1-PBase.ccf",
                "--regions",
                "HCT-116_SP1-PBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("HCT-116_SP1-PBase.ccf")
        os.remove("HCT-116_SP1-PBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 29.226377780159183)

    def test_HyPBase(self):
        HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471638&format=file&file=GSM4471638%5FHCT%2D116%5FHyPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(HyPBase_url) as response, open("HCT-116_HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_HyPBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_HyPBase.blocks", sep='\t', index=None, header=None)
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "HCT-116_HyPBase.ccf",
                "--regions",
                "HCT-116_HyPBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("HCT-116_HyPBase.ccf")
        os.remove("HCT-116_HyPBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 3.2922055992793187)

    def test_SP1_HyPBase(self):
        SP1_HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471639&format=file&file=GSM4471639%5FHCT%2D116%5FSP1%2DHyPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(SP1_HyPBase_url) as response, open("HCT-116_SP1-HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_SP1-HyPBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_SP1-HyPBase.blocks", sep='\t', index=None, header=None)
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "HCT-116_SP1-HyPBase.ccf",
                "--regions",
                "HCT-116_SP1-HyPBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("HCT-116_SP1-HyPBase.ccf")
        os.remove("HCT-116_SP1-HyPBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 23.6461239322913)


class TestNormalizationParameters(unittest.TestCase):
    def test_libraryFactor(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--output",
                "HCT-116_PBase.blocks",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_PBase.blocks", sep='\t', index=None, header=None)
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "HCT-116_PBase.ccf",
                "--regions",
                "HCT-116_PBase.blocks",
                "--libraryFactor",
                "0",
            ]
        )
        with self.assertRaises(AssertionError):
            normalization.normalize_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")
        os.remove("HCT-116_PBase.blocks")

    def test_lengthFactor(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--output",
                "HCT-116_PBase.blocks",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_PBase.blocks", sep='\t', index=None, header=None)
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "HCT-116_PBase.ccf",
                "--regions",
                "HCT-116_PBase.blocks",
                "--lengthFactor",
                "0",
            ]
        )
        with self.assertRaises(AssertionError):
            normalization.normalize_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")
        os.remove("HCT-116_PBase.blocks")


class TestNormalizationAPI(unittest.TestCase):
    def test_normalization_API(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--output",
                "HCT-116_PBase.blocks",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        result.df.to_csv("HCT-116_PBase.blocks", sep='\t', index=None, header=None)
        result = normalization.normalize(
            BedTool("HCT-116_PBase.ccf"),
            BedTool("HCT-116_PBase.blocks"),
            DEFAULT_NORMALIZATION_LIBRARY_FACTOR,
            DEFAULT_NORMALIZATION_LENGTH_FACTOR
        )
        os.remove("HCT-116_PBase.ccf")
        os.remove("HCT-116_PBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 44.04857703372938)


if __name__ == "__main__":
    unittest.main()
