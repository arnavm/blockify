from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_NORMALIZATION_LIBRARY_FACTOR, DEFAULT_NORMALIZATION_LENGTH_FACTOR
import blockify.normalization as normalization
import gzip
import os
from pybedtools import BedTool
import sys
import unittest
from urllib.request import Request, urlopen


# Enable warnings
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestNormalization(unittest.TestCase):
    def test_yeast(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        self.assertEqual(result.to_dataframe()["name"][0], 861.7222956281956)

    def test_PBase(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urlopen(PBase_url) as response, open("tests/data/HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        PBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_PBase.blocks"
        with urlopen(Request(PBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_PBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT-116_PBase.ccf",
                "--regions",
                "tests/data/HCT-116_PBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("tests/data/HCT-116_PBase.ccf")
        os.remove("tests/data/HCT-116_PBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 44.04857703372938)

    def test_SP1_PBase(self):
        SP1_PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471637&format=file&file=GSM4471637%5FHCT%2D116%5FSP1%2DPBase%2Eccf%2Etxt%2Egz"
        with urlopen(SP1_PBase_url) as response, open("tests/data/HCT-116_SP1-PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        SP1_PBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_SP1-PBase.blocks"
        with urlopen(Request(SP1_PBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_SP1-PBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT-116_SP1-PBase.ccf",
                "--regions",
                "tests/data/HCT-116_SP1-PBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("tests/data/HCT-116_SP1-PBase.ccf")
        os.remove("tests/data/HCT-116_SP1-PBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 29.226377780159183)

    def test_HyPBase(self):
        HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471638&format=file&file=GSM4471638%5FHCT%2D116%5FHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(HyPBase_url) as response, open("tests/data/HCT-116_HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        HyPBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_HyPBase.blocks"
        with urlopen(Request(HyPBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_HyPBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT-116_HyPBase.ccf",
                "--regions",
                "tests/data/HCT-116_HyPBase.blocks",
            ]
        ) 
        result = normalization.normalize_from_command_line(args)
        os.remove("tests/data/HCT-116_HyPBase.ccf")
        os.remove("tests/data/HCT-116_HyPBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 3.2922055992793187)

    def test_SP1_HyPBase(self):
        SP1_HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471639&format=file&file=GSM4471639%5FHCT%2D116%5FSP1%2DHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(SP1_HyPBase_url) as response, open("tests/data/HCT-116_SP1-HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        SP1_HyPBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_SP1-HyPBase.blocks"
        with urlopen(Request(SP1_HyPBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_SP1-HyPBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT-116_SP1-HyPBase.ccf",
                "--regions",
                "tests/data/HCT-116_SP1-HyPBase.blocks",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        os.remove("tests/data/HCT-116_SP1-HyPBase.ccf")
        os.remove("tests/data/HCT-116_SP1-HyPBase.blocks")
        self.assertEqual(result.to_dataframe()["name"][0], 23.6461239322913)


class TestNormalizationParameters(unittest.TestCase):
    def test_libraryFactor(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--libraryFactor",
                "0",
            ]
        )
        with self.assertRaises(AssertionError):
            normalization.normalize_from_command_line(args)

    def test_lengthFactor(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--lengthFactor",
                "0",
            ]
        )
        with self.assertRaises(AssertionError):
            normalization.normalize_from_command_line(args)


class TestNormalizationAPI(unittest.TestCase):
    def test_normalization_API(self):
        result = normalization.normalize(
            BedTool("tests/data/S288C_CBF1.qbed"),
            BedTool("tests/data/S288C_CBF1.blocks"),
            DEFAULT_NORMALIZATION_LIBRARY_FACTOR,
            DEFAULT_NORMALIZATION_LENGTH_FACTOR
        )
        self.assertEqual(result.to_dataframe()["name"][0], 861.7222956281956)


if __name__ == "__main__":
    unittest.main()
