from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_NORMALIZATION_LIBRARY_FACTOR, DEFAULT_NORMALIZATION_LENGTH_FACTOR
import blockify.normalization as normalization
from pybedtools import BedTool
import sys
import unittest


# Enable warnings
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestNormalization(unittest.TestCase):
    def test_PBase(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--output",
                "test.bedgraph",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        self.assertEqual(result.to_dataframe()["name"][0], 44.04857703372938)

    def test_SP1_PBase(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT116-SP1-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-SP1-PBase_sort_canon.blocks",
                "--output",
                "test.bedgraph",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        self.assertEqual(result.to_dataframe()["name"][0], 29.226377780159183)

    def test_HyPBase(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT116-HyPBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-HyPBase_sort_canon.blocks",
                "--output",
                "test.bedgraph",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        self.assertEqual(result.to_dataframe()["name"][0], 3.2922055992793187)

    def test_SP1_HyPBase(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT116-SP1-HyPBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-SP1-HyPBase_sort_canon.blocks",
                "--output",
                "test.bedgraph",
            ]
        )
        result = normalization.normalize_from_command_line(args)
        self.assertEqual(result.to_dataframe()["name"][0], 23.6461239322913)


class TestNormalizationParameters(unittest.TestCase):
    def test_libraryFactor(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--libraryFactor",
                "0",
                "--output",
                "test.bedgraph",
            ]
        )
        with self.assertRaises(AssertionError):
            normalization.normalize_from_command_line(args)

    def test_lengthFactor(self):
        args = blockify_parser.parse_args(
            [
                "normalize",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--lengthFactor",
                "0",
                "--output",
                "test.bedgraph",
            ]
        )
        with self.assertRaises(AssertionError):
            normalization.normalize_from_command_line(args)


class TestNormalizationAPI(unittest.TestCase):
    def test_normalization_API(self):
        result = normalization.normalize(
            BedTool("tests/data/HCT116-PBase_sort_canon.ccf"),
            BedTool("tests/data/HCT116-PBase_sort_canon.blocks"),
            DEFAULT_NORMALIZATION_LIBRARY_FACTOR,
            DEFAULT_NORMALIZATION_LENGTH_FACTOR
        )
        self.assertEqual(result.to_dataframe()["name"][0], 44.04857703372938)


if __name__ == "__main__":
    unittest.main()
