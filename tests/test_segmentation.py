from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_SEGMENTATION_P0
import blockify.segmentation as segmentation
from pybedtools import BedTool
import sys
import unittest

# Enable warnings (see https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestSegmentation(unittest.TestCase):
    def test_PBase(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "test.bed"
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 30172)
        self.assertEqual(result.total_fitness, 10271168.936047826)

    def test_SP1_PBase(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-SP1-PBase_sort_canon.ccf",
                "test.bed",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 21375)
        self.assertEqual(result.total_fitness, 2987508.577343301)

    def test_HyPBase(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-HyPBase_sort_canon.ccf",
                "test.bed",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 106211)
        self.assertEqual(result.total_fitness, 31079910.017576993)

    def test_SP1_HyPBase(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-SP1-HyPBase_sort_canon.ccf",
                "test.bed",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 68935)
        self.assertEqual(result.total_fitness, 11410811.054755056)

    def test_uniformity(self):
        # Test that the algorithm segments a truly uniform dataset into a single block
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/uniform_test.ccf",
                "test.bed"
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 1)

    def test_equivalence(self):
        # Test that segmenting using OP or PELT yield the same segmentation
        args_OP = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--method",
                "OP",
                "test.bed",
            ]
        )
        result_OP = segmentation.segment_from_command_line(args_OP)
        args_PELT = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--method",
                "PELT",
                "test.bed",
            ]
        )
        result_PELT = segmentation.segment_from_command_line(args_PELT)
        self.assertEqual(result_OP.total_blocks, result_PELT.total_blocks)
        # The fitness calculation is identical between the algorithms, except for a sign change
        self.assertEqual(result_OP.total_fitness, -result_PELT.total_fitness)


class TestSegmentationParameters(unittest.TestCase):
    def test_p0(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--p0",
                "2",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            segmentation.segment_from_command_line(args)

    def test_prior(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--prior",
                "-1",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            segmentation.segment_from_command_line(args)


class TestSegmentationAPI(unittest.TestCase):
    def test_segmentation_API(self):
        result = segmentation.segment(
            BedTool("tests/data/HCT116-PBase_sort_canon.ccf"),
            "PELT",
            p0=DEFAULT_SEGMENTATION_P0,
            prior=None
        )
        self.assertEqual(result.total_blocks, 30172)
        self.assertEqual(result.total_fitness, 10271168.936047826)


if __name__ == "__main__":
    unittest.main()
