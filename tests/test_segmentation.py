from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_SEGMENTATION_P0
import blockify.segmentation as segmentation
import gzip
import os
from pybedtools import BedTool
import sys
import unittest
from urllib.request import urlopen


# Enable warnings (see https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestSegmentation(unittest.TestCase):
    def test_yeast(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/S288C_CBF1.qbed",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 1408)
        self.assertEqual(result.total_fitness, 99777.41540439503)

    def test_PBase(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urlopen(PBase_url) as response, open("tests/data/HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT-116_PBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        os.remove("tests/data/HCT-116_PBase.ccf")
        self.assertEqual(result.total_blocks, 30172)
        self.assertEqual(result.total_fitness, 10271168.936047826)

    def test_SP1_PBase(self):
        SP1_PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471637&format=file&file=GSM4471637%5FHCT%2D116%5FSP1%2DPBase%2Eccf%2Etxt%2Egz"
        with urlopen(SP1_PBase_url) as response, open("tests/data/HCT-116_SP1-PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT-116_SP1-PBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        os.remove("tests/data/HCT-116_SP1-PBase.ccf")
        self.assertEqual(result.total_blocks, 21375)
        self.assertEqual(result.total_fitness, 2987508.577343301)

    def test_HyPBase(self):
        HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471638&format=file&file=GSM4471638%5FHCT%2D116%5FHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(HyPBase_url) as response, open("tests/data/HCT-116_HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT-116_HyPBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        os.remove("tests/data/HCT-116_HyPBase.ccf")
        self.assertEqual(result.total_blocks, 106211)
        self.assertEqual(result.total_fitness, 31079910.017576993)

    def test_SP1_HyPBase(self):
        SP1_HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471639&format=file&file=GSM4471639%5FHCT%2D116%5FSP1%2DHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(SP1_HyPBase_url) as response, open("tests/data/HCT-116_SP1-HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/HCT-116_SP1-HyPBase.ccf",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        os.remove("tests/data/HCT-116_SP1-HyPBase.ccf")
        self.assertEqual(result.total_blocks, 68935)
        self.assertEqual(result.total_fitness, 11410811.054755056)

    def test_uniformity(self):
        # Test that the algorithm segments a truly uniform dataset into a single block
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/test_uniform_sorted.qbed",
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
                "tests/data/S288C_CBF1.qbed",
                "--method",
                "OP",
            ]
        )
        result_OP = segmentation.segment_from_command_line(args_OP)
        args_PELT = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--method",
                "PELT",
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
                "tests/data/S288C_CBF1.qbed",
                "--p0",
                "2",
            ]
        )
        with self.assertRaises(AssertionError):
            segmentation.segment_from_command_line(args)

    def test_prior(self):
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--prior",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            segmentation.segment_from_command_line(args)


class TestSegmentationAPI(unittest.TestCase):
    def test_segmentation_API(self):
        result = segmentation.segment(
            BedTool("tests/data/S288C_CBF1.qbed"),
            "PELT",
            p0=DEFAULT_SEGMENTATION_P0,
            prior=None
        )
        self.assertEqual(result.total_blocks, 1408)
        self.assertEqual(result.total_fitness, 99777.41540439503)


if __name__ == "__main__":
    unittest.main()
