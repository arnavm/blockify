from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_SEGMENTATION_P0
import blockify.segmentation as segmentation
import gzip
import os
from pybedtools import BedTool
import re
import sys
import unittest
import urllib.request


# Enable warnings (see https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestSegmentation(unittest.TestCase):
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
        os.remove("HCT-116_PBase.ccf")
        self.assertEqual(result.total_blocks, 30172)
        self.assertEqual(result.total_fitness, 10271168.936047826)

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
        os.remove("HCT-116_SP1-PBase.ccf")
        self.assertEqual(result.total_blocks, 21375)
        self.assertEqual(result.total_fitness, 2987508.577343301)

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
        os.remove("HCT-116_HyPBase.ccf")
        self.assertEqual(result.total_blocks, 106211)
        self.assertEqual(result.total_fitness, 31079910.017576993)

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
        os.remove("HCT-116_SP1-HyPBase.ccf")
        self.assertEqual(result.total_blocks, 68935)
        self.assertEqual(result.total_fitness, 11410811.054755056)

    def test_uniformity(self):
        # Test that the algorithm segments a truly uniform dataset into a single block
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "data/test_uniform.qbed",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        self.assertEqual(result.total_blocks, 1)

    def test_equivalence(self):
        # Test that segmenting using OP or PELT yield the same segmentation
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args_OP = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--method",
                "OP",
            ]
        )
        result_OP = segmentation.segment_from_command_line(args_OP)
        args_PELT = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--method",
                "PELT",
            ]
        )
        result_PELT = segmentation.segment_from_command_line(args_PELT)
        os.remove("HCT-116_PBase.ccf")
        self.assertEqual(result_OP.total_blocks, result_PELT.total_blocks)
        # The fitness calculation is identical between the algorithms, except for a sign change
        self.assertEqual(result_OP.total_fitness, -result_PELT.total_fitness)


class TestSegmentationParameters(unittest.TestCase):
    def test_p0(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--p0",
                "2",
            ]
        )
        with self.assertRaises(AssertionError):
            segmentation.segment_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")

    def test_prior(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "HCT-116_PBase.ccf",
                "--prior",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            segmentation.segment_from_command_line(args)
        os.remove("HCT-116_PBase.ccf")


class TestSegmentationAPI(unittest.TestCase):
    def test_segmentation_API(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urllib.request.urlopen(PBase_url) as response, open("HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        result = segmentation.segment(
            BedTool("HCT-116_PBase.ccf"),
            "PELT",
            p0=DEFAULT_SEGMENTATION_P0,
            prior=None
        )
        os.remove("HCT-116_PBase.ccf")
        self.assertEqual(result.total_blocks, 30172)
        self.assertEqual(result.total_fitness, 10271168.936047826)


if __name__ == "__main__":
    unittest.main()
