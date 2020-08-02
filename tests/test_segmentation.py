from blockify.parsers import blockify_parser
from blockify.parsers import DEFAULT_SEGMENTATION_P0
import blockify.segmentation as segmentation
import gzip
import os
from pybedtools import BedTool
import shutil
import re
import sys
import unittest
import urllib.request


# Enable warnings (see https://docs.python.org/3/library/warnings.html#overriding-the-default-filter)
if not sys.warnoptions:
    import os
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
                # "test.bed"
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
                # "test.bed",
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
                # "test.bed",
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
                # "test.bed",
            ]
        )
        result = segmentation.segment_from_command_line(args)
        os.remove("HCT-116_SP1-HyPBase.ccf")
        self.assertEqual(result.total_blocks, 68935)
        self.assertEqual(result.total_fitness, 11410811.054755056)

    def test_uniformity(self):
        # Test that the algorithm segments a truly uniform dataset into a single block
        with open("uniform_test.ccf", 'w') as out_file:
            uniform_data = """
                chr1	1	2	1
                chr1	2	3	1
                chr1	3	4	1
                chr1	4	5	1
                chr1	5	6	1
                chr1	6	7	1
                chr1	7	8	1
                chr1	8	9	1
                chr1	9	10	1
                chr1	10	11	1
                chr1	11	12	1
                chr1	12	13	1
                chr1	13	14	1
                chr1	14	15	1
                chr1	15	16	1
                chr1	16	17	1
                chr1	17	18	1
                chr1	18	19	1
                chr1	19	20	1
                chr1	20	21	1
                chr1	21	22	1
                chr1	22	23	1
                chr1	23	24	1
                chr1	24	25	1
                chr1	25	26	1
                chr1	26	27	1
                chr1	27	28	1
                chr1	28	29	1
                chr1	29	30	1
                chr1	30	31	1
                chr1	31	32	1
                chr1	32	33	1
                chr1	33	34	1
                chr1	34	35	1
                chr1	35	36	1
                chr1	36	37	1
                chr1	37	38	1
                chr1	38	39	1
                chr1	39	40	1
                chr1	40	41	1
                chr1	41	42	1
                chr1	42	43	1
                chr1	43	44	1
                chr1	44	45	1
                chr1	45	46	1
                chr1	46	47	1
                chr1	47	48	1
                chr1	48	49	1
                chr1	49	50	1
                chr1	50	51	1
                chr1	51	52	1
                chr1	52	53	1
                chr1	53	54	1
                chr1	54	55	1
                chr1	55	56	1
                chr1	56	57	1
                chr1	57	58	1
                chr1	58	59	1
                chr1	59	60	1
                chr1	60	61	1
                chr1	61	62	1
                chr1	62	63	1
                chr1	63	64	1
                chr1	64	65	1
                chr1	65	66	1
                chr1	66	67	1
                chr1	67	68	1
                chr1	68	69	1
                chr1	69	70	1
                chr1	70	71	1
                chr1	71	72	1
                chr1	72	73	1
                chr1	73	74	1
                chr1	74	75	1
                chr1	75	76	1
                chr1	76	77	1
                chr1	77	78	1
                chr1	78	79	1
                chr1	79	80	1
                chr1	80	81	1
                chr1	81	82	1
                chr1	82	83	1
                chr1	83	84	1
                chr1	84	85	1
                chr1	85	86	1
                chr1	86	87	1
                chr1	87	88	1
                chr1	88	89	1
                chr1	89	90	1
                chr1	90	91	1
                chr1	91	92	1
                chr1	92	93	1
                chr1	93	94	1
                chr1	94	95	1
                chr1	95	96	1
                chr1	96	97	1
                chr1	97	98	1
                chr1	98	99	1
                """
            out_file.write(re.sub(r'^\s+', '', uniform_data, flags=re.M))
        args = blockify_parser.parse_args(
            [
                "segment",
                "--input",
                "uniform_test.ccf",
                # "test.bed"
            ]
        )
        result = segmentation.segment_from_command_line(args)
        os.remove("uniform_test.ccf")
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
                # "test.bed",
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
                # "test.bed",
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
                # "test.bed",
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
                # "test.bed",
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
