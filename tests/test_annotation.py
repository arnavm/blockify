from blockify.parsers import blockify_parser
import blockify.annotation as annotation
import gzip
import lzma
import os
from pybedtools import BedTool
import sys
import unittest
from urllib.request import Request, urlopen


# Enable warnings
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


# Global variable
hg38_TTAA_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/TTAA/hg38_TTAA.bed.xz"


class TestAnnotation(unittest.TestCase):
    def test_yeast(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--alpha",
                "0.05",
                "--correction",
                "bonferroni",
                "--distance",
                "0",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        self.assertEqual(len(result), 236)

    def test_PBase(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urlopen(PBase_url) as response, open("tests/data/HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        PBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_PBase.blocks"
        with urlopen(Request(PBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_PBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        with urlopen(Request(hg38_TTAA_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/hg38_TTAA.bed", 'wb') as out_file:
            out_file.write(lzma.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT-116_PBase.ccf",
                "--regions",
                "tests/data/HCT-116_PBase.blocks",
                "--background",
                "tests/data/hg38_TTAA.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-30",
                "--pseudocount",
                "0",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        os.remove("tests/data/HCT-116_PBase.ccf")
        os.remove("tests/data/HCT-116_PBase.blocks")
        os.remove("tests/data/hg38_TTAA.bed")
        self.assertEqual(len(result), 1939)

    def test_SP1_PBase(self):
        PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
        with urlopen(PBase_url) as response, open("tests/data/HCT-116_PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        SP1_PBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471637&format=file&file=GSM4471637%5FHCT%2D116%5FSP1%2DPBase%2Eccf%2Etxt%2Egz"
        with urlopen(SP1_PBase_url) as response, open("tests/data/HCT-116_SP1-PBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        SP1_PBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_SP1-PBase.blocks"
        with urlopen(Request(SP1_PBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_SP1-PBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT-116_SP1-PBase.ccf",
                "--regions",
                "tests/data/HCT-116_SP1-PBase.blocks",
                "--background",
                "tests/data/HCT-116_PBase.ccf",
                "--distance",
                "250",
                "--pValueCutoff",
                "1e-6",
                "--pseudocount",
                "0.1",
                "--tight",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        os.remove("tests/data/HCT-116_SP1-PBase.ccf")
        os.remove("tests/data/HCT-116_SP1-PBase.blocks")
        os.remove("tests/data/HCT-116_PBase.ccf")
        self.assertEqual(len(result), 5067)

    def test_HyPBase(self):
        HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471638&format=file&file=GSM4471638%5FHCT%2D116%5FHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(HyPBase_url) as response, open("tests/data/HCT-116_HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        HyPBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_HyPBase.blocks"
        with urlopen(Request(HyPBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_HyPBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        with urlopen(Request(hg38_TTAA_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/hg38_TTAA.bed", 'wb') as out_file:
            out_file.write(lzma.decompress(response.read()))
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT-116_HyPBase.ccf",
                "--regions",
                "tests/data/HCT-116_HyPBase.blocks",
                "--background",
                "tests/data/hg38_TTAA.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-62",
                "--pseudocount",
                "0",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        os.remove("tests/data/HCT-116_HyPBase.ccf")
        os.remove("tests/data/HCT-116_HyPBase.blocks")
        os.remove("tests/data/hg38_TTAA.bed")
        self.assertEqual(len(result), 1956)

    def test_SP1_HyPBase(self):
        HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471638&format=file&file=GSM4471638%5FHCT%2D116%5FHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(HyPBase_url) as response, open("tests/data/HCT-116_HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        SP1_HyPBase_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471639&format=file&file=GSM4471639%5FHCT%2D116%5FSP1%2DHyPBase%2Eccf%2Etxt%2Egz"
        with urlopen(SP1_HyPBase_url) as response, open("tests/data/HCT-116_SP1-HyPBase.ccf", 'wb') as out_file:
            out_file.write(gzip.decompress(response.read()))
        SP1_HyPBase_blocks_URL = "https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/Blocks/HCT-116_SP1-HyPBase.blocks"
        with urlopen(Request(SP1_HyPBase_blocks_URL, headers={'User-Agent': 'Mozilla/5.0'})) as response, open("tests/data/HCT-116_SP1-HyPBase.blocks", 'wb') as out_file:
            out_file.write(response.read())
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT-116_SP1-HyPBase.ccf",
                "--regions",
                "tests/data/HCT-116_SP1-HyPBase.blocks",
                "--background",
                "tests/data/HCT-116_HyPBase.ccf",
                "--distance",
                "250",
                "--pValueCutoff",
                "1e-22",
                "--pseudocount",
                "0.1",
                "--tight",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        os.remove("tests/data/HCT-116_SP1-HyPBase.ccf")
        os.remove("tests/data/HCT-116_SP1-HyPBase.blocks")
        os.remove("tests/data/HCT-116_HyPBase.ccf")
        self.assertEqual(len(result), 5049)


class TestAnnotationParameters(unittest.TestCase):
    def test_pseudocount(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--alpha",
                "0.05",
                "--correction",
                "bonferroni",
                "--distance",
                "0",
                "--pseudocount",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_alpha(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--alpha",
                "0.05",
                "--correction",
                "bonferroni",
                "--distance",
                "0",
                "--alpha",
                "1.1",
                "--correction",
                "fdr_bh",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_pValueCutoff(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--distance",
                "0",
                "--pValueCutoff",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_distance(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--alpha",
                "0.05",
                "--correction",
                "bonferroni",
                "--distance",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_min(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--alpha",
                "0.05",
                "--correction",
                "bonferroni",
                "--distance",
                "0",
                "--min",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_max(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/S288C_CBF1.qbed",
                "--regions",
                "tests/data/S288C_CBF1.blocks",
                "--background",
                "tests/data/S288C_dSIR4.qbed",
                "--alpha",
                "0.05",
                "--correction",
                "bonferroni",
                "--distance",
                "0",
                "--max",
                "-1",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)


class TestAnnotationAPI(unittest.TestCase):
    def test_annotation_API(self):
        result, intermediate = annotation.annotate(
            BedTool("tests/data/S288C_CBF1.qbed"),
            BedTool("tests/data/S288C_CBF1.blocks"),
            BedTool("tests/data/S288C_dSIR4.qbed"),
            measure="enrichment",
            intermediate=True,
            distance=0,
            alpha=0.05,
            correction="bonferroni",
        )
        self.assertEqual(len(result), 236)
        self.assertEqual(len(intermediate), 1408)


if __name__ == "__main__":
    unittest.main()
