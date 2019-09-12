from blockify.parsers import blockify_parser
import blockify.annotation as annotation
from pybedtools import BedTool
import sys
import unittest


# Enable warnings
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestAnnotation(unittest.TestCase):
    def test_PBase(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-30",
                "--pseudocount",
                "0",
                "test.bed",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        self.assertEqual(len(result), 1939)

    def test_SP1_PBase(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-SP1-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-SP1-PBase_sort_canon.blocks",
                "--background",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--distance",
                "250",
                "--pValueCutoff",
                "1e-6",
                "--pseudocount",
                "0.1",
                "--tight",
                "test.bed",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        self.assertEqual(len(result), 5067)

    def test_HyPBase(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-HyPBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-HyPBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-62",
                "--pseudocount",
                "0",
                "test.bed",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        self.assertEqual(len(result), 1956)

    def test_SP1_HyPBase(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-SP1-HyPBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-SP1-HyPBase_sort_canon.blocks",
                "--background",
                "tests/data/HCT116-HyPBase_sort_canon.ccf",
                "--distance",
                "250",
                "--pValueCutoff",
                "1e-22",
                "--pseudocount",
                "0.1",
                "--tight",
                "test.bed",
            ]
        )
        result, intermediate = annotation.annotate_from_command_line(args)
        self.assertEqual(len(result), 5049)


class TestAnnotationParameters(unittest.TestCase):
    def test_pseudocount(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-3",
                "--pseudocount",
                "-1",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_alpha(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--alpha",
                "1.1",
                "--correction",
                "fdr_bh",
                "--pseudocount",
                "0",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_pValueCutoff(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "-1",
                "--pseudocount",
                "0",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_distance(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "-1",
                "--pValueCutoff",
                "1e-3",
                "--pseudocount",
                "0",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_min(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-3",
                "--pseudocount",
                "0",
                "--min",
                "-1",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)

    def test_max(self):
        args = blockify_parser.parse_args(
            [
                "call",
                "--input",
                "tests/data/HCT116-PBase_sort_canon.ccf",
                "--regions",
                "tests/data/HCT116-PBase_sort_canon.blocks",
                "--background",
                "tests/data/hg38_TTAA_canon.bed",
                "--distance",
                "12500",
                "--pValueCutoff",
                "1e-3",
                "--pseudocount",
                "0",
                "--max",
                "-1",
                "test.bed",
            ]
        )
        with self.assertRaises(AssertionError):
            annotation.annotate_from_command_line(args)


class TestAnnotationAPI(unittest.TestCase):
    def test_annotation_API(self):
        result, intermediate = annotation.annotate(
            BedTool("tests/data/HCT116-PBase_sort_canon.ccf"),
            BedTool("tests/data/HCT116-PBase_sort_canon.blocks"),
            BedTool("tests/data/hg38_TTAA_canon.bed"),
            intermediate=True,
            distance=12500,
            p_value=1e-30,
            pseudocount=0
        )
        self.assertEqual(len(result), 1939)
        self.assertEqual(len(intermediate), 30172)


if __name__ == "__main__":
    unittest.main()
