import blockify.utilities as utilities
from pybedtools import BedTool
import sys
import unittest


# Enable warnings
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


class TestUtilities(unittest.TestCase):
    def test_sorted(self):
        self.assertTrue(utilities.isSortedBEDObject(BedTool("tests/data/HCT116-PBase_sort_canon.ccf")))

    def test_unsorted(self):
        self.assertFalse(utilities.isSortedBEDObject(BedTool("tests/data/HCT116-PBase_unsorted.ccf")))


if __name__ == "__main__":
    unittest.main()
