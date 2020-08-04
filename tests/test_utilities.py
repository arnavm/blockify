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
        self.assertTrue(utilities.isSortedBEDObject(BedTool("tests/data/test_uniform_sorted.qbed")))

    def test_unsorted(self):
        self.assertFalse(utilities.isSortedBEDObject(BedTool("tests/data/test_uniform_unsorted.qbed")))


if __name__ == "__main__":
    unittest.main()
