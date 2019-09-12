import tests.test_segmentation as test_segmentation
import tests.test_normalization as test_normalization
import tests.test_annotation as test_annotation
import tests.test_utilities as test_utilities
import sys
import unittest


# Enable warnings
if not sys.warnoptions:
    import os
    import warnings
    warnings.simplefilter("default") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses


# Create an advanced test suite
advanced_tests = unittest.TestSuite()
loader = unittest.TestLoader()

# Add segmentation tests
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_segmentation.TestSegmentation)
)
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_segmentation.TestSegmentationParameters)
)
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_segmentation.TestSegmentationAPI)
)

# Add normalization tests
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_normalization.TestNormalization)
)
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_normalization.TestNormalizationParameters)
)
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_normalization.TestNormalizationAPI)
)

# Add annotation (peak calling) tests
advanced_tests.addTests(loader.loadTestsFromTestCase(test_annotation.TestAnnotation))
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_annotation.TestAnnotationParameters)
)
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_annotation.TestAnnotationAPI)
)

# Add utilities tests
advanced_tests.addTests(
    loader.loadTestsFromTestCase(test_utilities.TestUtilities)
)

# Run tests
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(advanced_tests)
