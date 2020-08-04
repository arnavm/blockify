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


# Create a basic test suite
basic_tests = unittest.TestSuite()
loader = unittest.TestLoader()

# Add segmentation tests
basic_tests.addTest(test_segmentation.TestSegmentation("test_yeast"))
basic_tests.addTest(test_segmentation.TestSegmentation("test_uniformity"))
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_segmentation.TestSegmentationParameters)
)
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_segmentation.TestSegmentationAPI)
)

# Add normalization tests
basic_tests.addTest(test_normalization.TestNormalization("test_yeast"))
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_normalization.TestNormalizationParameters)
)
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_normalization.TestNormalizationAPI)
)

# Add annotation (peak calling) tests
basic_tests.addTest(test_annotation.TestAnnotation("test_yeast"))
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_annotation.TestAnnotationParameters)
)
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_annotation.TestAnnotationAPI)
)

# Add utilities tests
basic_tests.addTests(
    loader.loadTestsFromTestCase(test_utilities.TestUtilities)
)

# Run tests
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(basic_tests)
