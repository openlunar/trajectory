import unittest
from test import test_patched_conic, test_propagate

def trajectory_test_suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromModule(test_patched_conic))
    suite.addTests(loader.loadTestsFromModule(test_propagate))

    return suite


if __name__ == '__main__':
    unittest.main()
