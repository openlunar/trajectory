import unittest
from test import test_patched_conic

def trajectory_test_suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromModule(test_patched_conic))

    return suite


if __name__ == '__main__':
    unittest.main()
