import unittest
from testytools import result
from os.path import isfile, join, splitext
from os import listdir

class TestBoundsClass(unittest.TestCase):

    def test_valid_bounds(self):
        result.Bounds(2,8)
        result.Bounds(1412, 124544, 4000)
        result.Bounds(123.456, 782)
        result.Bounds(3.14, 6.23)
        result.Bounds('a','z', 'c')
        result.Bounds(3.14, 6.23, None)

    @unittest.expectedFailure
    def test_invalid_bounds(self):
        result.Bounds(123512, 435)

    @unittest.expectedFailure
    def test_invalid_point(self):
        result.Bounds(6.253112, 9.231521, 12.3)

    def test_bounds_equality(self):
        self.assertTrue(result.Bounds(5.12, 10.452) == \
                        result.Bounds(7.1298, 15.21))
        self.assertTrue(result.Bounds(100.13, 3.1e5, 600) == \
                        result.Bounds(55, 100.13, 56))
        self.assertFalse(result.Bounds(1,6) == \
                         result.Bounds(7,1e8))

class TestResultClass(unittest.TestCase):
    def test_results(self):
        dirname =\
        '/usr/project/dlab/Users/gth/projects/SHARKStar/debug_tests/200204_tests'
        files = [join(dirname,fp) for fp in listdir(dirname)
                 if isfile(join(dirname, fp)) and splitext(fp)[-1] == '.json']
        results = [result.Result.from_file(open(f)) for f in files]



if __name__ == '__main__':
    unittest.main()


