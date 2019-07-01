#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
import pandas as pd
from sphere import sphere_filter
from sphere.sphere_utils import get_logger


class SphereFilterTest(unittest.TestCase):
    logger = get_logger(__name__)

    def setUp(self):
        self.maxDiff = None
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_filter"
        self.__input1 = d_dir + "/input1.tsv"
        self.__input2 = d_dir + "/input2.tsv"
        self.__input3 = d_dir + "/input3.tsv"
        self.__input4 = d_dir + "/input4.tsv"
        self.__answer6 = d_dir + "/answer6.tsv"
        self.__answer7 = d_dir + "/answer7.tsv"
        self.__answer8 = d_dir + "/answer8.tsv"
        self.__answer9 = d_dir + "/answer9.tsv"
        self.__answer10 = d_dir + "/answer10.tsv"
        self.__answer11 = d_dir + "/answer11.tsv"
        self.__output = d_dir + "/output.tsv"

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_filter_main(self):
        args = {
            "depth_file_path": self.__input1,
            "output_dest": self.__output,
            "s": 10,
            "w": 10,
            "t": "median"
        }
        sphere_filter.main(args, SphereFilterTest.logger)

    def test_sphere_filter_main_non_devidable(self):
        args = {
            "depth_file_path": self.__input1,
            "output_dest": self.__output,
            "s": 7,
            "w": 11,
            "t": "median"
        }
        sphere_filter.main(args, SphereFilterTest.logger)

    def test_sphere_filter_command1(self):
        argv_str = "{0} {1} -s 10 -w 10".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

    def test_sphere_filter_command2(self):
        argv_str = "{0} {1} -s 1 -w 6".format(self.__output, self.__input2)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

    def test_sphere_filter_command3(self):
        argv_str = "{0} {1} -s 2 -w 3".format(self.__output, self.__input2)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

    def test_sphere_filter_command4(self):
        argv_str = "{0} {1} -s 1 -w 1".format(self.__output, self.__input2)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

    def test_sphere_filter_command5(self):
        argv_str = "{0} {1} -t variance".format(self.__output, self.__input2)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__input2, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        self.assertEqual(df1.shape, df2.shape)

    def test_sphere_filter_command6(self):
        argv_str = "{0} {1} -t percentile -m 99999".format(
            self.__output,
            self.__input2
        )
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__answer6, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        pd.testing.assert_frame_equal(df1, df2)

    def test_sphere_filter_command7(self):
        argv_str = "{0} {1} -t percentile -m 0".format(
            self.__output,
            self.__input2
        )
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__answer7, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        pd.testing.assert_frame_equal(df1, df2)

    def test_sphere_filter_command8(self):
        argv_str = "{0} {1} -t fill -m 99999".format(
            self.__output,
            self.__input3
        )
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__answer8, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        pd.testing.assert_frame_equal(df1, df2)

    def test_sphere_filter_command9(self):
        argv_str = "{0} {1} -s 1 -w 3".format(self.__output, self.__input3)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__answer9, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        pd.testing.assert_frame_equal(df1, df2)

    def test_sphere_filter_command10(self):
        argv_str = "{0} {1} -s 1 -w 2".format(self.__output, self.__input4)
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__answer10, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        pd.testing.assert_frame_equal(df1, df2)

    def test_sphere_filter_command11(self):
        argv_str = "{0} {1} -s 10 -w 10 -t sum".format(
            self.__output,
            self.__input1
        )
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)

        df1 = pd.read_csv(self.__answer11, sep="\t", header=None)
        df2 = pd.read_csv(self.__output, sep="\t", header=None)
        pd.testing.assert_frame_equal(df1, df2)

    def test_sphere_filter_command12(self):
        argv_str = "{0} {1} -s 1 -w 3 -t mvariance".format(
            self.__output,
            self.__input2
        )
        argv = argv_str.split()
        args = sphere_filter.argument_parse(argv)
        sphere_filter.main(args, SphereFilterTest.logger)


if __name__ == '__main__':
    unittest.main()
