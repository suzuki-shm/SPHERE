#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26


import unittest
import os
from sphere import sphere_dplot
from sphere.sphere_utils import get_logger


class SphereDplotTest(unittest.TestCase):
    logger = get_logger(__name__)

    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_dplot"
        self.__input1 = d_dir + "/input1.tsv"
        self.__input2 = d_dir + "/input2.tsv"
        self.__output = d_dir + "/output.png"

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_dplot_main_np5(self):
        args = {
            "depth_file_path": self.__input1,
            "output_dest": self.__output,
            "np": 5,
            "fs": 30,
        }
        sphere_dplot.main(args, SphereDplotTest.logger)

    def test_sphere_dplot_main_np29(self):
        args = {
            "depth_file_path": self.__input1,
            "output_dest": self.__output,
            "np": 29,
            "fs": 30,
        }
        sphere_dplot.main(args, SphereDplotTest.logger)

    def test_sphere_dplot_main_cl101(self):
        args = {
            "depth_file_path": self.__input1,
            "output_dest": self.__output,
            "np": 10,
            "fs": 30,
        }
        sphere_dplot.main(args, SphereDplotTest.logger)

    def test_sphere_dplot_argument_parse(self):
        argv_str = "{0} {1}".format(self.__input1, self.__output)
        argv = argv_str.split()
        args = sphere_dplot.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input1,
            "output_dest": self.__output,
            "np": 30,
            "fs": 30,
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_dplot_command1(self):
        argv_str = "{0} {1}".format(self.__input1, self.__output)
        argv = argv_str.split()
        args = sphere_dplot.argument_parse(argv)
        sphere_dplot.main(args, SphereDplotTest.logger)

    def test_sphere_dplot_command2(self):
        argv_str = "{0} {1}".format(self.__input2, self.__output)
        argv = argv_str.split()
        args = sphere_dplot.argument_parse(argv)
        sphere_dplot.main(args, SphereDplotTest.logger)


if __name__ == '__main__':
    unittest.main()
