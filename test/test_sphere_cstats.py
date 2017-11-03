#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28


import unittest
import os
from sphere import sphere_cstats
from sphere.sphere_utils import get_logger


class SpherecstatsTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_cstats"
        self.__output = d_dir + "/output.tsv"
        self.__input1 = d_dir + "/input1.tsv"
        self.__input2 = d_dir + "/input2.tsv"
        self.__input = [self.__input1, self.__input2]
        self.__logger = get_logger(__name__)

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_cstats_main_multi(self):
        args = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
        }
        sphere_cstats.main(args, self.__logger)

    def test_sphere_cstats_main_single(self):
        args = {
            "depth_file_path": [self.__input1],
            "output_dest": self.__output,
        }
        sphere_cstats.main(args, self.__logger)

    def test_sphere_cstats_argument_parse_multi(self):
        argv_str = "{0} {1} {2}".format(self.__output,
                                        self.__input1,
                                        self.__input2)
        argv = argv_str.split()
        args = sphere_cstats.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_cstats_argument_parse_single(self):
        argv_str = "{0} {1}".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_cstats.argument_parse(argv)
        args_answer = {
            "depth_file_path": [self.__input1],
            "output_dest": self.__output
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_cstats_command_multi(self):
        argv_str = "{0} {1} {2}".format(self.__output,
                                        self.__input1,
                                        self.__input2)
        argv = argv_str.split()
        args = sphere_cstats.argument_parse(argv)
        sphere_cstats.main(args, self.__logger)

    def test_sphere_cstats_command_single(self):
        argv_str = "{0} {1}".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_cstats.argument_parse(argv)
        sphere_cstats.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
