#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
from sphere import sphere_rdump
from sphere.sphere_utils import get_logger


class SphereRdumpTest(unittest.TestCase):
    logger = get_logger(__name__)

    def setUp(self):
        self.maxDiff = None
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_rdump"
        self.__input = [d_dir+"/input_1.tsv", d_dir+"/input_2.tsv"]
        self.__output = d_dir + "/output.tsv"

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    # argument parse evaluation
    def test_sphere_rdump_argument_parse_single(self):
        argv_str = """{0} {1}""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_rdump.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "depth_file_path": [self.__input[0]],
            "nmix": 1,
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_rdump_argument_parse_multiple(self):
        argv_str = """{0} {1} {2}""".format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_rdump.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "depth_file_path": self.__input,
            "nmix": 1,
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_rdump_command_single(self):
        argv_str = """{0} {1}""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_rdump.argument_parse(argv)
        sphere_rdump.main(args, SphereRdumpTest.logger)

    def test_sphere_rdump_command_multiple(self):
        argv_str = """{0} {1} {2} """.format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_rdump.argument_parse(argv)
        sphere_rdump.main(args, SphereRdumpTest.logger)


if __name__ == '__main__':
    unittest.main()
