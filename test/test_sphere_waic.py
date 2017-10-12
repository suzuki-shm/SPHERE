#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28


import unittest
import os
from sphere import sphere_waic


class SphereWaicTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_waic"
        self.__output = d_dir + "/output.tsv"
        self.__input1 = d_dir + "/input1.tsv"
        self.__input2 = d_dir + "/input2.tsv"
        self.__input = [self.__input1, self.__input2]
        self.__logger = sphere_waic.get_logger()

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_waic_main(self):
        args = {
            "log_lik_files": self.__input,
            "output_dest": self.__output,
            "t": "bda3"
        }
        sphere_waic.main(args, self.__logger)

    def test_sphere_waic_main_original(self):
        args = {
            "log_lik_files": self.__input,
            "output_dest": self.__output,
            "t": "original"
        }
        sphere_waic.main(args, self.__logger)

    def test_sphere_waic_main_both(self):
        args = {
            "log_lik_files": self.__input,
            "output_dest": self.__output,
            "t": "both"
        }
        sphere_waic.main(args, self.__logger)

    def test_sphere_waic_argument_parse_single(self):
        argv_str = "{0} {1}".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_waic.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "log_lik_files": [self.__input1],
            "t": "both"
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_waic_argument_parse_multi(self):
        argv_str = "{0} {1} {2}".format(self.__output,
                                        self.__input1,
                                        self.__input2)
        argv = argv_str.split()
        args = sphere_waic.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "log_lik_files": self.__input,
            "t": "both"
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_cstats_command_multi(self):
        argv_str = "{0} {1} {2}".format(self.__output,
                                        self.__input1,
                                        self.__input2)
        argv = argv_str.split()
        args = sphere_waic.argument_parse(argv)
        sphere_waic.main(args, self.__logger)

    def test_sphere_cstats_command_single(self):
        argv_str = "{0} {1}".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_waic.argument_parse(argv)
        sphere_waic.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
