#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
from sphere import sphere_estimate


class SphereEstimateTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_estimate"
        self.__input = d_dir + "/input.tsv"
        self.__output = d_dir + "/output.tsv"
        self.__output_model = d_dir + "/model.pkl"
        self.__output_fit = d_dir + "/fit.pkl"
        self.__output_ll = d_dir + "/log_lik.tsv"
        self.__logger = sphere_estimate.get_logger()

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)
        if os.path.exists(self.__output_model):
            os.remove(self.__output_model)
        if os.path.exists(self.__output_fit):
            os.remove(self.__output_fit)
        if os.path.exists(self.__output_ll):
            os.remove(self.__output_ll)

    def test_sphere_estimate_main(self):
        args = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
            "pmd": self.__output_model,
            "pmp": None,
            "m": "trigonal",
            "fod": self.__output_fit,
            "lld": self.__output_ll,
            "cl": 100,
            "si": 50,
            "sw": 20,
            "sc": 1,
            "st": 1,
            "ss": 1234,
            "ff": True,
            "p": None
        }
        sphere_estimate.main(args, self.__logger)

    def test_sphere_estimate_main_linear(self):
        args = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
            "pmd": self.__output_model,
            "pmp": None,
            "m": "linear",
            "fod": self.__output_fit,
            "lld": self.__output_ll,
            "cl": 100,
            "si": 50,
            "sw": 20,
            "sc": 1,
            "st": 1,
            "ss": 1234,
            "ff": True,
            "p": None
        }
        sphere_estimate.main(args, self.__logger)

    def test_sphere_estimate_main_mix(self):
        args = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
            "pmd": self.__output_model,
            "pmp": None,
            "m": "mix",
            "fod": self.__output_fit,
            "lld": self.__output_ll,
            "cl": 100,
            "si": 50,
            "sw": 20,
            "sc": 1,
            "st": 1,
            "ss": 1234,
            "ff": True,
            "p": None
        }
        sphere_estimate.main(args, self.__logger)

    def test_sphere_estimate_argument_parse(self):
        argv_str = "{0} {1} -pmd {2} -fod {3} -lld {4}".format(
            self.__input,
            self.__output,
            self.__output_model,
            self.__output_fit,
            self.__output_ll
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
            "pmd": self.__output_model,
            "pmp": None,
            "fod": self.__output_fit,
            "lld": self.__output_ll,
            "m": "trigonal",
            "cl": 10000,
            "si": 3000,
            "sw": 1000,
            "sc": 3,
            "st": 1,
            "ss": 1234,
            "ff": False,
            "p": None
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_estimate_command(self):
        argv_str = """{0} {1} -pmd {2} -fod {3} -lld {4}
                       -sc 1 -si 30 -sw 20""".format(
            self.__input,
            self.__output,
            self.__output_model,
            self.__output_fit,
            self.__output_ll
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
