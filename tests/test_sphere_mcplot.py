#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
from sphere import sphere_mcplot
from sphere.sphere_utils import get_logger


class SphereMcplotTest(unittest.TestCase):
    logger = get_logger(__name__)

    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_mcplot"
        self.__input_d = d_dir + "/input_depth.tsv"
        self.__input_e_lc = d_dir + "/input_estimated_lc.tsv"
        self.__input_e_c = d_dir + "/input_estimated_c.tsv"
        self.__input_e_wc = d_dir + "/input_estimated_wc.tsv"
        self.__input_e_sswc = d_dir + "/input_estimated_sswc.tsv"
        self.__input_e_vm = d_dir + "/input_estimated_vm.tsv"
        self.__output = d_dir + "/output.png"

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_dplot_main_linearcardioid(self):
        args = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_lc,
            "output_dest": self.__output,
            "index": 0,
            "model_type": "linearcardioid",
            "fs": 18,
        }
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_dplot_main_cardioid(self):
        args = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_c,
            "output_dest": self.__output,
            "index": 0,
            "model_type": "cardioid",
            "fs": 18,
        }
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_dplot_main_wrappedcauchy(self):
        args = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_wc,
            "output_dest": self.__output,
            "index": 0,
            "model_type": "wrappedcauchy",
            "fs": 18,
        }
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_dplot_main_sswrappedcauchy(self):
        args = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_sswc,
            "output_dest": self.__output,
            "index": 0,
            "model_type": "sswrappedcauchy",
            "fs": 18,
        }
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_dplot_main_vonmises(self):
        args = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_vm,
            "output_dest": self.__output,
            "index": 0,
            "model_type": "vonmises",
            "fs": 18,
        }
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_argument_parse(self):
        argv_str = "{0} {1} {2} 0".format(self.__input_d,
                                          self.__input_e_wc,
                                          self.__output)
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_wc,
            "output_dest": self.__output,
            "index": 0,
            "model_type": "wrappedcauchy",
            "fs": 18
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_mcplot_command(self):
        argv_str = "{0} {1} {2} 0 -m wrappedcauchy".format(self.__input_d,
                                                           self.__input_e_wc,
                                                           self.__output)
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)


if __name__ == '__main__':
    unittest.main()
