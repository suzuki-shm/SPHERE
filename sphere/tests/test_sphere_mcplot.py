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
        self.__input_e_vm = d_dir + "/input_estimated_vm.tsv"
        self.__input_e_vm_K2 = d_dir + "/input_estimated_vm_K2.tsv"
        self.__input_e_vm_opt = d_dir + "/input_estimated_vm_opt.tsv"
        self.__input_e_ssc = d_dir + "/input_estimated_ssc.tsv"
        self.__input_e_sswc = d_dir + "/input_estimated_sswc.tsv"
        self.__input_e_ssvm = d_dir + "/input_estimated_ssvm.tsv"
        self.__input_e_sst = d_dir + "/input_estimated_sst.tsv"
        self.__input_e_ssl = d_dir + "/input_estimated_ssl.tsv"
        self.__input_e_t = d_dir + "/input_estimated_t.tsv"
        self.__input_e_l = d_dir + "/input_estimated_l.tsv"
        self.__output = d_dir + "/output.png"

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_mcplot_argument_parse(self):
        argv_str = "{0} {1} {2} 0".format(
            self.__output,
            self.__input_d,
            self.__input_e_wc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e_wc,
            "output_dest": self.__output,
            "index": 0,
            "pn": 50,
            "model_type": "vonmises",
            "fs": 18,
            "M": "sampling"
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_mcplot_command_lc(self):
        argv_str = "{0} {1} {2} 0 -m linearcardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_lc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_c(self):
        argv_str = "{0} {1} {2} 0 -m cardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_c
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_wc(self):
        argv_str = "{0} {1} {2} 0 -m wrappedcauchy".format(
            self.__output,
            self.__input_d,
            self.__input_e_wc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_vm(self):
        argv_str = "{0} {1} {2} 0 -m vonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_vm
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_vm_K2(self):
        argv_str = "{0} {1} {2} 0 -m vonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_vm_K2
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_vm_opt(self):
        argv_str = "{0} {1} {2} 0 -m vonmises -M optimizing".format(
            self.__output,
            self.__input_d,
            self.__input_e_vm_opt
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_ssc(self):
        argv_str = "{0} {1} {2} 0 -m sscardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_ssc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_sswc(self):
        argv_str = "{0} {1} {2} 0 -m sswrappedcauchy".format(
            self.__output,
            self.__input_d,
            self.__input_e_sswc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_ssvm(self):
        argv_str = "{0} {1} {2} 0 -m ssvonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_ssvm
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_sst(self):
        argv_str = "{0} {1} {2} 0 -m statespacetrigonal".format(
            self.__output,
            self.__input_d,
            self.__input_e_sst
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_ssl(self):
        argv_str = "{0} {1} {2} 0 -m statespacelinear".format(
            self.__output,
            self.__input_d,
            self.__input_e_ssl
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_t(self):
        argv_str = "{0} {1} {2} 0 -m trigonal".format(
            self.__output,
            self.__input_d,
            self.__input_e_t
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_l(self):
        argv_str = "{0} {1} {2} 0 -m linear".format(
            self.__output,
            self.__input_d,
            self.__input_e_l
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)


if __name__ == '__main__':
    unittest.main()
