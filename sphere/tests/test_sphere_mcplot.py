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
        self.__input_e_aec = d_dir + "/input_estimated_aec.tsv"
        self.__input_e_aewc = d_dir + "/input_estimated_aewc.tsv"
        self.__input_e_aevm = d_dir + "/input_estimated_aevm.tsv"
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

    def test_sphere_mcplot_command_aec(self):
        argv_str = "{0} {1} {2} 0 -m aecardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_aec
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_aewc(self):
        argv_str = "{0} {1} {2} 0 -m aewrappedcauchy".format(
            self.__output,
            self.__input_d,
            self.__input_e_aewc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_aevm(self):
        argv_str = "{0} {1} {2} 0 -m aevonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_aevm
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dlc(self):
        argv_str = "{0} {1} {2} 0 -m dlinearcardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_lc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dc(self):
        argv_str = "{0} {1} {2} 0 -m dcardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_c
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dwc(self):
        argv_str = "{0} {1} {2} 0 -m dwrappedcauchy".format(
            self.__output,
            self.__input_d,
            self.__input_e_wc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dvm(self):
        argv_str = "{0} {1} {2} 0 -m dvonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_vm
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dvm_K2(self):
        argv_str = "{0} {1} {2} 0 -m dvonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_vm_K2
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dvm_opt(self):
        argv_str = "{0} {1} {2} 0 -m dvonmises -M optimizing".format(
            self.__output,
            self.__input_d,
            self.__input_e_vm_opt
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_aedc(self):
        argv_str = "{0} {1} {2} 0 -m aedcardioid".format(
            self.__output,
            self.__input_d,
            self.__input_e_aec
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_aedwc(self):
        argv_str = "{0} {1} {2} 0 -m aedwrappedcauchy".format(
            self.__output,
            self.__input_d,
            self.__input_e_aewc
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_aedvm(self):
        argv_str = "{0} {1} {2} 0 -m aedvonmises".format(
            self.__output,
            self.__input_d,
            self.__input_e_aevm
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)


if __name__ == '__main__':
    unittest.main()
