#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
import pandas as pd
from sphere import sphere_estimate
from sphere.sphere_utils import get_logger


class SphereEstimateTest(unittest.TestCase):
    logger = get_logger(__name__)

    def setUp(self):
        self.maxDiff = None
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_estimate"
        self.__input = [d_dir+"/input_1.tsv", d_dir+"/input_2.tsv"]
        self.__input_mix = [d_dir + "/input_3.tsv"]
        self.__output = d_dir + "/output.tsv"
        self.__output_fit = d_dir + "/fit.pkl"
        self.__output_ll = d_dir + "/log_lik.tsv"

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)
        if os.path.exists(self.__output_fit):
            os.remove(self.__output_fit)
        if os.path.exists(self.__output_ll):
            os.remove(self.__output_ll)

    # argument parse evaluation
    def test_sphere_estimate_argument_parse_vonmises(self):
        argv_str = """{0} {1} -fod {2}
                       -lld {3} -sc 1 -si 50 -sw 20 -ff""".format(
            self.__output,
            self.__input[0],
            self.__output_fit,
            self.__output_ll
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "depth_file_path": [self.__input[0]],
            "fod": self.__output_fit,
            "lld": self.__output_ll,
            "m": "vonmises",
            "M": "sampling",
            "nmix": 1,
            "si": 50,
            "sw": 20,
            "sc": 1,
            "st": 1,
            "ss": 1234,
            "ff": True,
            "p": None,
            "ll": False,
            "j": -1,
            "om": "LBFGS"
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_estimate_argument_parse_cardioid(self):
        argv_str = """{0} {1} -m cardioid""".format(
            self.__output,
            self.__input[0]
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "depth_file_path": [self.__input[0]],
            "fod": None,
            "lld": None,
            "m": "cardioid",
            "M": "sampling",
            "nmix": 1,
            "si": 3000,
            "sw": 1000,
            "sc": 1,
            "st": 1,
            "ss": 1234,
            "ff": False,
            "p": None,
            "ll": False,
            "j": -1,
            "om": "LBFGS"
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_estimate_argument_parse_wrappedcauchy(self):
        argv_str = """{0} {1} -m wrappedcauchy""".format(
            self.__output,
            self.__input[0]
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        args_answer = {
            "output_dest": self.__output,
            "depth_file_path": [self.__input[0]],
            "fod": None,
            "lld": None,
            "m": "wrappedcauchy",
            "M": "sampling",
            "nmix": 1,
            "si": 3000,
            "sw": 1000,
            "sc": 1,
            "st": 1,
            "ss": 1234,
            "ff": False,
            "p": None,
            "ll": False,
            "j": -1,
            "om": "LBFGS"
        }
        self.assertDictEqual(args, args_answer)

    # Full test for sampling
    def test_sphere_estimate_command_sampling_vm_single(self):
        argv_str = """{0} {1} -m vonmises -sc 1 -si 30 -sw 20 -ff""".format(
            self.__output,
            self.__input[0],
            self.__output_fit,
            self.__output_ll
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_sampling_lc_single(self):
        argv_str = """{0} {1} -m linearcardioid -sc 1 -si 30 -sw 20
                      -ff""".format(
            self.__output,
            self.__input[0],
            self.__output_fit,
            self.__output_ll
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    # Check if log_lik destination is used but -ll frag is not used.
    def test_sphere_estimate_command_sampling_vm_single_ll_lld(self):
        argv_str = """{0} {1} -lld {2} -sc 1 -si 30 -sw 20 -ff""".format(
            self.__output,
            self.__input[0],
            self.__output_ll
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    # Full test for optimizing; single
    def test_sphere_estimate_command_optimizing_lc_single(self):
        argv_str = """{0} {1} -M optimizing -m linearcardioid
                      -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_c_single(self):
        argv_str = """{0} {1} -M optimizing -m cardioid -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_wc_single(self):
        argv_str = """{0} {1} -M optimizing -m wrappedcauchy -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_vm_single(self):
        argv_str = """{0} {1} -M optimizing -m vonmises -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_jp_single(self):
        argv_str = """{0} {1} -M optimizing -m jonespewsey -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaec_single(self):
        argv_str = """{0} {1} -M optimizing -m miaecardioid -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaewc_single(self):
        argv_str = """{0} {1} -M optimizing -m miaewrappedcauchy -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaevm_single(self):
        argv_str = """{0} {1} -M optimizing -m miaevonmises -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaejp_single(self):
        argv_str = """{0} {1} -M optimizing -m miaejonespewsey -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_dvm_single(self):
        argv_str = """{0} {1} -M optimizing -m dvonmises -ff""".format(
            self.__output,
            self.__input[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    # Full test for optimizing; multiple
    def test_sphere_estimate_command_optimizing_vm_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m vonmises -ff""".format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_lc_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m linearcardioid
                      -ff""".format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaevm_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m miaevonmises -ff""".format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaec_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m miaecardioid
                      -ff""".format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaewc_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m miaewrappedcauchy
            -ff""".format(
                self.__output,
                self.__input[0],
                self.__input[1],
            )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_miaejp_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m miaejonespewsey
            -ff""".format(
                self.__output,
                self.__input[0],
                self.__input[1],
            )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_invmiaevm_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m invmiaevonmises
            -ff""".format(
                self.__output,
                self.__input[0],
                self.__input[1],
            )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_invmiaec_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m invmiaecardioid
                      -ff""".format(
            self.__output,
            self.__input[0],
            self.__input[1],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_invmiaewc_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m invmiaewrappedcauchy
            -ff""".format(
                self.__output,
                self.__input[0],
                self.__input[1],
            )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def test_sphere_estimate_command_optimizing_invmiaejp_multiple(self):
        argv_str = """{0} {1} {2} -M optimizing -m invmiaejonespewsey
            -ff""".format(
                self.__output,
                self.__input[0],
                self.__input[1],
            )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    # Full test for optimizing; single mix
    def test_sphere_estimate_command_optimizing_vm_single_mix(self):
        argv_str = """{0} {1} -M optimizing -m vonmises -nmix 2 -ff""".format(
            self.__output,
            self.__input_mix[0],
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereEstimateTest.logger)

    def assertStanConvergence(self, args):
        df = pd.read_csv(args["output_dest"],
                         sep="\t",
                         index_col=0)
        n_not_converted = len(df[df["Rhat"] >= 1.1])
        self.assertEqual(n_not_converted, 0)


if __name__ == '__main__':
    unittest.main()
