#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
from sphere import sphere_mcplot
from sphere import sphere_estimate
from sphere.sphere_utils import get_logger


class SphereMcplotTest(unittest.TestCase):
    logger = get_logger(__name__)

    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_mcplot"
        self.__input_d = d_dir + "/input_depth.tsv"
        self.__input_e = d_dir + "/input_estimated.tsv"
        self.__output = d_dir + "/output.png"

    def tearDown(self):
        if os.path.exists(self.__input_e):
            os.remove(self.__input_e)
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_mcplot_argument_parse(self):
        argv_str = "{0} {1} {2} 0".format(
            self.__output,
            self.__input_d,
            self.__input_e
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e,
            "output_dest": self.__output,
            "index": 0,
            "pn": 50,
            "model_type": "vonmises",
            "fs": 18,
            "M": "sampling"
        }
        self.assertDictEqual(args, args_answer)

    # Test for sampling result
    def test_sphere_mcplot_command_lc(self):
        model = "linearcardioid"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_c(self):
        model = "cardioid"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_wc(self):
        model = "wrappedcauchy"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_vm(self):
        model = "vonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dvm(self):
        model = "dvonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_sevm(self):
        model = "sevonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_invsevm(self):
        model = "invsevonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_miaec(self):
        model = "miaecardioid"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_miaewc(self):
        model = "miaewrappedcauchy"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_miaevm(self):
        model = "miaevonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_miaejp(self):
        model = "miaejonespewsey"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_invmiaec(self):
        model = "invmiaecardioid"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_invmiaewc(self):
        model = "invmiaewrappedcauchy"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_invmiaevm(self):
        model = "invmiaevonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_invmiaejp(self):
        model = "invmiaejonespewsey"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    # Test for sampling result of multimodal model
    def test_sphere_mcplot_command_vm_K2(self):
        model = "vonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1 -nmix 2".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dvm_K2(self):
        model = "dvonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1 -nmix 2".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3}".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    # Test for optimizing result
    def test_sphere_mcplot_command_vm_opt(self):
        model = "vonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1 -M optimizing".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3} -M optimizing".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)

    def test_sphere_mcplot_command_dvm_opt(self):
        model = "dvonmises"
        # Generate input data
        argv_str = "{0} {1} -m {2} -si 50 -sw 10 -sc 1 -M optimizing".format(
            self.__input_e,
            self.__input_d,
            model
        )
        argv = argv_str.split()
        args = sphere_estimate.argument_parse(argv)
        sphere_estimate.main(args, SphereMcplotTest.logger)

        # Visualize plot
        argv_str = "{0} {1} {2} 0 -m {3} -M optimizing".format(
            self.__output,
            self.__input_d,
            self.__input_e,
            model
        )
        argv = argv_str.split()
        args = sphere_mcplot.argument_parse(argv)
        sphere_mcplot.main(args, SphereMcplotTest.logger)


if __name__ == '__main__':
    unittest.main()
