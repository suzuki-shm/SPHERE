#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28


import unittest
import os
import numpy as np
from sphere import sphere_symtest
from scipy.stats import norm


class SphereSymtestTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_symtest"
        self.__output = d_dir + "/output.tsv"
        self.__input1 = d_dir + "/input1.tsv"
        self.__input2 = d_dir + "/input2.tsv"
        self.__input = [self.__input1, self.__input2]

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_symtest_main_multi(self):
        args = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
        }
        sphere_symtest.main(args)

    def test_sphere_symtest_main_single(self):
        args = {
            "depth_file_path": [self.__input1],
            "output_dest": self.__output,
        }
        sphere_symtest.main(args)

    def test_sphere_symtest_argument_parse_multi(self):
        argv_str = "{0} {1} {2}".format(self.__output,
                                        self.__input1,
                                        self.__input2)
        argv = argv_str.split()
        args = sphere_symtest.argument_parse(argv)
        args_answer = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_symtest_argument_parse_single(self):
        argv_str = "{0} {1}".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_symtest.argument_parse(argv)
        args_answer = {
            "depth_file_path": [self.__input1],
            "output_dest": self.__output
        }
        self.assertDictEqual(args, args_answer)

    def test_sphere_symtest_command_multi(self):
        argv_str = "{0} {1} {2}".format(self.__output,
                                        self.__input1,
                                        self.__input2)
        argv = argv_str.split()
        args = sphere_symtest.argument_parse(argv)
        sphere_symtest.main(args)

    def test_sphere_symtest_command_single(self):
        argv_str = "{0} {1}".format(self.__output, self.__input1)
        argv = argv_str.split()
        args = sphere_symtest.argument_parse(argv)
        sphere_symtest.main(args)

    def test_perwey_test(self):
        theta = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        depth = np.array([1, 2, 10, 2, 1])
        thetaj = []
        for d, t in zip(depth, theta):
            for r in ([t] * d):
                thetaj.append(r)
        thetaj = np.array(thetaj)
        n_t = np.size(thetaj)
        md = 2.5132741228718345
        mrl = 0.6011271242968684
        b2 = 0
        a2 = np.mean(np.cos(2*(thetaj - md)))
        a3 = np.mean(np.cos(3*(thetaj - md)))
        a4 = np.mean(np.cos(4*(thetaj - md)))
        varb2 = ((1.0-a4)/2.0-2.0*a2+2.0*a2/mrl*(a3+(a2*(1.0-a2))/mrl))/n_t
        z = b2 / np.sqrt(varb2)
        p = 1 - norm.cdf(abs(z))
        p = 2 * np.min([p, 1-p])
        answer = p

        result = sphere_symtest.perwey_test(theta, depth)[1]

        self.assertEqual(result, answer)


if __name__ == '__main__':
    unittest.main()
