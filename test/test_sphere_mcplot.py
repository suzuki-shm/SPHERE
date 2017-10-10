#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27


import unittest
import os
from sphere import sphere_mcplot


class SphereMcplotTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_mcplot"
        self.__input_d = d_dir + "/input_depth.tsv"
        self.__input_e = d_dir + "/input_estimated.tsv"
        self.__output = d_dir + "/output.png"
        self.__logger = sphere_mcplot.get_logger()

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_dplot_main(self):
        args = {
            "depth_file_path": self.__input_d,
            "estimated_tsv": self.__input_e,
            "output_dest": self.__output,
            "fs": 18,
            "cl": 100
        }
        sphere_mcplot.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
