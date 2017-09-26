#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26


import unittest
import os
from sphere import depth_plot


class DepthPlotTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_depth_plot"
        self.__input = d_dir + "/input.tsv"
        self.__output = d_dir + "/output.png"
        self.__logger = depth_plot.get_logger()

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_depth_plot_main_np5(self):
        args = {
            "depth_file_path": self.__input,
            "output_file_dest": self.__output,
            "np": 5,
            "fs": 30
        }
        depth_plot.main(args, self.__logger)

    def test_depth_plot_main_np29(self):
        args = {
            "depth_file_path": self.__input,
            "output_file_dest": self.__output,
            "np": 29,
            "fs": 30
        }
        depth_plot.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
