#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28


import unittest
import os
from sphere import sphere_cstats


class SpherecstatsTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_cstats"
        self.__input = [d_dir + "/input.tsv"]
        self.__output = d_dir + "/output.tsv"
        self.__logger = sphere_cstats.get_logger()

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_cstats_main(self):
        args = {
            "depth_file_path": self.__input,
            "output_dest": self.__output,
        }
        sphere_cstats.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
