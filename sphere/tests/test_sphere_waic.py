#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28


import unittest
import os
from sphere import sphere_waic


class SphereWaicTest(unittest.TestCase):
    def setUp(self):
        d_dir = os.path.dirname(__file__) + "/data/test_sphere_waic"
        self.__input = [d_dir + "/input1.tsv"]
        self.__output = d_dir + "/output.tsv"
        self.__logger = sphere_waic.get_logger()

    def tearDown(self):
        if os.path.exists(self.__output):
            os.remove(self.__output)

    def test_sphere_waic_main(self):
        args = {
            "log_lik_files": self.__input,
            "o": self.__output,
            "t": "bda3"
        }
        sphere_waic.main(args, self.__logger)

    def test_sphere_waic_main_original(self):
        args = {
            "log_lik_files": self.__input,
            "o": self.__output,
            "t": "original"
        }
        sphere_waic.main(args, self.__logger)

    def test_sphere_waic_main_both(self):
        args = {
            "log_lik_files": self.__input,
            "o": self.__output,
            "t": "both"
        }
        sphere_waic.main(args, self.__logger)


if __name__ == '__main__':
    unittest.main()
