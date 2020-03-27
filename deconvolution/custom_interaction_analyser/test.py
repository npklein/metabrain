"""
File:         test.py
Created:      2020/03/26
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""
# Standard imports.
import unittest

# Third party imports.

# Local application imports.
from .src.workers import calc_f_value, get_p_value
from .src.main import get_z_score


class Test(unittest.TestCase):
    """
    Test: unit tester for the programming1.py program.
    """
    @classmethod
    def setUpClass(cls):
        """"
        Method executed at class setup.
        """
        print('setUpClass\n')

    @classmethod
    def tearDownClass(cls):
        """"
        Method executed at class tear down.
        """
        print('tearDownClass\n')

    def setUp(self):
        """
        Method that sets up the test input: runs before each test.
        """
        print('setUp\n')

    def tearDown(self):
        """
        Method that cleans up; runs after each test
        """
        print('tearDown\n')

    def test_calc_f_value(self):
        """
        Test the calculate_ratios method in the main class.
        """
        print('calc_f_value\n')


    def test_get_p_value(self):
        """
        Test the calculate_ratios method in the main class.
        """
        print('get_p_value\n')

    def test_get_z_score(self):
        """
        Test the calculate_ratios method in the main class.
        """
        print('get_z_score\n')


if __name__ == '__main__':
    unittest.main()
