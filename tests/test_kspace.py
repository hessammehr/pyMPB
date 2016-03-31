    #Copyright 2016 Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import unittest

import sys
sys.path.append('../')
import numpy as np
from kspace import KSpace, KSpaceTriangular
from kspace import KSpaceRectangular, KSpaceRectangularGrid

class TestKSpaces(unittest.TestCase):

    def test_KSpaceTriangular_string_representation(self):
        k_interpolation = 8
        target = '\n'.join(['(interpolate %i (list',
                            '    (vector3 0 0 0)',
                            '    (vector3 0 0.5 0)',
                            '    (vector3 (/ -3) (/ 3) 0)',
                            '    (vector3 0 0 0)))']) % k_interpolation
        test_kspace = KSpaceTriangular(k_interpolation=k_interpolation)
        self.assertEqual(str(test_kspace), target)

    def test_KSpaceTriangular_string_representation_no_interpolation(self):
        k_interpolation = 0
        target = '\n'.join(['(list',
                            '    (vector3 0 0 0)',
                            '    (vector3 0 0.5 0)',
                            '    (vector3 (/ -3) (/ 3) 0)',
                            '    (vector3 0 0 0))'])
        test_kspace = KSpaceTriangular(k_interpolation=k_interpolation)
        self.assertEqual(str(test_kspace), target)

    def test_KSpaceTriangular_count_interpolated(self):
        k_interpolation = 10
        target = 34
        test_kspace = KSpaceTriangular(k_interpolation=k_interpolation)
        self.assertEqual(test_kspace.count_interpolated(), target)

    def test_KSpaceRectangular_string_representation(self):
        k_interpolation = 15
        target = '\n'.join(['(interpolate %i (list',
                            '    (vector3 0 0 0)',
                            '    (vector3 0.5 0 0)',
                            '    (vector3 0.5 0.5 0)',
                            '    (vector3 0 0 0)))']) % k_interpolation
        test_kspace = KSpaceRectangular(k_interpolation=k_interpolation)
        self.assertEqual(str(test_kspace), target)

    def test_KSpaceRectangular_string_representation_no_interpolation(self):
        k_interpolation = 0
        target = '\n'.join(['(list',
                            '    (vector3 0 0 0)',
                            '    (vector3 0.5 0 0)',
                            '    (vector3 0.5 0.5 0)',
                            '    (vector3 0 0 0))'])
        test_kspace = KSpaceRectangular(k_interpolation=k_interpolation)
        self.assertEqual(str(test_kspace), target)

    def test_KSpaceRectangular_count_interpolated(self):
        k_interpolation = 16
        target = 3 * 16 + 4
        test_kspace = KSpaceRectangular(k_interpolation=k_interpolation)
        self.assertEqual(test_kspace.count_interpolated(), target)

    def test_KSpaceRectangularGrid_points(self):
        x_steps = 25
        y_steps = 20
        testpoints = []
        for y in np.linspace(-0.5, 0.5, y_steps):
            for x in np.linspace(-0.5, 0.5, x_steps):
                testpoints.append((x, y, 0.0))
        test_kspace = KSpaceRectangularGrid(x_steps=x_steps, y_steps=y_steps)
        self.assertEqual(test_kspace.points(), testpoints)

    def test_KSpaceRectangularGrid_count_interpolated(self):
        """Test that there is not interpolation with KSpaceRectangularGrid."""
        x_steps = 34
        y_steps = 21
        test_kspace = KSpaceRectangularGrid(x_steps=x_steps, y_steps=y_steps)
        self.assertEqual(
            len(test_kspace.points()),
            test_kspace.count_interpolated())

    def test_KSpace_custom_points_list(self):
        input_points = [
            1, 2.2, '(/ 3)',
            [1], (2.2,), ['(/ 4)'], np.array([5.5]),
            (1, 2), [2.2, '3'], np.array([6.6, 7.7]), np.array([1.0, 7.7]),
            (2, 3.3, '4.5'), (0, 0.5, 0), np.array([6.6, 7.7, 8]),
            (1, 2, 3, 4), (1, 2, 3, 4, 5), np.arange(4)]
        corrected_points = [
            (1, 0, 0), (2.2, 0, 0), ('(/ 3)', 0, 0),
            (1, 0, 0), (2.2, 0, 0), ('(/ 4)', 0, 0), (5.5, 0, 0),
            (1, 2, 0), (2.2, '3', 0), (6.6, 7.7, 0), (1.0, 7.7, 0),
            (2, 3.3, '4.5'), (0, 0.5, 0), (6.6, 7.7, 8.0),
            (1, 2, 3), (1, 2, 3), (0, 1, 2)]
        test_kspace = KSpace(points_list=input_points)
        self.assertEqual(test_kspace.points(), corrected_points)

if __name__ == '__main__':
    unittest.main()