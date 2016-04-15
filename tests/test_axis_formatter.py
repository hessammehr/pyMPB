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
import defaults
import axis_formatter

class TestKSpaces(unittest.TestCase):

    def test_infer_k_axis_label_from_format_string(self):
        testdummies = [
            ('{0}', '$k_x$'),
            ('{}', '$k_x$'),
            ('{0:.{1}f}',  '$k_x$'),
            ('{:.{}f}', '$k_x$'),
            ('{0:{1}.{2}f}',  '$k_x$'),
            ('{:{}.{}f}',  '$k_x$'),
            ('{3:.3f}', r'$|k|$'),
            ('({0:.3f}, {1:.3f})', r'$(k_x, k_y)$'),
            ('({:.3f}, {:.3f})', r'$(k_x, k_y)$'),
            ('({1:.3f}, {2:.3f})', r'$(k_y, k_z)$'),
            ('({0:.3f}, {1:.3f}, {2:.3f})', r'$\vec{k}$'),
            ('({:.3f}, {:.3f}, {:.3f})', r'$\vec{k}$'),
            ('({0:.{3}f}, {1:.{3}f}, {2:.{3}f})',  r'$\vec{k}$'),
            ('({:.{}f}, {:.{}f}, {:.{}f})',  r'$\vec{k}$'),
            ('{{{0:.{3}f}, {1:.{3}f}, {2:.{3}f}}}', '${k_x, k_y, k_z}$'),
            ('{{abc:def}}{0:.{3}f}, {1:.{3}f}, {2:.{3}f}{{qrs:xyz}}',
                '${abc:def}k_x, k_y, k_z{qrs:xyz}$'),
            ('{{{0}', '${k_x$'),
            ('{0}{{', '$k_x{$'),
            ('}}{0}', '$}k_x$'),
            ('{0}}}', '$k_x}$'),
        ]
        for format_str, target in testdummies:
            result = axis_formatter.infer_k_axis_label_from_format_string(
                format_str)
            self.assertEqual(
                result,
                defaults.default_x_axis_label.format(target)
            )


if __name__ == '__main__':
    unittest.main()