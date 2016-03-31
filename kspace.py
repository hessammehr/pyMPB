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
from math import sqrt
from numpy import linspace
from defaults import default_k_interpolation
import log

class KSpace(object):
    def __init__(self, points_list, **kwargs):
        """Setup a k-space for the simulation with custom k-points.

        Keyword arguments:
        points_list -- the list of critical k-points:
                       Usually, a list of 3-tuples, but can also be a list of
                       2-tuples or a list of numbers, from which the 3-tuples
                       will be built internally by expanding the missing
                       dimensions with zeros,
                       e.g. [0.5, (0.5, 1)] -> [(0.5, 0, 0), (0.5, 1, 0)].
        k_interpolation -- in the simulation, the points_list will be expanded
                       by linearly interpolating between every consecutive pair
                       of points, by adding k_interpolation points between each
                       pair. (default: defaults.default_k_interpolation)
        """
        # build list of 3-tuples:
        three_list = []
        for item in points_list:
            try:
                length = len(item)
                if hasattr(item, 'isalnum'):
                    # This is a string. That is a single item,
                    # even though it has a length:
                    length = 0
            except TypeError:
                # This item is not a list, tuple or similar.
                length = 0
            if length == 0:
                three_list.append((item, 0, 0))
            elif length == 1:
                three_list.append((item[0], 0, 0))
            elif length == 2:
                three_list.append((item[0], item[1], 0))
            elif length == 3:
                three_list.append(tuple(item))
            else:
                three_list.append(tuple(item[0:3]))
                log.warning(
                    'KSpace: a point has been supplied with length > 3. '
                    'I will only use the first 3 entries.')

        default_dict = {'k_interpolation':default_k_interpolation,
                        'points_list':three_list}
        default_dict.update(kwargs)
        self.__dict__.update(default_dict)

    def __str__(self):
        vector3 = '\n    (vector3 %s %s %s)'
        vectors = ''.join(vector3 % (x,y,z) for x,y,z in self.points())
        if self.k_interpolation:
            return ('(interpolate %i (list%s))' %
            #return ('(kinterpolate-uniform %i (list %s))' %
                    (self.k_interpolation, vectors))
        else:
            return '(list%s)'%vectors

    def count_interpolated(self):
        """Return total number of k-vecs after interpolation."""
        return (len(self.points()) - 1) * (self.k_interpolation + 1) + 1

    def points(self):
        """Return the bare list of k-points, before any k_interpolation is
        applied.

        """
        return self.points_list

class KSpaceTriangular(KSpace):
    def __init__(self, k_interpolation=default_k_interpolation):
        """Setup a k-space for the simulation with critical k-points along the
        boundary of the irreducible brillouin zone of the triangular/hexagonal
        lattice, i.e.: [Gamma, M, K, Gamma].

        """
        KSpace.__init__(self,
            points_list=[(0, 0, 0), (0, 0.5, 0), ('(/ -3)', '(/ 3)', 0),
                         (0, 0, 0)],
            k_interpolation=k_interpolation)

class KSpaceRectangular(KSpace):
    def __init__(self, k_interpolation=default_k_interpolation):
        """Setup a k-space for the simulation with critical k-points along the
        boundary of the irreducible brillouin zone of the rectangular lattice,
        i.e.: [Gamma, X, M, Gamma].

        """
        KSpace.__init__(self,
            points_list=[(0, 0, 0), (0.5, 0, 0), (0.5, 0.5, 0),
                         (0, 0, 0)],
            k_interpolation=k_interpolation)

class KSpaceRectangularGrid(KSpace):
    def __init__(self, x_steps, y_steps):
        """Setup a k-space with k-points distributed on a rectangular grid.

        The k-points are distributed in the k_x-k_y-plane over the smallest
        rectangular brillouin zone of the rectangular lattice, i.e. k_x (k_y)
        varies in x_steps (y_steps) from -0.5 to 0.5 (inclusive), respectively.

        """
        grid = [(x, y, 0.0)
                for y in linspace(-0.5, 0.5, y_steps)
                for x in linspace(-0.5, 0.5, x_steps)]
        # x_steps and y_steps needed in bandstructure plot in graphics.py:
        KSpace.__init__(self, points_list=grid, k_interpolation=0,
                        x_steps=x_steps, y_steps=y_steps)