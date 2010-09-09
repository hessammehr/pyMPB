	#Copyright 2009 Seyed Hessam Moosavi Mehr
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

class KSpace(object):
	vector3 = '(vector3 %s %s %s)\n'
	def __init__(self,dimensions,**kwargs):
		default_dict = {'resolution':32,'x_res':20,'y_res':20}
		default_dict.update(kwargs)
		self.dimensions = dimensions
		self.__dict__.update(default_dict)
	def __str__(self):
		vectors = ''.join(KSpace.vector3 % (x,y,z) \
		for x,y,z in self.points())
		return '(list %s)'%vectors
	def points(self):
		if self.dimensions==1:
			resolution = self.resolution
			hyp_resol = int(round(resolution * sqrt(2)))
			side1 = [(0.5 * n/resolution,0,0) for n in \
			range(resolution+1)]
			side2 = [(0.5,0.5 * n/resolution,0) for n in \
			range(1,resolution+1)]
			side3 = [(0.5 - 0.5*n/hyp_resol,0.5 - 0.5*n/hyp_resol,0) \
			for n in range(1,hyp_resol+1)]
			return side1 + side2 + side3
		elif self.dimensions==2:
			# Very unpythonic, sorry!
			grid = []
			x_res, y_res = self.x_res, self.y_res
			for y in range(y_res):
				for x in range(x_res):
					grid.append((-0.5+x/(x_res-1),-0.5+y/(y_res-1),0))
			return grid
