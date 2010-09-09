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
from math import sqrt,pi,sin,cos
from geometry import Geometry
from objects import Rod
from copy import copy

def occupancy_radius(occupancy,n,cell_area = 1.0):
	return sqrt(cell_area*occupancy/n/pi)

def wheel(width,height,n,occupancy,separation,material,priority='None'):
	'''Produces a geometry consisting of n rods of specified material
	with specified cell occupancy in a cell of default area = 1.
	The distance between the center of the cell and the body of each
	rod is adjustable using separation.'''
	distance = lambda point1,point2 : sqrt((point1[0]-point2[0])**2+\
	(point1[1]-point2[1])**2)
	r = occupancy_radius(occupancy,n,width*height)
	if n==1 : return Geometry(width,height,[Rod(0,0,material,r)])
	R = r + separation
	wheel_point = lambda N : (R*cos(N*2*pi/n),R*sin(N*2*pi/n))
	wheel_priority = {\
	'None':lambda R,r,n:(R,r,n)\
	,'Occupancy':lambda R,r,n:(distance(wheel_point(0),wheel_point(1))>2*r\
	and (R,r,n) or (r/sin(2*pi/(2*n)),r,n))\
	,'Distance':lambda R,r,n:(distance(wheel_point(0),wheel_point(1))>2*r\
	and (R,r,n) or (R,R*sin(2*pi/(2*n)),n))\
	}
	R,r,n = wheel_priority[priority](R,r,n)
	return Geometry(width,height,[Rod(*wheel_point(N),material=\
	copy(material),radius=r) for N in range(n)])

def max_epsilon(geometry):
	return max(obj.material.epsilon for obj in geometry.objects)
