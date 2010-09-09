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
import objects
from objects import *

class Geometry(object):
	def __init__(self,width,height,objects):
		self.width = width
		self.height = height
		self.objects = objects
		
	def get_area(self):
		return self.width*self.height		
	cell_area = property(get_area)
	
	def get_lattice(self):
		return '(make lattice (size %s %s no-size))' % \
		(self.width,self.height)
	lattice = property(get_lattice)
	
	def get_max_epsilon(self):
		return max(a.material.epsilon for a in self.objects)
	max_epsilon = property(get_max_epsilon)
	
	def __repr__(self):
		object_template = {objects.Rod : '''EdgeForm[Directive[Dashed,Darker[\
Green,0.6]],Disk[{%(x)s,%(y)s},%(radius)s'''}
		
		maxeps = self.max_epsilon
		hue_function = lambda epsilon:(maxeps-epsilon)/epsilon
		return('Graphics[{'+','.join(object_template[a.__class__]\
		%a.__dict__ for a in self.objects)+'}]')
	def __str__(self):
		return '(list' + ''.join(str(a) for a in self.objects) + ')'
	def __iter__(self):
		return self
	#~ def geom_objects(self):
		#~ for obj in self.objects:
			#~ yield 
