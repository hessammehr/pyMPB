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
from data import dielectrics
class Object(object):
	template_str = '''(make %(shape)s (center %(x)s %(y)s %(z)s) (height %(height)s)\
(material %(material)s) %(other)s)'''
	def __init__(self,x,y,z,material,height,shape='',**others):
		self.x = x
		self.y = y
		self.z = z
		self.material = material
		self.height = height
		self.others = others
		self.shape = shape
		
	def __str__(self):
		self.other = ''.join('(%s %s)'%(a,self.others[a]) for a in self.others)
		return Object.template_str % self.__dict__
		del self.other

class Rod(Object):
	def __init__(self,x,y,material,radius):
		self.radius = radius
		super(Rod,self).__init__(x,y,0,material,'infinity','cylinder',\
		radius=radius)

class Dielectric(object):
	def __init__(self,dielectric):
		if hasattr(dielectric,'__int__'):
			self.epsilon = dielectric
		else:
			self.epsilon = dielectrics[dielectric]
	def __str__(self):
		return '(make dielectric (epsilon %(epsilon)s))' % self.__dict__
	def __repr__(self):
		return 'Dielectric with epsilon = %(epsilon)s' % self.__dict__
