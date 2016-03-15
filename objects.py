    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the imKvecFormatterplied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
from data import dielectrics, material_names
class Object(object):
    template_str = (
         "\n    (make %(shape)s\n"
         "        (center %(x)s %(y)s %(z)s)\n"
         "        %(other)s\n"
         "        (material %(material)s) )")
    def __init__(self,x,y,z,material,shape='',**others):
        self.x = x
        self.y = y
        self.z = z
        self.material = material
        self.others = others
        self.shape = shape
        
    def __str__(self):
        self.other = '\n        '.join(
            '(%s %s)'%(a,self.others[a]) for a in self.others)
        return Object.template_str % self.__dict__
        del self.other


class Rod(Object):
    
    def __init__(self, x, y, material, radius):
        self.radius = radius
        super(Rod,self).__init__(x, y, 0, material, 'cylinder', 
            height='infinity', radius=radius)


class Block(Object):
    
    def __init__(self, x, y, z, material, size):
        self.size = size
        super(Block, self).__init__(x, y, z, material, 'block',
            size=' '.join([str(s) for s in size]))


class Dielectric(object):
    
    def __init__(self,dielectric):
        if hasattr(dielectric,'__int__'):
            self.epsilon = dielectric
            self.name = 'eps={0:.3f}'.format(dielectric) 
        else:
            self.epsilon = dielectrics[dielectric]
            self.name = material_names[dielectric]
            
    def __str__(self):
        if isinstance(self.epsilon, (list, tuple)) and len(self.epsilon) == 3:
            return ('(make dielectric-anisotropic\n'
                    '            (epsilon-diag '
                    '{0[0]} {0[1]} {0[2]}))'.format(self.epsilon))
        else:            
            return '(make dielectric (epsilon %(epsilon)s))' % self.__dict__
            
    def __repr__(self):
        return 'Dielectric with epsilon = %(epsilon)s' % self.__dict__
