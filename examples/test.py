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

import sys
sys.path.append('../')
from objects import Rod,Dielectric
from simulation import Simulation
from geometry import Geometry
from kspace import KSpace
from utility import wheel,occupancy_radius
from graphics import draw_geometry
from objects import Dielectric,Rod

def basic_test():
	test_geom = Geometry(1,1,[Rod(0,0,Dielectric(11.8), \
	occupancy_radius(0.3,1))])
	sim = Simulation('basic_test', test_geom)
	sim.run_simulation()
	sim.post_process()
	return True
	
def geometry_test():
	test_geom = wheel(1,1,5,0.3,0,Dielectric(11.8)\
	,priority='Occupancy')
	format = 'pdf'
	draw_geometry(test_geom, 'geometry_test',format)
	return True

def bs2d_test():
	test_geom = wheel(1,1,5,0.3,0,Dielectric(11.8)\
	,priority='Occupancy')
	draw_geometry(test_geom,'test_2d')
	test_kspace = KSpace(2,x_res=50,y_res=50)
	sim = Simulation(
            'test_2d', test_geom, test_kspace, numbands=5, resolution=16)
	sim.run_simulation()
	sim.post_process()
	sim.draw_bandstructure_2D(5, filled=False)
	return True
    
tests = [geometry_test, basic_test]
long_tests = [bs2d_test]

results = [test() for test in tests]
