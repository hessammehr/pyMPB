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
from pylab import figure,show,linspace,savefig,text,\
griddata,plot,contour,clf,clabel,colorbar
from matplotlib.patches import Ellipse
from numpy import loadtxt
from utility import max_epsilon
from defaults import fig_size,contour_lines,contour_filled,contour_plain,\
colorbar_style
import objects

def draw_geometry(geometry,jobname,format='pdf'):
	global maxeps
	clf()
	maxeps = max_epsilon(geometry)
	drawing_dict = {objects.Rod : draw_rod}
	fig = figure(figsize=fig_size)
	ax = fig.add_subplot(111,aspect='equal')
	X = geometry.width/2
	Y = geometry.height/2
	ax.set_xlim(-X,X)
	ax.set_ylim(-Y,Y)
	
	graphic_objs = [drawing_dict[obj.__class__](index,obj) for index,obj \
	in enumerate(geometry.objects)]
	for graph in graphic_objs:
		for elem in graph:
			ax.add_artist(elem)
	show()
	savefig(jobname+'.'+format,format=format,transparent=True)

def draw_rod(index,rod):
	return (Ellipse((rod.x,rod.y),rod.radius*2,rod.radius*2,alpha=\
	rod.material.epsilon/maxeps,lw=2.0,ls='dashed',ec='darkblue',fc=\
	'lightblue'),text(rod.x,rod.y,'#%s\n$\epsilon=%s$'%\
	(index,rod.material.epsilon),
	ha='center',va='center',family='sans-serif'))

def draw_bandstructure(jobname,kspace,band,ext='.csv',format='pdf',\
filled=True,levels=15,lines=False,labeled=False,legend=False):
	#clf()
	fig = figure(figsize=fig_size)
	ax = fig.add_subplot(111,aspect='equal')
	x,y,z = loadtxt(jobname+ext,delimiter=', ',skiprows=1,usecols=\
	(1,2,4+band),unpack=True)
	if kspace.dimensions == 1:
		pylab.plot(x,y,z)
	elif kspace.dimensions ==2:
		xi = linspace(-0.5,0.5,kspace.x_res)
		yi = linspace(-0.5,0.5,kspace.y_res)
		zi = griddata(x,y,z,xi,yi)
		if filled:
			cs = ax.contourf(xi,yi,zi,levels,**contour_filled)
			legend and colorbar(cs,**colorbar_style)
			cs = lines and ax.contour(xi,yi,zi,levels,**contour_lines)
			labeled and lines and clabel(cs,fontsize=8,inline=1)
		else:
			cs = ax.contour(xi,yi,zi,levels,**contour_plain)
			legend and colorbar(cs,**colorbar_style)
			labeled and clabel(cs,fontsize=8,inline=1)	
		ax.set_xlim(-0.5,0.5)
		ax.set_ylim(-0.5,0.5)
	savefig(jobname+format,format=format,transparent=True)
