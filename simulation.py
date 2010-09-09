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
from os import path,system,environ,unlink
from tempfile import NamedTemporaryFile as ntempfile
from re import findall,sub,search
from math import pi,sqrt,sin,cos
import re
from defaults import *
from pygments import highlight
from pygments.lexers import SchemeLexer
from pygments.formatters import TerminalFormatter,HtmlFormatter, \
Terminal256Formatter
from graphics import draw_bandstructure as draw_bs

default_output_pattern = 'output%s'

class Simulation(object): 
	def __init__(self,jobname,geometry,kspace=default_kspace \
	,initcode=default_initcode,mode=default_mode,resolution= \
	default_resolution,numbands=default_numbands,quiet=isQuiet, \
	log_format = 'terminal'):
		self.jobname = jobname
		self.geometry = geometry
		self.kspace = kspace
		self.initcode = initcode
		self.mode = mode
		self.resolution = resolution
		self.numbands = numbands
		self.quiet = quiet
		self.ctl_file = jobname + '.ctl'
		self.out_file = jobname + '.out'
		self.pp_file = jobname + '.csv'
		self.h5_file =  'epsilon.h5'
		self.log_format = log_format
	def __str__(self):
		temp_dict = self.__dict__.copy()
		temp_dict['geometry'] = ''.join(str(a) for \
		 a in self.geometry.objects)
		temp_dict['lattice'] = self.geometry.lattice
		return (template%temp_dict)
	
	def runSimulation(self):
		self.log(str(self))
		input_file = open(self.ctl_file,'w')
		input_file.write(str(self))
		input_file.close()
		environ['GUILE_WARN_DEPRECATED']='no'
		command = 'mpb %s > %s'%(self.ctl_file,self.out_file)
		self.log(command)
		system(command)
		unlink(self.ctl_file)
	
	def postProcess(self):
		output_file = open(self.out_file,'r')
		pp_buffer = output_file.read()
		output_file.close()
		unlink(self.out_file)
		output_file = open(self.pp_file,'w')
		pattern = self.mode == 'tm' and r'''tmfreqs:, (.+)''' \
		or r'''tefreqs:, (.+)'''
		output_lines = [x+'\n' for x in findall(pattern,pp_buffer,re.M)]
		output_file.writelines(output_lines)
		output_file.close()
		unlink(self.h5_file)
	def log(self,text):
		if not self.quiet:
			print highlight(text,SchemeLexer(), \
			Terminal256Formatter(style='pastie'))
	def draw_bandstructure(self,band,filled=True,levels=15,lines=False, \
	labeled=False,legend=False):
		draw_bs(self.jobname,self.kspace,band,filled=filled,levels= \
		levels,lines=lines,labeled=labeled,legend=legend)

