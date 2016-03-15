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

from __future__ import division, print_function
from os import path, environ, remove, rename, mkdir
import sys
from shutil import rmtree
import subprocess as sp
#from tempfile import NamedTemporaryFile as ntempfile
from re import findall
#from math import pi,sqrt,sin,cos
import re
from defaults import *
#from pygments import highlight
#from pygments.lexers import SchemeLexer
#from pygments.formatters import TerminalFormatter,HtmlFormatter, \
#Terminal256Formatter
import graphics
from datetime import datetime
from glob import glob1
from utility import distribute_pattern_images

default_output_pattern = 'output%s'

class Simulation(object): 
    def __init__(
            self, jobname, geometry, kspace=default_kspace,
            initcode=default_initcode, postcode=default_postcode,
            runcode=default_runcode, resolution=default_resolution, 
            numbands=default_numbands, work_in_subfolder=True,
            quiet=isQuiet, log_format='terminal'):
                
        self.jobname = jobname
        self.geometry = geometry
        self.kspace = kspace
        self.initcode = initcode
        self.postcode = postcode
        self.resolution = resolution
        self.numbands = numbands
        self.quiet = quiet
        self.log_format = log_format
        self.runcode = runcode
        
        self.work_in_subfolder = work_in_subfolder
        if work_in_subfolder:
            self.workingdir = path.abspath(path.join(path.curdir, jobname))
        else:
            self.workingdir = path.abspath(path.curdir)
        self.log('Working in directory ' + self.workingdir)
        
        self.ctl_file = jobname + '.ctl'
        self.out_file = path.join(self.workingdir, jobname + '.out')
        self.eps_file =  path.join(self.workingdir, 'epsilon.h5')

        # get modes from runcode:
        self.modes = re.findall("\(run[-]?(.*?)[\s\)]", runcode, re.MULTILINE)

        # In 3D, there are no pure tm or te modes. MPB renames them 
        # automatically to zodd and zeven, respectively. Do the same:
        if self.geometry.is3D:
            for i, mode in enumerate(self.modes):
                if mode == 'te':
                    self.log('In 3D, there is no pure TE mode. '
                             'I will change it to zeven.')
                    self.modes[i] = 'zeven'
                if mode == 'tm':
                    self.log('In 3D, there is no pure TM mode. '
                             'I will change it to zodd.')
                    self.modes[i] = 'zodd'
                    
        self.log("calculating modes: " + str(self.modes) + '\n')
        self.number_of_tiles_to_output = default_number_of_tiles_to_output

        environ['GUILE_WARN_DEPRECATED'] = 'no'
        # load scheme files also from pyMPB directory (e.g. dos.scm): 
        environ['GUILE_LOAD_PATH'] = \
            path.dirname(path.abspath(graphics.__file__))
            
        # in these fields, we will later (in postProcess) save all results:
        # NOT IMPLEMENTED YET    
        # dictionary (entry for each mode) of numpy-arrays with band data;
        #   shape of arrays (number k vecs, number of bands + 5);
        #   the columns of these arrays are: 
        #   k_index, k_x, k_y, k_z, k_magnitude and frequencies for each band 
        #self.banddata = dict()
        
        # dictionary (entry for each mode) of (numfreqs, 2)-numpy-arrays:
        #     columns: frequency, density of states
        #self.dosdata = dict()
                
    def __str__(self):
        temp_dict = self.__dict__.copy()
        temp_dict['geometry'] = ''.join(str(a) for \
         a in self.geometry.objects)
        temp_dict['lattice'] = self.geometry.lattice
        return (template%temp_dict)
    
    def write_ctl_file(self, where='./'):
        self.log("### ctl file ###")
        self.log(str(self), end='\n\n')
        with open(path.join(where, self.ctl_file),'w') as input_file:
            input_file.write(str(self))
        
    def runSimulation(self, num_processors=2):
        if self.work_in_subfolder:
            if path.exists(self.workingdir):
                # directory exists, make backup
                if path.exists(self.workingdir + '_bak'):
                    # previous backup exists already, remove old backup:
                    rmtree(self.workingdir + '_bak')
                    self.log(self.workingdir + '_bak removed')
                rename(self.workingdir, self.workingdir + '_bak')
                self.log(self.workingdir + ' renamed to ' + 
                        self.workingdir + '_bak')
            # make new working directory:
            mkdir(self.workingdir)
            self.log('created directory ' + self.workingdir + '\n')
            
        self.write_ctl_file(self.workingdir)
        
        mpb_call_str = mpb_call % dict(num_procs=num_processors)        
        
        with open(self.out_file, 'w') as outputFile:
            self.log("Running the MPB-computation using the following call:")
            self.log(" ".join([mpb_call_str, self.ctl_file]))
            # write Time and ctl as reference:     
            outputFile.write("This is a simulation started by pyMPB\n")
            starttime = datetime.now()
            outputFile.write("Date: " + str(starttime) + "\n")
            outputFile.write(
                ["2D-Simulation\n", "3D-Simulation\n"]
                [int(self.geometry.is3D)])
            outputFile.write("\n=================================\n")
            outputFile.write("=========== CTL INPUT ===========\n")
            outputFile.write("=================================\n\n")
            outputFile.write(str(self))
            outputFile.write("\n\n==================================\n")
            outputFile.write("=========== MPB OUTPUT ===========\n")
            outputFile.write("==================================\n\n")
            outputFile.flush()
            # run MPB, write output to outputFile:
            p = sp.Popen(mpb_call_str.split() + [self.ctl_file], 
                               stdout=outputFile,
                               stderr=sp.STDOUT,
                               cwd=self.workingdir)
            retcode = p.wait()
            endtime = datetime.now()
            outputFile.write("finished on: %s (duration: %s)\n" % 
                             (str(endtime), str(endtime - starttime)))
            outputFile.write("returncode: " + str(retcode))
            self.log("Simulation finished, returncode: " + str(retcode))
        return retcode
    
    def epsilon_to_png(self):
        """Convert epsilon.h5 to epsilon.png. """
    
        # make rectangular cell etc:
        callstr = mpbdata_call % dict(
                    self.__dict__, 
                    h5_file=self.eps_file + ':data', 
                    output_file=temporary_epsh5)
        self.log("calling: {0}".format(callstr))
        if not sp.call(callstr.split(), cwd=self.workingdir):
            # no error, continue:
            dct = dict(self.__dict__, h5_file=temporary_epsh5)
            # save dielectric to png:
            if self.geometry.is3D:
                callstr = (epsh5topng_call_3D % dct,
                           epsh5topng_call_3D_cross_sect % dct)
            else:
                callstr = (epsh5topng_call_2D % dct,)
            retcode = 0
            for s in callstr:
                if not retcode:
                    self.log("calling: {0}".format(s))
                    retcode = retcode or sp.call(s.split(), 
                                                 cwd=self.workingdir)
        else:
            return 1
                
        return retcode
        
    def fieldpatterns_to_png(self):
        """Convert all field patterns (saved during simulation in h5-files) to 
        png-files. Move them to subdirectories. Move the h5-files to the
        subdirectory 'patterns_h5~'. epsilon_to_png must be called before!
        
        """
        # prepare temporary folder:
        # (the h5 files will be moved here after conversion)
        if not path.isdir(path.join(self.workingdir, temporary_h5_folder)):
            self.log("creating subdirectory: " + temporary_h5_folder)
            mkdir(path.join(self.workingdir, temporary_h5_folder))

        # make list of all field pattern h5 files:
        filenames = glob1(self.workingdir, "*[edh].*.h5")
        if not filenames:
            return 0
            
        self.log("will now convert following files to png: %s" % filenames)
        self.log("on all these files, mpb-data will be called like so:")
        self.log(mpbdata_call % dict(self.__dict__, output_file=temporary_h5,
                                     h5_file='<file.h5>'))
        self.log("and then 4 times h5topng:")
        for comp in ['.r', '.i']:            
            dct = dict(self.__dict__, 
               h5_file=temporary_h5 + ':' + 
                   default_field_component_to_export + comp,
               eps_file=temporary_epsh5,
               output_file='<mode>/<filename>' + comp + '.png',
               output_file_no_ovl='<mode>_no_ovl/<filename>' + 
                   comp + '.png')
            if self.geometry.is3D:
                self.log(fieldh5topng_call_3D % dct)
                self.log(fieldh5topng_call_3D_no_ovl % dct)
            else:
                self.log(fieldh5topng_call_2D % dct)
                self.log(fieldh5topng_call_2D_no_ovl % dct)
        self.log("and finally move the h5 file to temporary folder " + 
                            temporary_h5_folder)
                
        self.log('|' + "-" * len(filenames) + '|\n|', end='')
        
        # i will later try to extract the mode from file name:
        findmode_re = re.compile(r'.*[.](.+?)[.]h5')
        
        for fname in filenames:
            # use mode name as folder name: 
            # TODO: this does not work if run with (run) only!
            match = findmode_re.match(fname)
            foldername = match.groups()[-1] + '/'
            foldername_no_ovl = match.groups()[-1] + '_no_ovl/'
            if not path.isdir(path.join(self.workingdir, foldername)):
                #self.log("creating subdirectory: " + foldername)
                mkdir(path.join(self.workingdir, foldername))
            if not path.isdir(path.join(self.workingdir, foldername_no_ovl)):
                #self.log("creating subdirectory: " + foldername_no_ovl)
                mkdir(path.join(self.workingdir, foldername_no_ovl))
        
            # make rectangular cell etc:
            callstr = mpbdata_call % dict(
                        self.__dict__, 
                        h5_file=fname, 
                        output_file=temporary_h5)
            # note: never include :dataset here (or as -d),
            # because then mpb-data will only see the real
            # or imaginary part, not both, and will not 
            # properly apply the exponential phase shift
            # if multiple tiles are exported.
            #self.log("calling: {0}".format(callstr))
            if not sp.call(callstr.split(), cwd=self.workingdir):
                # no error, continue:
                # show progress:
                self.log('=', end='')
                sys.stdout.flush()
                for comp in ['.r', '.i']:            
                    dct = dict(self.__dict__, 
                        h5_file=temporary_h5 + ':' + 
                            default_field_component_to_export + comp,
                        eps_file=temporary_epsh5,
                        output_file=foldername + 
                            fname.rstrip('.h5') + comp + '.png',
                        output_file_no_ovl=foldername_no_ovl + 
                            fname.rstrip('.h5') + comp + '.png')
                    # save mode pattern to png:
                    if self.geometry.is3D:
                        callstr = (fieldh5topng_call_3D % dct,
                                   fieldh5topng_call_3D_no_ovl % dct)
                    else:
                        callstr = (fieldh5topng_call_2D % dct,
                                   fieldh5topng_call_2D_no_ovl % dct)
                    
                    retcode = 0
                    for s in callstr:
                        if not retcode:
                            #self.log("calling: {0}".format(s))
                            retcode = retcode or sp.call(s.split(), 
                                                         cwd=self.workingdir)
                if retcode:
                    return retcode
            else:
                return 1
                
            if not retcode:
                    # move h5 file to temporary folder:
                    rename(path.join(self.workingdir, fname), 
                           path.join(
                               self.workingdir, temporary_h5_folder, fname))
        
        self.log('|')      
        return retcode
    
    def postProcess(self, convert_field_patterns=True):
        # make csv files for all band information:
        output_file = open(self.out_file,'r')
        pp_buffer = output_file.read()
        output_file.close()
        for mode in self.modes:
            # export frequencies: (for band diagrams)
            if mode:
                pp_file_name = '{0}_{1}.csv'.format(self.jobname, mode)
                pattern = r'^{0}freqs:, (.+)'.format(mode.lower())
                self.log("postprocessing mode: {0}".format(mode))
            else:
                pp_file_name = '{0}.csv'.format(self.jobname)
                pattern = r'^freqs:, (.+)'
                self.log("postprocessing all modes")

            #pattern = mode == 'tm' and r'tmfreqs:, (.+)' or r'tefreqs:, (.+)'
            self.log("writing data to {0}".format(pp_file_name))
            output_lines = [
                x + '\n' for x in findall(pattern, pp_buffer, re.MULTILINE)]
            with open(
                    path.join(self.workingdir, pp_file_name), 'w') as pp_file:
                pp_file.writelines(output_lines)
                
            # export group velocities (if available):
            if mode:
                pp_file_name = '{0}_velocity_{1}.csv'.format(
                    self.jobname, mode)
                pattern = r'^{0}velocity:, (.+)'.format(mode.lower())
            else:
                pp_file_name = '{0}_velocity.csv'.format(self.jobname)
                pattern = r'^velocity:, (.+)'
            output_lines = [
                x + '\n' for x in findall(pattern, pp_buffer, re.MULTILINE)]
            if output_lines:
                self.log("writing group velocity data to {0}".format(
                    pp_file_name))
                with open(
                        path.join(
                            self.workingdir, pp_file_name), 'w') as pp_file:
                    pp_file.writelines(output_lines)                
                
            # export density of states (if available):    
            if mode:
                pp_file_name = '{0}_dos_{1}.csv'.format(self.jobname, mode)
                pattern = r'^{0}dos:, (.+)'.format(mode.lower())
            else:
                pp_file_name = '{0}_dos.csv'.format(self.jobname)
                pattern = r'^dos:, (.+)'
            output_lines = [
                x + '\n' for x in findall(pattern, pp_buffer, re.MULTILINE)]
            if output_lines:
                self.log("writing dos data to {0}".format(pp_file_name))
                with open(
                        path.join(
                            self.workingdir, pp_file_name), 'w') as pp_file:
                    pp_file.writelines(output_lines)
                    
        # following does not find all Gaps in output, only one for each mode  
        #gaps = findall(r'^(.*)freqs:, .+$(?:\n.*)+?\n(Gap from band (.+)\((.+)\) to band (.+)\((.+)\), (.+)%\n)+', 
        #              pp_buffer, re.MULTILINE)
        #for line in result:
        #    print line
         
        if not path.exists(self.eps_file) and path.exists(self.eps_file + '~'):
            # The epsilon.h5 file was renamed before to mark it as temporary.
            # Name it back, otherwise h5topng can't handle the file:
            rename(self.eps_file + '~', self.eps_file)
            self.log("renamed {0} to {1}".format(
                        self.eps_file + '~', self.eps_file)) 
         
        self.epsilon_to_png()
        if convert_field_patterns:
            self.fieldpatterns_to_png()
        
        # delete temporary files:
        if path.isfile(path.join(self.workingdir, temporary_epsh5)):
            remove(path.join(self.workingdir, temporary_epsh5))       
        if path.isfile(path.join(self.workingdir, temporary_h5)):
            remove(path.join(self.workingdir, temporary_h5))        
        
        # rename the epsilon.h5 file so my system knows it is a temporary file:        
        if path.exists(self.eps_file + '~'):
            # but delete old temporary file first:
            remove(self.eps_file + '~')
        rename(self.eps_file, self.eps_file + '~')
        self.log("renamed {0} to {1}".format(
                    self.eps_file, self.eps_file + '~'))
        
        return
        
    def display_epsilon(self):
        if self.geometry.is3D:
            sp.call(['display', 'epsilon.png', 'epsilonslab.png'], 
                    cwd=self.workingdir)
        else:
            sp.call(['display', 'epsilon.png'], cwd=self.workingdir)     

    def log(self, text, end='\n'):
        if not self.quiet:
            print(text, end=end)
            #print(highlight(text, SchemeLexer(), 
            #                Terminal256Formatter(style='pastie')))
            
##    def draw_bandstructure(
##            self, band, filled=True, levels=15, lines=False, labeled=False,
##            legend=False):
##        """Draw 2D band contour map of one band"""
##        jobname = path.join(self.workingdir, self.jobname)
##        graphics.draw_bandstructure(
##            jobname, self.kspace, band, filled=filled, levels=levels,
##            lines=lines, labeled=labeled, legend=legend)

    def draw_field_patterns(self, title='', only_k_slice=None, show=False):
        """ Place all field pattern pngs in one diagram and save it to file.
        If only_k_slice is None (default) all found images at all k vec 
        numbers will be added. Specify a tuple (from, to) to only include 
        those (indices into the list of found k-vecs, inclusive).
        Specify show to also show the figure (script will not block) or
        show='block' to show and block.
        
        """
        for mode in self.modes:
            dstfile = path.join(self.workingdir, 
                      self.jobname + '_{0}_patterns'.format(mode))
            dirname = path.join(self.workingdir, mode)
            distribute_pattern_images(
                imgfolder=dirname,  
                dstfile=[dstfile + ext for ext in ['.pdf', '.png']], 
                only_k_slice=only_k_slice, 
                title=title, 
                show=show)

    def draw_bands(
            self, title='', comparison_files=[], show=False, block=True,
            crop_y=True):
        """Draw all bands calculated along all k vecs in a band diagram.
        If crop_y is true (default), the y-axis (frequency) will be limited 
        so that only frequency values are shown where all bands are known.
        Alternatively, a numeric value of crop_y denotes the upper frequency
        value where the plot will be cropped.
        
        TODO: add diagrams for each file in comparison_files
        
        """
        jobname = path.join(self.workingdir, self.jobname)
        # draw data with matplotlib in one subplot:
        plotter = graphics.draw_bands(jobname, self.modes, title=title, 
                                      crop_y=crop_y, 
                                      light_cone=self.geometry.is3D)
        # use returned plotter to add to figure:
        #graphics.draw_dos(jobname, self.modes, custom_plotter=plotter)                                      
        
        plotter.savefig(
            jobname + '_bands.pdf', transparent=True, 
            bbox_inches='tight', pad_inches=0)
        plotter.savefig(
            jobname + '_bands.png', transparent=False, 
            bbox_inches='tight', pad_inches=0)
        
        if show:
            plotter.show(block=block)
        else:
            # I don't want the old figures to show later
            del plotter
        
        
