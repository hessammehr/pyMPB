    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
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

from __future__ import division, print_function
from os import path, environ, remove, rename, mkdir
import sys
from shutil import rmtree
import subprocess as sp
import re
import numpy as np
import defaults
import graphics
from datetime import datetime
import time
from glob import glob1
from utility import distribute_pattern_images
from kspace import KSpaceRectangular
import log

class Simulation(object): 
    def __init__(
            self, jobname, geometry, kspace=KSpaceRectangular(),
            resolution=defaults.default_resolution,
            mesh_size=defaults.default_mesh_size,
            numbands=defaults.default_numbands,
            initcode=defaults.default_initcode,
            runcode=defaults.default_runcode,
            postcode=defaults.default_postcode,
            work_in_subfolder=True, clear_subfolder=True,
            logger=True, quiet=defaults.isQuiet):
        """Create a simulation object with all parameters describing the
        simulation, including a unique jobname (all generated filenames
        will include this name), the geometry (pyMPB Geometry object),
        the kspace (a pyMPB KSpace object), the resolution, mesh_size
        (see MPB docs), number of bands to calculate and some optional
        strings with Scheme code which will be added to the MPB .ctl
        file as initialization code (initcode), as run commands
        (runcode) and as code executed after the simulation (postcode).

        If work_in_subfolder is True (default), all simulation and log
        output will be placed in a separate subdirectory under the
        current working directory called like the jobname.
        work_in_subfolder can also be a custom subfolder name,
        including any path, if not to be placed in current working
        directory.

        Set clear_subfolder to True (default) if you want this
        subfolder to be emptied (will make backup if there is an old
        folder with the same name). clear_subfolder should be False if
        you want to do postprocessing on existing simulation data.

        If logger is True (default), a jobname.log file will be created
        with all pyMPB output and errors. This output will also go to
        stdout if quiet is False. Alternatively, set logger to a
        customized logger (any object with log(level, msg, *args,
        **kwargs) method).

        """

        self.jobname = jobname
        self.geometry = geometry
        self.kspace = kspace
        self.initcode = initcode
        self.postcode = postcode
        self.resolution = resolution
        self.meshsize = mesh_size
        self.numbands = numbands
        self.quiet = quiet
        self.runcode = runcode

        self.work_in_subfolder = work_in_subfolder
        self.clear_subfolder = clear_subfolder
        if isinstance(work_in_subfolder, bool):
            if work_in_subfolder:
                # create default subfolder from jobname:
                self.workingdir = path.abspath(
                    path.join(path.curdir, jobname))
            else:
                # work here, no subfolder:
                self.workingdir = path.abspath(path.curdir)
        else:
            # hopefully a string
            self.workingdir = path.abspath(
                    path.join(path.curdir, work_in_subfolder))

        # the .ctl file that MPB will use:
        self.ctl_file = jobname + '.ctl'
        # a date & time stamp added to log and output filenames:
        dtstamp = ('_{0.tm_year}-{0.tm_mon:02}-{0.tm_mday:02}'
                   '_{0.tm_hour:02}-{0.tm_min:02}-'
                   '{0.tm_sec:02}').format(time.localtime())
        # the output file, where all MPB output will go:
        self.out_file = path.join(self.workingdir, jobname + dtstamp + '.out')
        # a log file, where information from pyMPB will go:
        self.log_file = path.join(self.workingdir, jobname + dtstamp + '.log')
        # the file where MPB usually saves the dielectric:
        self.eps_file =  path.join(self.workingdir, 'epsilon.h5')

        # logger is not setup yet, because the log file might be placed in a
        # subfolder that still needs to be created. But, I want to log that
        # I created a new directory. So make a simple log buffer:
        to_log = []

        to_log.append('Working in directory ' + self.workingdir)
        if self.work_in_subfolder:
            if path.exists(self.workingdir):
                to_log.append('directory exists already: ' + self.workingdir)
                if self.clear_subfolder:
                    # directory exists, make backup
                    backupdir = self.workingdir + '_bak'
                    if path.exists(backupdir):
                        # previous backup exists already, remove old
                        # backup, but keep .log and .out files (they have
                        # unique names):
                        keepers = (glob1(self.workingdir + '_bak', '*.log') +
                                glob1(self.workingdir + '_bak', '*.out'))
                        to_log.append(
                            ('removing existing backup {0}, but keeping {1}'
                            ' old log and output files').format(
                                backupdir, len(keepers)))
                        for f in keepers:
                            rename(path.join(backupdir, f),
                                path.join(self.workingdir, f))
                        rmtree(backupdir)
                        to_log.append(backupdir + ' removed')
                    # rename current (old) dir to backup:
                    rename(self.workingdir, backupdir)
                    to_log.append('existing ' + self.workingdir +
                                  ' renamed to ' + backupdir)
                    # make new empty working directory:
                    mkdir(self.workingdir)
                    to_log.append(
                        'created directory ' + self.workingdir + '\n')
                else:
                    to_log.append('working in existing directory.')
            else:
                # make new empty working directory:
                mkdir(self.workingdir)
                to_log.append('created directory ' + self.workingdir + '\n')

        if logger:
            if hasattr(logger, 'log') and callable(logger.log):
                # a custom logger was given as parameter, use it:
                log.logger = logger
            else:
                # Create the logger. Afterwards, we can also use
                # log.info() etc. in other modules. All status, logging
                # and stderr output will go through this logger (except
                # MPB's output during simulation):
                log.setup_logger(
                    'root.' + self.jobname, self.log_file, self.quiet,
                    redirect_stderr=True)

        # now we can log the stuff from before:
        if to_log:
            log.info('\n' + '\n'.join(to_log))
        del to_log

        # get modes from runcode:
        self.modes = re.findall("\(run[-]?(.*?)[\s\)]", runcode, re.MULTILINE)

        self.number_of_tiles_to_output = defaults.number_of_tiles_to_output

        # In 3D, there are no pure tm or te modes. MPB renames them 
        # automatically to zodd and zeven, respectively. Do the same:
        if self.geometry.is3D:
            for i, mode in enumerate(self.modes):
                if mode == 'te':
                    log.info('In 3D, there is no pure TE mode. '
                             'I will change it to zeven.')
                    self.modes[i] = 'zeven'
                if mode == 'tm':
                    log.info('In 3D, there is no pure TM mode. '
                             'I will change it to zodd.')
                    self.modes[i] = 'zodd'

        log.info("working with modes: " + str(self.modes) + '\n')

        new_environ_dict = {
            'GUILE_WARN_DEPRECATED': 'no',
            # load scheme files also from pyMPB directory (e.g. dos.scm): 
            'GUILE_LOAD_PATH' : path.dirname(path.abspath(graphics.__file__))}
        environ.update(new_environ_dict)
        log.info('added to environment:' + 
                 ''.join(['\n  {0}={1}'.format(key, environ[key]) for key in 
                         new_environ_dict.keys()]))

        log.info(
            'pyMPB Simulation created with following properties:' + 
            ''.join(['\npyMPBprop: {0}={1!r}'.format(key, val) for key, val in
                self.__dict__.items()]) + '\n\n')
        # TODO log all parameters of Simulation object in such a way
        # that it can be recreated exactly.
        # Maybe even add relevant parts of data.py and defaults.py

    def __str__(self):
        temp_dict = self.__dict__.copy()
        temp_dict['geometry'] = ''.join(str(a) for \
         a in self.geometry.objects)
        temp_dict['lattice'] = self.geometry.lattice
        return (defaults.template%temp_dict)

    def write_ctl_file(self, where='./'):
        filename = path.join(where, self.ctl_file)
        log.info("writing ctl file to %s" % filename)
        log.info("### ctl file for reference: ###\n" + 
            str(self) + '\n### end of ctl file ###\n\n')
        with open(filename,'w') as input_file:
            input_file.write(str(self))

    def run_simulation(self, num_processors=2):
        self.write_ctl_file(self.workingdir)

        mpb_call_str = defaults.mpb_call % dict(num_procs=num_processors)

        with open(self.out_file, 'w') as outputFile:
            log.info("Using MPB " + defaults.mpbversion)
            log.info("Running the MPB-computation using the following "
                     "call:\n" +
                " ".join([mpb_call_str, self.ctl_file]))
            log.info("Writing MPB output to %s" % self.out_file)
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
            log.info('MPB simulation is running... To see progress, please '
                'check the output file %s' % self.out_file)
            # run MPB, write output to outputFile:
            # TODO can we also pipe MPB output to stdout, so the user can
            # see progress?
            p = sp.Popen(mpb_call_str.split() + [self.ctl_file], 
                               stdout=outputFile,
                               stderr=sp.STDOUT,
                               cwd=self.workingdir)
            retcode = p.wait()
            endtime = datetime.now()
            outputFile.write("finished on: %s (duration: %s)\n" % 
                             (str(endtime), str(endtime - starttime)))
            outputFile.write("returncode: " + str(retcode))
            log.info("Simulation finished, returncode: " + str(retcode))

        return retcode

    def epsilon_to_png(self):
        """Convert epsilon.h5 to epsilon.png. """

        if not path.isfile(self.eps_file):
            log.info('epsilon file {0} does not exist, '
                'will not create epsilon PNG.'.format(self.eps_file))
            return

        # make rectangular cell etc:
        callstr = defaults.mpbdata_call % dict(
                    self.__dict__, 
                    h5_file=self.eps_file + ':data', 
                    output_file=defaults.temporary_epsh5)
        log.info("calling: {0}".format(callstr))
        
        try:
            sp.call(callstr.split(), cwd=self.workingdir)
        except OSError as err:
            log.warning('Command could not be executed. Will continue ' +
                        'without converting epsilon.h5 to png.'
                        '\n\tOSError message: {}\n'.format(err))
            return 1
        
        # no error, continue:
        dct = dict(self.__dict__, h5_file=defaults.temporary_epsh5)
        # save dielectric to png:
        if self.geometry.is3D:
            callstr = (defaults.epsh5topng_call_3D % dct,
                       defaults.epsh5topng_call_3D_cross_sect % dct)
        else:
            callstr = (defaults.epsh5topng_call_2D % dct,)
        retcode = 0
        for s in callstr:
            if not retcode:
                log.info("calling: {0}".format(s))
                retcode = retcode or sp.call(s.split(), 
                                             cwd=self.workingdir)

        return retcode

    def fieldpatterns_to_png(self):
        """Convert all field patterns (saved during simulation in h5-files)
        to png-files. Move them to subdirectories. Move the h5-files to the
        subdirectory 'patterns_h5~'. epsilon_to_png must be called before!

        """
        # make list of all field pattern h5 files:
        filenames = glob1(self.workingdir, "*.h5")
        for exclude in ["epsilon.h5", defaults.temporary_epsh5, 'foo']:
            d, f = path.split(path.join(self.workingdir, exclude))
            to_remove = glob1(d, f)
            for fname in to_remove:
                filenames.remove(fname)
        if not filenames:
            return 0

        if not defaults.delete_h5_after_postprocessing:
            # prepare temporary folder:
            # (the h5 files will be moved here after conversion)
            if not path.isdir(path.join(
                        self.workingdir,
                        defaults.temporary_h5_folder)):
                log.info(
                    "creating subdirectory: " +
                    defaults.temporary_h5_folder)
                mkdir(path.join(
                    self.workingdir,
                    defaults.temporary_h5_folder))

        log.info("Will now convert following files to png: %s" % filenames)
        log.info("On all these files, mpb-data will be called with the "
                 "command:")
        log.info(defaults.mpbdata_call % dict(
            self.__dict__,
            output_file=defaults.temporary_h5,
            h5_file='<file.h5>'))
        # NOTE: it would be really nice if we could use h5topng's -R
        # parameter also for multiple datasets, i.e. call h5topng on
        # multiple datasets (rather than on multiple files) and have it
        # automatically append the dataset name to each png file name.
        # That way, if all vector components of a field are exported
        # in one file, the exported components' absolute values will be
        # comparable.
        # Possible workaround: make multiple copies of the h5 file
        # from mpb-data (one cannot specify -d in mpb-data, see below)
        # with filenames denoting the different components, then call
        # h5topng with -R and multiple h5 files while specifying
        # :dataset for each one.
        # -> not very nice. Maybe it is time to implement a Python
        # function to directly read and export the h5 files.
        log.info("and then 4 times h5topng for each field component:")
        for ri in ['.r', '.i']:
            dct = dict(self.__dict__, 
               h5_file=defaults.temporary_h5 + ':<xyz>' + ri,
               eps_file=defaults.temporary_epsh5,
               output_file='pngs_<mode>/<filename>[.<xyz>]' + ri + '.png',
               output_file_no_ovl='pngs_<mode>_no_ovl/<filename>[.<xyz>]' +
                   ri + '.png')
            if self.geometry.is3D:
                log.info(defaults.fieldh5topng_call_3D % dct)
                log.info(defaults.fieldh5topng_call_3D_no_ovl % dct)
            else:
                log.info(defaults.fieldh5topng_call_2D % dct)
                log.info(defaults.fieldh5topng_call_2D_no_ovl % dct)
        if defaults.delete_h5_after_postprocessing:
            log.info("and finally delete the h5 file.")
        else:
            log.info("and finally move the h5 file to temporary folder " +
                     defaults.temporary_h5_folder)

        # Build the regular expression pattern for parsing filenames:

        # re that matches the output, i.e. field (e, d or h) or 'dpwr' etc.:
        f = r'(?P<field>[edh]|hpwr|dpwr)'
        # re that matches the k number part, starting with '.':
        k = r'[.]k\d+'
        # re that matches the band number part, starting with '.':
        b = r'[.]b\d+'
        # re that matches the field component (.x, .y or .z) or nothing:
        c = '(?:[.](?P<comp>[xyz]))?'
        # re that matches anything following '.', which does not contain
        # another period (this should be the mode: te, tm, zodd etc.):
        m = '(?:[.](?P<mode>[^.]+))?'
        # The final re pattern matches MPB output hdf-filenames:
        retest = re.compile(
            ''.join(['(?P<filenamebase>', f, k, b, ')', c, m, '.h5']))

        for fname in filenames:
            # parse the filename to get mode and component(s):
            m = retest.match(fname)
            if m is None:
                # this is strange
                log.warning('Convert field patterns to png: Could not parse '
                            'the file name: {0}'.format(fname))
                continue
            redict = m.groupdict()
            datasets = [redict.get('comp', None)]
            if redict['field'] in 'edh':
                if datasets[0] is None:
                    datasets = ['x', 'y', 'z']
                datasets = [ds + ri for ri in ['.r', '.i'] for ds in datasets]
            else:
                datasets = ['data']
            if redict.get('mode', None) is not None:
                # use mode name in folder name:
                mode = redict['mode']
                foldername = defaults.field_output_folder_prefix + '_' + mode
            else:
                mode = ''
                foldername = defaults.field_output_folder_prefix
            foldername_no_ovl = foldername + '_no_ovl/'
            if not path.isdir(path.join(self.workingdir, foldername)):
                log.info("creating subdirectory: " + foldername)
                mkdir(path.join(self.workingdir, foldername))
            if not path.isdir(path.join(self.workingdir, foldername_no_ovl)):
                log.info("creating subdirectory: " + foldername_no_ovl)
                mkdir(path.join(self.workingdir, foldername_no_ovl))

            # make rectangular cell etc:
            callstr = defaults.mpbdata_call % dict(
                        self.__dict__, 
                        h5_file=fname, 
                        output_file=defaults.temporary_h5)
            # note: never include :dataset here (or as -d),
            # because then mpb-data will only see the real
            # or imaginary part, not both, and will not 
            # properly apply the exponential phase shift
            # if multiple tiles are exported.
            log.debug("calling: {0}".format(callstr))
            if not sp.call(callstr.split(), cwd=self.workingdir):
                #log.debug("success")
                # no error, continue:
                # show some progress:
                print('.', end='')
                sys.stdout.flush()
                retcode = 0
                for dataset in datasets:
                    # make png file name:
                    if mode:
                        png_fname = '.'.join(
                            [redict['filenamebase'], dataset, mode, 'png'])
                    else:
                        png_fname = '.'.join(
                            [redict['filenamebase'], dataset, 'png'])
                    dct = dict(
                        self.__dict__,
                        h5_file=defaults.temporary_h5 + ':' + dataset,
                        eps_file=defaults.temporary_epsh5,
                        output_file=path.join(foldername, png_fname),
                        output_file_no_ovl=path.join(
                            foldername_no_ovl, png_fname))
                    # save mode pattern to png:
                    if self.geometry.is3D:
                        callstr = (defaults.fieldh5topng_call_3D % dct,
                                   defaults.fieldh5topng_call_3D_no_ovl % dct)
                    else:
                        callstr = (defaults.fieldh5topng_call_2D % dct,
                                   defaults.fieldh5topng_call_2D_no_ovl % dct)

                    retcode = 0
                    for s in callstr:
                        if not retcode:
                            log.debug("calling: {0}".format(s))
                            retcode = retcode or sp.call(
                                s.split(),
                                cwd=self.workingdir)
                            if retcode:
                                log.error('error calling {0}'.format(s))
            else:
                log.error('error calling {0}'.format(callstr))
                return 1

            if not retcode:
                if defaults.delete_h5_after_postprocessing:
                    remove(path.join(self.workingdir, fname))
                    log.debug('deleted {0}'.format(fname))
                else:
                    # move h5 file to temporary folder:
                    rename(path.join(self.workingdir, fname), 
                           path.join(
                               self.workingdir,
                               defaults.temporary_h5_folder,
                               fname))
        return 0


    def _export_data_helper(self, output_buffer, dataname):
        """grep for *dataname* in  *output_buffer* and save the data following
        it to a .csv file.

        *output_buffer* is a multiline string returned from outputfile.read(),
        where outputfile is the file object of the MPB output opened in read
        mode.

        *dataname* is something like 'tefreqs', 'zevenfreqs', 'freqs',
        'yparity' or anything else that can be found in the output file, where
        the line starts with *dataname*, is followed by ':' and the data to be
        exported in the same line.

        """
        pp_file_name = '{0}_{1}.csv'.format(self.jobname, dataname)
        pattern = r'^{0}:, (.+)'.format(dataname.lower())
        output_lines = [
            x + '\n' for x in re.findall(pattern, output_buffer, re.MULTILINE)]
        if output_lines:
            log.info("writing {0} data to {1}".format(dataname, pp_file_name))
            with open(
                    path.join(
                        self.workingdir, pp_file_name), 'w') as pp_file:
                pp_file.writelines(output_lines)
        else:
            log.info("No {0} data found in output".format(dataname))


    def post_process(
            self, convert_field_patterns=True, project_bands_list=None):
        """Make csv files for all band information. Make png of epsilon
        file.
        :param convert_field_patterns: If True, also make pngs of all
        fields exported during simulation.
        :param project_bands_list: a list of simulation folders
        (strings), with previously run simulations containing the bands
        to be projected. The list must have exactly one entry for each
        k-vector of the current simulation, or only one entry if the gap
        to plot stays the same for all k-vectors.

        In both these cases, (and when the list is empty), a
        jobname_projected.csv file will be created (empty if empty
        list), which will be read when the bands are plotted in
        draw_bands.

        Leave this None if there are no bands to be projected, then no
        jobname_projected.csv file will be created (in this case the
        bands plot made in draw_bands will contain automatically found
        band gaps.

        :return: None
        """
        #
        try:
            output_file = open(self.out_file,'r')
        except IOError:
            # Could not open output file. This is normal if the
            # simulation was run earlier and now only the
            # postprocessing needs to be done, in which case
            # self.out_file has a different timestamp in the filename
            # than the output file from the earlier simulation run, on
            # which the postprocessing must be done.
            # So just look for the latest .out file in the folder:
            canditates = glob1(self.workingdir, self.jobname + '*.out')
            if len(canditates) > 0:
                # make newest first:
                canditates.sort(reverse=True)
                self.out_file = path.join(self.workingdir, canditates[0])
                log.info('Post-processing output file from previous '
                         'simulation run: {0}'.format(self.out_file))
                output_file=open(self.out_file, 'r')
            else:
                log.exception('Cannot post-process, no simulation output '
                              'file found!')
                return
        output_buffer = output_file.read()
        output_file.close()
        for mode in self.modes:
            if mode:
                log.info("post-processing mode: {0}".format(mode))
            else:
                log.info(
                    'post-processing '
                    '(unrestricted modes, simulated with (run))')

            # export frequencies: (for band diagrams)
            # try to export all possible data:
            datanames = ['freqs', 'velocity', 'dos', 'yparity', 'zparity']
            for dataname in datanames:
                self._export_data_helper(output_buffer, mode + dataname)

            # Save band frequency ranges to csv, from the just generated
            # freqs.csv. Needed e.g. if these bands are going to be
            # projected in another simulation.
            fnamebase = path.join(
                self.workingdir,
                '{0}_{1}{{0}}.csv'.format(self.jobname, mode))
            data = np.loadtxt(
                fnamebase.format('freqs'), delimiter=',', skiprows=1)
            assert (self.numbands == data.shape[1] - 5)
            bandsmax = np.amax(data[:, 5:], axis=0)
            bandsmin = np.amin(data[:, 5:], axis=0)
            # format is %.6f, because MPB only outputs so many digits:
            np.savetxt(
                fnamebase.format('_ranges'),
                np.array(
                    [np.arange(1, self.numbands + 1),
                     bandsmin,
                     bandsmax
                     ]).transpose(),
                header='bandnum, min, max',
                fmt=['%.0f', '%.6f', '%.6f'],
                delimiter=', ')

            # if project_bands_list is supplied, a csv with the continuum
            # band ranges is created:
            if project_bands_list is not None:
                if (len(project_bands_list) > 1 and
                    self.kspace.count_interpolated() !=
                            len(project_bands_list)):
                    log.warning(
                        'project_bands_list supplied to '
                         'Simulation.postprocess does not have the same '
                         'amount of entries than there are k-vectors in '
                         'this simulation.')
                else:
                    # load all ranges files:
                    ranges = []
                    # minimum amount of bands all simulations share:
                    numbands = float('inf')
                    for folder in project_bands_list:
                        jobname = path.basename(path.normpath(folder))
                        filename = path.join(
                            folder,
                            jobname + '_' + mode + '_ranges.csv')
                        try:
                            rng = np.loadtxt(filename, delimiter=',', ndmin=2)
                            if rng.shape[1] == 3:
                                # drop band numbers:
                                rng = rng[:, 1:]
                            if rng.shape[1] != 2:
                                log.warning(
                                    'file "{0}" is malformed.'.format(
                                        filename) +
                                    'Will not handle projected bands.'
                                )
                                break
                            ranges.append(rng)
                            numbands = min(rng.shape[0], numbands)
                        except IOError:
                            # file not found
                            log.warning(
                                'entry "{0}" in project_bands_list supplied '
                                'to Simulation.postprocess does not exist. '
                                'Will not handle projected bands.'.format(
                                    filename
                                )
                            )
                            break

                    # only continue if we did not abort because we could
                    # not load a file:
                    if len(ranges) == len(project_bands_list):
                        # make all the same size and flat (alternating
                        # min/max):
                        for i in range(len(ranges)):
                            ranges[i] = (ranges[i][:numbands]).flatten()
                        if len(ranges) == 0:
                            # No projected bands supplied. This is OK if we
                            # don't want to plot any band gaps or projected
                            # continuum bands. In this case, we write an
                            # empty _projected file:
                            np.savetxt(
                                fnamebase.format('_projected'),
                                [],
                                header='no projected bands')
                        else:
                            contibands = np.empty(
                                (len(ranges), 1 + 2 * numbands))
                            contibands[:, 0] = np.arange(1, len(ranges) + 1)
                            contibands[:, 1:] = np.array(ranges)
                            # format is %.6f, because MPB only outputs so
                            # many digits:
                            np.savetxt(
                                fnamebase.format('_projected'),
                                contibands,
                                header=', '.join(
                                    ['knum'] +
                                    ['band{0} {1}'.format(i, m)
                                     for i in range(1, numbands + 1)
                                     for m in ['min', 'max']]),
                                fmt=['%.0f'] + ['%.6f'] * 2 * numbands,
                                delimiter=', ')

        if not path.exists(self.eps_file) and path.isfile(self.eps_file + '~'):
            # The epsilon.h5 file was renamed before to mark it as temporary.
            # Name it back, otherwise h5topng can't handle the file:
            rename(self.eps_file + '~', self.eps_file)
            log.info("renamed {0} to {1}".format(
                        self.eps_file + '~', self.eps_file)) 

        if not self.epsilon_to_png() == 1:
            if convert_field_patterns:
                self.fieldpatterns_to_png()
    
            # delete temporary files:
            if path.isfile(path.join(self.workingdir, 
                                     defaults.temporary_epsh5)):
                remove(path.join(self.workingdir, defaults.temporary_epsh5))
            if path.isfile(path.join(self.workingdir, defaults.temporary_h5)):
                remove(path.join(self.workingdir, defaults.temporary_h5))
            if not defaults.delete_h5_after_postprocessing:
                # rename the epsilon.h5 file so my system knows it is a
                # temporary file:
                if path.isfile(self.eps_file + '~'):
                    # but delete old temporary file first:
                    remove(self.eps_file + '~')
                if path.isfile(self.eps_file):
                    rename(self.eps_file, self.eps_file + '~')
                    log.info("renamed {0} to {1}".format(
                        self.eps_file, self.eps_file + '~'))
            else:
                # delete all .h5 files:
                if path.isfile(self.eps_file):
                    remove(self.eps_file)
                    log.info("deleted {0}".format(self.eps_file))
        return

    def display_epsilon(self):
        if not path.isfile(path.join(self.workingdir, 'epsilon.png')):
            return
        if self.geometry.is3D:
            files = ['epsilon.png', 'epsilonslab.png']
        else:
            files = ['epsilon.png']
        call = defaults.display_png_call % {'files': ' '.join(files)}
        sp.call(call.split(), cwd=self.workingdir)

    def draw_bandstructure_2D(
            self, band, mode=None, filled=True, levels=15, lines=False,
            labeled=False, legend=False):
        """Draw 2D band contour map of one band"""
        jobname = path.join(self.workingdir, self.jobname)
        if mode is None:
            # default. Draw all modes:
            modes = self.modes
        elif isinstance(mode, (tuple, list)):
            modes = mode
        else:
            modes = [mode]
        for mode in modes:
            graphics.draw_bandstructure_2D(
                jobname, mode, self.kspace, band, filled=filled, levels=levels,
                lines=lines, labeled=labeled, legend=legend)


    def draw_bands(
            self, title='', crop_y=True,
            x_axis_hint=defaults.default_x_axis_hint,
            show=False, block=True, save=True,
            add_epsilon_as_inset=False,
            color_by_parity=False):
        """Plot dispersion relation of all bands calculated along all
        k vectors.

        *x_axis_formatter* is an object with the method
        'apply_to_axis(axis, **kwargs)' which sets the x-axis' tickpositions,
        major ticklabels and label.
        *title* is the subplot's title.
        If *crop_y* is true (default), the y-axis (frequency) will be limited
        so that only frequency values are shown where all bands are known.
        Alternatively, a numeric value of *crop_y* denotes the upper frequency
        value where the plot will be cropped, or if *crop_y* is a
        2-tuple, it denotes the minimum and maximum y-value.

        *x_axis_hint* gives a hint on which kind of ticks and labels should be
        shown on the x-axis and provides the data needed.
        *x_axis_hint* can be one of the following:
        -- integer number: The axis' labels will be the 3D k-vectors. The
               number denotes the number of major ticks and labels distributed
               on the axis.
        -- [integer, format-string]: Same as above, but the labels are
               formatted with the format-string - this gives the possibility
               to only show one of the three vector components, e.g. the
               string "{2}" to only show the k-vector's z-component. The axis
               title will be inferred from the format-string.
        -- KSpace object: This must be a KSpace object created with
               point_labels. These labels usually denote the high symmetry or
               crititical points, and they will be shown on the axis.
        -- CustomAxisFormatter object: This gives the possibility to completely
               customize the x-axis' tick positions, tick labels and axis
               label. If the CustomAxisFormatter's hover data have not been
               set, it will be set here with the k-vectors read from the
               .csv file.

        If *show*, the plot is shown in a window. In this case, set *block*
        to True if the script should wait, otherwise the script might end
        and close the figure. Set *block* to False if you want to display
        other figures.
        If *save* the figure is saved to 'png' and 'pdf' files.

        The band data is loaded from previously saved .csv files, usually
        done in post_process().

        If a .csv file with projected band data exists, projected bands
        will be plotted, band gaps not.

        :param add_epsilon_as_inset: epsilon,png will be added as inset. See
        defaults.py for parameters like size and location.
        :param color_by_parity: Specify 'y' or 'z' to color the plot
        lines with the data taken from the parity files
        <jobname>_<mode>[z/y]parity.csv.

        TODO: add subplots for each file in argument 'comparison_files=[]'

        """
        jobname = path.join(self.workingdir, self.jobname)
        # see if projected bands were calculated:
        projected = False not in [
            path.isfile(jobname + '_{0}_projected.csv'.format(mode)) for
            mode in self.modes]
        # draw data with matplotlib in one subplot:
        plotter = graphics.draw_bands(
            jobname,
            self.modes,
            x_axis_hint=x_axis_hint,
            title=title,
            crop_y=crop_y,
            band_gaps=not projected,
            light_cone=(
                self.geometry.substrate_index if self.geometry.is3D
                else False),
            projected_bands=projected,
            add_epsilon_as_inset=add_epsilon_as_inset,
            color_by_parity=color_by_parity
        )
        # use returned plotter to add to figure:
        #graphics.draw_dos(jobname, self.modes, custom_plotter=plotter)

        if save:
            filename = jobname + '_bands.pdf'
            log.info('saving band diagram to file %s' % filename)
            plotter.savefig(
                filename, transparent=True,
                bbox_inches='tight', pad_inches=0)
            filename = jobname + '_bands.png'
            log.info('saving band diagram to file %s' % filename)
            plotter.savefig(
                filename, transparent=False,
                bbox_inches='tight', pad_inches=0)

        if show:
            plotter.show(block=block)
        else:
            # I don't want the old figures to show later
            del plotter

    def draw_field_patterns(
            self, title='', only_k=None, show=False,
            filetype='pdf'):
        """ Place all field pattern pngs in one diagram and save it to file.
        If only_k is None (default) all found images at all k-vec
        numbers will be added. Specify a tuple (from, to) to only
        include those (indexes into the list of found k-vecs,
        inclusive), OR a list of indexes (also indexes of found
        k-vectors, not all k-vectors simulated!)
        Specify show to also show the figure (script will not block) or
        show='block' to show and block.

        """
        for mode in self.modes:
            dstfile_prefix = path.join(
                self.workingdir, self.jobname + '_patterns')
            dirname = path.join(
                self.workingdir,
                '{0}_{1}'.format(
                    defaults.field_output_folder_prefix, mode) if mode else
                    defaults.field_output_folder_prefix
            )
            distribute_pattern_images(
                imgfolder=dirname,
                dstfile_prefix=dstfile_prefix,
                dstfile_type=filetype,
                vertical_complex_pairs=defaults.
                    field_dist_vertical_cmplx_comps,
                only_k=only_k,
                title=title,
                show=show)
