# -*- coding:utf-8 -*-
# ----------------------------------------------------------------------
# Copyright 2016 Juergen Probst
#
# This file is part of pyMPB.
#
# pyMPB is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMPB is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyMPB.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------

from __future__ import division

import sys
sys.path.append('../../')
from os import path

import numpy as np
import matplotlib.pyplot as plt

from phc_simulations import TriHoles2D_yWaveguide
from utility import get_gap_bands, sum_of_squares
import log


def main():
    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'
    minstep = 1
    maxstep = 13
    stepsize = 2
    numsteps = int((maxstep - minstep) / stepsize + 1.5)
    steps = np.linspace(minstep, maxstep, num=numsteps, endpoint=True)
    # save previous step's data, needed for comparison with current data:
    prev_step_data = None

    for i, step in enumerate(steps):

        log.info(
            "running simulation with "
            "{0:n} supercell height steps:\n{1}".format(
                numsteps, steps) +
            '\n  ### current step: #{0} ({1}) ###\n'.format(i+1, step))


        ### create and run simulation (according to mode) ###

        ksteps = 17

        # The bigger the unit cell, the more bands fold in below the
        # waveguide band:
        first_wg_band = 3 + 2 * (step - 1)

        sim = TriHoles2D_yWaveguide(
            material='SiN',
            radius=0.380,
            mode='te',
            numbands=first_wg_band + 5,
            k_steps=ksteps,
            supercell_x=step,
            resolution=32,
            mesh_size=7,
            runmode=mode,
            num_processors=8,
            projected_bands_folder='./projected_bands_repo',
            #save_field_patterns_kvecs=[
            #    (0, x, 0) for x in np.linspace(0, 0.5, num=ksteps)],
            #save_field_patterns_bandnums=[11, 12, 13, 22, 23],
            #convert_field_patterns=True
            job_name_suffix='_scx{0:02.0f}'.format(step),
            bands_title_appendix=', supcell width={0:02.0f}'.format(step),
            plot_crop_y=0.6
        )

        if not sim:
            log.error('an error occurred during simulation. See the .out file')
            return


        ### load some data ###

        # load te mode band data:
        fname = path.join(sim.workingdir, sim.jobname + '_tefreqs.csv')
        data = np.loadtxt(fname, delimiter=',', skiprows=1)
        gapbands = get_gap_bands(data[:, 5:], light_line=data[:, 4])


        ### save comparison data to file ###

        # If this is not the first step, we have previous data to compare
        # with:
        if prev_step_data is not None:
            sums = []
            # compare the waveguided bands:
            for j in [0, first_wg_band - 1, first_wg_band, first_wg_band + 1]:
                sums.append(
                    sum_of_squares(
                        prev_step_data[:, 5 + j],
                        data[:, 5 + j]
                    )
                )
            with open("sum_of_squares.dat", "a") as f:
                f.write('\t'.join([str(step)] + map(str, sums)) + '\n')
        # save data for next iteration:
        prev_step_data = data

        log.info(' ##### step={0} - success! #####\n\n'.format(step))

        # reset logger; the next stuff logged is going to next step's file:
        log.reset_logger()

    ### finally, save some plots ###

    data = np.loadtxt('sum_of_squares.dat')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # make double-logarithmic plot for every band:
    for i in range(data.shape[1] -1):
        ax.loglog(data[:,0], data[:,i + 1], 'o-',
                  label='band {0}'.format(i+1))

    compdata = np.power(data[:,0], -2)
    ax.loglog(data[:,0], compdata, '.:', label='$x^{-2}$')
    compdata = np.power(data[:,0], -4)
    ax.loglog(data[:,0], compdata, '.:', label='$x^{-4}$')
    compdata = np.power(data[:,0], -6)
    ax.loglog(data[:,0], compdata, '.:', label='$x^{-6}$')
    plt.legend()
    plt.title('sum of squared differences between simulation steps')
    fig.savefig('sum_of_squares.png')

if __name__ == '__main__':
    main()
