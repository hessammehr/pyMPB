    #Copyright 2016 Juergen Probst
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

import sys
sys.path.append('../')
from os import path

import numpy as np
import matplotlib.pyplot as plt

from phc_simulations import TriHoles2D
from utility import get_gap_bands
import log

def main():
    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'

    minrad = 0.2
    maxrad = 0.4
    radstep = 0.05
    numsteps = int((maxrad - minrad) / radstep + 1.5)
    steps = np.linspace(minrad, maxrad, num=numsteps, endpoint=True)
   
    for i, radius in enumerate(steps):
        
        log.info("running simulation with {0:n} radius steps:\n{1}".format(
            numsteps, steps) +
            '\n  ### current step: #{0} ({1}) ###\n'.format(i+1, radius))

        sim = TriHoles2D(
            material='SiN',
            radius=radius,
            numbands=4,#8,
            k_interpolation=5,#31, 
            resolution=16, 
            mesh_size=7,
            runmode=mode,
            num_processors=2,
            convert_field_patterns=True)
        
        if not sim:
            log.error('an error occured during simulation. See the .out file')
            return
    
        # load te mode band data:
        fname = path.join(sim.workingdir, sim.jobname + '_te.csv')
        data = np.loadtxt(fname, delimiter=',', skiprows=1)
        gapbands = get_gap_bands(data[:, 5:])
        
        # maybe there is no gap?
        if len(gapbands) == 0:
            gap = 0
        elif gapbands[0][0] != 1:
            # it must be a gap between band 1 and 2
            gap = 0
        else:
            gap = gapbands[0][3]

        # save gap sizes to file (first TE gap):
        with open("gaps.dat", "a") as f:
            f.write("{0}\t{1}\n".format(radius, gap))
    
        log.info(' ##### radius={0} - success! #####\n\n'.format(radius))
        
        # reset logger; the next stuff logged is going to next step's file:
        log.reset_logger()

    data = np.loadtxt('gaps.dat')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data[:,0], data[:,1], 'o-')
    fig.savefig('gaps.png')

if __name__ == '__main__':
    main()
