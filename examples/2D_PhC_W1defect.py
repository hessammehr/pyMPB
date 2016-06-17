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
sys.path.append('../')
from os import path, makedirs

import numpy as np
import matplotlib.pyplot as plt

from phc_simulations import TriHoles2D, TriHoles2D_yWaveguide
from kspace import KSpace
from utility import get_gap_bands
import log


def main():
    if len(sys.argv) > 1:
        mode = sys.argv[1]
    else:
        mode = 'sim'

    sim = TriHoles2D_yWaveguide(
            material='SiN',
            radius=0.375,
            numbands=10,
            k_steps=17,
            resolution=32,
            mesh_size=7,
            runmode=mode,
            num_processors=2,
            save_field_patterns=False,
            convert_field_patterns=False
    )

    if not sim:
        log.error('an error occured during simulation. See the .out file')
        return


def main2():
    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'

    resolution = 16
    material = 'SiN'
    radius = 0.375
    containing_folder = (
                    r'../yWaveguide_projections/'
                    r'2DtriH_{0}_res{1:03.0f}/'.format(material, resolution) +
                    r'radius0p{0:03.0f}'.format(radius * 1000))
    if not path.exists(path.abspath(containing_folder)):
        makedirs(path.abspath(containing_folder))

    # In the triangular lattice, in the basis of its reciprocal basis
    # vectors, this is the K' point, i.e. die boundary of the first
    # brillouin zone in the rectangular lattice, onto which we need to
    # project (see also : Steven G. Johnson et al., "Linear waveguides
    # in photonic-crystal slabs", Phys. Rev. B, Vol. 62, Nr.12,
    # 8212-8222 (2000); page 8216 & Fig. 8):
    rectBZ_K = np.array((0.25, -0.25))
    # the M point in the triangular lattice reciprocal basis, which points
    # along +X:
    # (note: if k_y is greater than 1/3, we leave the 1st BZ in +x
    # direction. But this is OK and we calculate it anyway, because it
    # does not change the projection. If we want to optimize
    # calculation time some time, we could limit this.)
    triBZ_M = np.array((0.5, 0.5))

    for k_y in np.linspace(0, 0.5, num=7):
        kspace = KSpace(
            points_list=[
                rectBZ_K * k_y * 2,
                rectBZ_K * k_y * 2 + triBZ_M
            ],
            k_interpolation=15,)

        sim = TriHoles2D(
                material=material,
                radius=radius,
                custom_k_space=kspace,
                numbands=10,
                resolution=resolution,
                mesh_size=7,
                runmode=mode,
                num_processors=2,
                containing_folder=containing_folder,
                save_field_patterns=False,
                convert_field_patterns=False,
                job_name_suffix='_ky{0:03.0f}'.format(k_y*1000),
                bands_title_appendix=', at ky={0:0.3f}'.format(k_y)
        )

        
        if not sim:
            log.error('an error occured during simulation. See the .out file')
            return
    
        # load te mode band data:
        fname = path.join(sim.workingdir, sim.jobname + '_tefreqs.csv')
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
