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
import numpy as np
from phc_simulations import TriHoles2D_Waveguide
import log
import defaults


def main():
    if len(sys.argv) > 1:
        mode = sys.argv[1]
    else:
        mode = 'sim'

    # monkey patching:
    # (don't output multiple tiles along y)
    # (x and y are switched because of -T)
    defaults.mpbdata_call = (
        'mpb-data -T -rn%(resolution)s '
        '-y%(number_of_tiles_to_output)s '
        '-o%(output_file)s '
        '%(h5_file)s')
    defaults.number_of_tiles_to_output=5

    ksteps = 17

    # width of the supercell:
    supercell_size = 9

    # The bigger the unit cell, the more bands fold in below the
    # waveguide band:
    first_wg_band = int(3 + 2 * (supercell_size - 1))

    sim = TriHoles2D_Waveguide(
        material='SiN',
        radius=0.380,
        mode='te',
        numbands=40,
        k_steps=ksteps,
        supercell_size=supercell_size,
        resolution=32,
        mesh_size=7,
        runmode=mode,
        num_processors=2,
        projected_bands_folder='./projected_bands_repo',
        save_field_patterns_kvecs=[
            (x, 0, 0) for x in np.linspace(0, 0.5, num=ksteps)],
        save_field_patterns_bandnums=[
            1, 2,
            first_wg_band, first_wg_band + 1, first_wg_band + 2],
        convert_field_patterns=True,
        field_pattern_plot_k_selection=[0, 5, 7, 9, 11, 13, 15, 16]
    )

    if not sim:
        log.error('an error occurred during simulation. See the .out file')
        return


if __name__ == '__main__':
    main()
