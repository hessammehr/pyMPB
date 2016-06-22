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
from phc_simulations import TriHoles2D_yWaveguide
import log
import defaults


def main():
    if len(sys.argv) > 1:
        mode = sys.argv[1]
    else:
        mode = 'sim'

    # monkey patching:
    # (don't output multiple tiles along x)
    defaults.mpbdata_call = (
        'mpb-data -T -rn%(resolution)s '
        '-y%(number_of_tiles_to_output)s '
        '-o%(output_file)s '
        '%(h5_file)s')
    defaults.number_of_tiles_to_output=5

    ksteps = 17

    sim = TriHoles2D_yWaveguide(
        material='SiN',
        radius=0.375,
        mode='te',
        numbands=24,
        k_steps=ksteps,
        supercell_x=5,
        resolution=32,
        mesh_size=7,
        runmode=mode,
        num_processors=2,
        projected_bands_folder='./projected_bands_repo',
        save_field_patterns_kvecs=[
            (0, x, 0) for x in np.linspace(0, 0.5, num=ksteps)],
        save_field_patterns_bandnums=[1, 2, 11, 12, 13, 22, 23],
        convert_field_patterns=True
    )

    if not sim:
        log.error('an error occurred during simulation. See the .out file')
        return


if __name__ == '__main__':
    main()
