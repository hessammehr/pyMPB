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

from phc_simulations import TriHoles2D_yWaveguide
import log


def main():
    if len(sys.argv) > 1:
        mode = sys.argv[1]
    else:
        mode = 'sim'

    sim = TriHoles2D_yWaveguide(
            material='SiN',
            radius=0.375,
            mode='te',
            numbands=24,
            k_steps=17,
            resolution=32,
            mesh_size=7,
            runmode=mode,
            num_processors=2,
            projected_bands_folder='./projected_bands_repo',
            save_field_patterns=False,
            convert_field_patterns=False
    )

    if not sim:
        log.error('an error occurred during simulation. See the .out file')
        return


if __name__ == '__main__':
    main()
