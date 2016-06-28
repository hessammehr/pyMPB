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

from phc_simulations import TriHolesSlab3D
import data
import log
import defaults


# no inversion symmetry!
defaults.mpb_call = defaults.mpb_call.replace('mpbi', 'mpb')

# add Silicon and glass:
data.refr_index['Si'] = 3.48627606
data.refr_index['glass'] = 1.5078
data.update_material_names()
data.update_dielectrics()

# pitch in nm:
pitch = 600.0


def main():
    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'

    sim = TriHolesSlab3D(
        material='Si',
        substrate_material='glass',
        radius=0.5 * 367 / pitch,
        thickness=116 / pitch,
        numbands=16,
        k_interpolation=31,
        resolution=32,
        mesh_size=7,
        supercell_z=10,
        runmode=mode,
        num_processors=8,
        save_field_patterns=False,
        convert_field_patterns=False,
        job_name_suffix='_on_glass',
        bands_title_appendix=' (on glass)',
        modes=[''] # use (run), not (run-zeven) etc.
    )
        
    if not sim:
        log.error('an error occurred during simulation. See the .out file')
    else:
        log.info(' ##### Si PhC slab on glass - success! #####\n\n')

if __name__ == '__main__':
    main()
