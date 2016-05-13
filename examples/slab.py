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
from simulation import Simulation
from geometry import Geometry
from kspace import KSpace
from objects import Dielectric, Block
from utility import do_runmode
import defaults
import log

def UniformSlab3D(
        material, numbands=8, k_interpolation=11,
        resolution=32, mesh_size=3, supercell_z=6, runmode='sim',
        num_processors=2, convert_field_patterns=True,
        job_name_suffix='', bands_title_appendix=''):
    """Create a 3D MPB Simulation of an unperturbed slab with thickness 1.

    material can be a string (e.g. SiN, 4H-SiC-anisotropic_c_in_z;
    defined in data.py) or just the epsilon value (float),
    numbands is the number of bands to calculate.
    k_interpolation is the number of k-vectors between every two of
    the used high symmetry points Gamma, X, M and Gamma again, so the
    total number of simulated k-vectors will be 3*k_interpolation + 4.
    resolution and mesh_size are as described in MPB documentation.
    supercell_z is the height of the supercell in units of lattice constant;
    runmode can be one of the following:
        ''       : just create and return the simulation object
        'ctl'    : create the sim object and save the ctl file
        'sim' (default): run the simulation and do all postprocessing
        'postpc' : do all postprocessing; simulation should have run before!
        'display': display all pngs done during postprocessing. This is the
                   only mode that is interactive.
    convert_field_patterns indicates whether field pattern h5 files
    should be converted to png (only when postprocessing).
    Optionally specify a job_name_suffix (appendix to the folder name etc.)
    which will be appended to the jobname created automatically from
    the most important parameters.
    Specify bands_title_appendix, which will be added to the title
    of the bands diagram.

    """
    #global geom # for debugging

    mat = Dielectric(material)
    # place the slab sideways, i.e the surface normal in y direction
    geom = Geometry(
        # I want to hide that there is actually a periodicity in x. The
        # bands are always cut of at kvec_x=0.5 and folded back. If I
        # decrease the size in x (e.g. by 1/8), the kvec_x at 0.5 will
        # actually be greater (here, 0.5 * 8). This is true, as can also
        # be seen by looking at k_magnitude in the simulation output. To
        # make sure there is no mixing between bands, I will only look
        # at kvec_x up to 0.25. At kvec=(1/8, 0, 0), the band
        # frequencies are the same as it would be at k_x*d/pi=1 (d is
        # slab thickness), when we had no periodic boundaries (compare
        # Sakoda - Optical Properties of Photonic Crystals 2nd Edition
        # 2005, page 177, Fig. 8.2 - actually, the frequencies are not
        # _exactly_ the same, why? Wrong refractive index?):
        width=0.125,
        height=1,
        depth=supercell_z,
        triangular=False,
        objects=[
            Block(
                x=0, y=0, z=0,
                material=mat,
                #make it bigger than computational cell, just in case:
                size=(1, 2, 1))])
    # Just below kx=0.5, the bands still mix due to the periodicity.
    # To hide this, just cut at kx=0.25:
    max_kx = 0.25
    # I actually want to see TM-like and TE-like modes, like described
    # in Joannopoulos et al., Photonic Crystals; Molding the Flow of Light,
    # Princeton University Press 2008:
    # Fig. 7 on page 130 & text on pages 128 and 136.
    # I try to look everywhere in kspace to find them, but it seems there
    # are only pure TE and TM modes.
    # (These ...-like modes are neither mentioned in
    # Sakoda - Optical Properties of Photonic Crystals 2nd Edition,
    # Springer 2005,pages 176ff.)
    kspace = KSpace(
        points_list=[
            (0, 0, 0), (max_kx, 0, 0),  # Gamma -> +x
            (max_kx, 0, 0.5),           # +x -> shifted +x (no bands change)
            (0, 0, 0.5), (0, 0, 0),     # -> +z -> Gamma
            (max_kx, 0, 0.5)       # Gamma -> shifted +x (same as Gamma->+x)
        ],
        k_interpolation=k_interpolation)

    # points of interest: (output mode patterns at these points)
    poi = [
        (0, 0, 0), (0.125, 0, 0), 
        (0.25, 0, 0), (0.25, 0, 0.25), 
        (0.25, 0, 0.5), (0.125, 0, 0.5),
        (0, 0, 0.5), (0, 0, 0.25),
        (0.125, 0, 0.25) 
    ]
    band_funcs = (
        '\n    display-zparities display-yparities' +
        ''.join([(
            '\n    (output-at-kpoint (vector3 %s) '
            'fix-efield-phase output-efield output-hfield)'
        ) % ' '.join(str(c) for c in vec) for vec in poi ])
    )

    jobname = 'Slab3D_{0}_res{1}_supcell{2}'.format(
        mat.name, resolution, supercell_z).replace('.', 'p')

    sim = Simulation(
        jobname=jobname + job_name_suffix,
        geometry=geom,
        kspace=kspace,
        numbands=numbands,
        resolution=resolution,
        mesh_size=mesh_size,
        initcode=defaults.default_initcode,
        postcode='',
        runcode='(run %s)\n' % band_funcs,
        clear_subfolder=runmode.startswith('s') or runmode.startswith('c'))

    draw_bands_title = ('Uniform slab; '
                        '{0}, resolution={1}, supercell_z={2}'.format(
                            mat.name,
                            resolution,
                            supercell_z) +
                        bands_title_appendix)

    ### some monkey patching: >:O  <-evil monkey
    # only make cross-section png:
    defaults.epsh5topng_call_3D = defaults.epsh5topng_call_3D_cross_sect
    # export field patterns along another slice (default: -0z0)
    # NOTE, for mpbdata_call: 
    # -T because we use mpi-mpb, then -y is actually for -x
    defaults.mpbdata_call = ('mpb-data -T -rn%(resolution)s '
                '-y8 '
                '-o%(output_file)s '
                '%(h5_file)s')
    # For the individual field components to be comparable, don't use
    # automatic scaling, but set range with -m and -M:
    defaults.fieldh5topng_call_3D = (
        'h5topng -0y0 -S3 -Zcbluered -C%(eps_file)s -m-1 -M1 '
        '-o%(output_file)s %(h5_file)s')
    defaults.fieldh5topng_call_3D_no_ovl = (
        'h5topng -0y0 -S3 -Zcbluered -m-1 -M1 '
        '-o%(output_file_no_ovl)s %(h5_file)s')

    return do_runmode(
        sim, runmode, num_processors, draw_bands_title,
        plot_crop_y=1,
        convert_field_patterns=convert_field_patterns,
        field_pattern_plot_k_slice=None,
        x_axis_hint=11)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'
    sim = UniformSlab3D(
        material=2.86**2,
        numbands=10,
        k_interpolation=31,
        resolution=64,
        supercell_z=8,
        runmode=mode,
        num_processors=8)
    if sim:
        log.info(' ##### uniform slab - success! #####\n\n')
