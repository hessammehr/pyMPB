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


from simulation import Simulation
from geometry import Geometry
from kspace import KSpaceTriangular, KSpace
from objects import Dielectric, Rod, Block
import defaults
import log
from utility import do_runmode
from os import path, makedirs
import numpy as np


def TriHoles2D(
        material, radius, numbands=8, k_interpolation=11, 
        resolution=32, mesh_size=7,
        runmode='sim', num_processors=2,
        save_field_patterns=True, convert_field_patterns=True,
        containing_folder='./',
        job_name_suffix='', bands_title_appendix='',
        custom_k_space=None, modes=('te', 'tm')):
    """Create a 2D MPB Simulation of a triangular lattice of holes.

    :param material: can be a string (e.g. SiN,
    4H-SiC-anisotropic_c_in_z; defined in data.py) or just the epsilon
    value (float)
    :param radius: the radius of holes in units of the lattice constant
    :param numbands: number of bands to calculate
    :param k_interpolation: number of the k-vectors between every two of
    the used high symmetry points Gamma, M, K and Gamma again, so the
    total number of simulated k-vectors will be 3*k_interpolation + 4.
    Only used if no custom_custom_k_space is provided.
    :param resolution: described in MPB documentation
    :param mesh_size: described in MPB documentation
    :param runmode: can be one of the following:
        ''       : just create and return the simulation object
        'ctl'    : create the sim object and save the ctl file
        'sim' (default): run the simulation and do all postprocessing
        'postpc' : do all postprocessing; simulation should have run
                   before!
        'display': display all pngs done during postprocessing. This is
                   the only mode that is interactive.
    :param num_processors: number of processors used during simulation
    :param save_field_patterns: indicates whether field pattern h5 files
    are generated during the simulation (at points of high symmetry)
    :param convert_field_patterns: indicates whether field pattern h5
    files should be converted to png (only when postprocessing)
    :param containing_folder: the path to the folder which will contain
    the simulation subfolder.
    :param job_name_suffix: Optionally specify a job_name_suffix
    (appendix to the folder name etc.) which will be appended to the
    jobname created automatically from the most important parameters.
    :param bands_title_appendix: will be added to the title of the bands
    diagram.
    :param custom_k_space: By default, KSpaceTriangular with
    k_interpolation interpolation steps are used. Provide any KSpace
    object here to customize this. k_interpolation will then be ignored.
    :param modes: a list of modes to run. Possible are 'te' and 'tm'.
    Default: both
    :return: the Simulation object

    """
    mat = Dielectric(material)    

    geom = Geometry(
        width=1,
        height=1,
        triangular=True,
        objects=[
            Rod(
                x=0,
                y=0,
                material='air',
                radius=radius)])

    if isinstance(custom_k_space, KSpace):
        kspace = custom_k_space
    else:
        kspace = KSpaceTriangular(k_interpolation=k_interpolation)

    # points of interest: (output mode patterns at these points)
    if save_field_patterns:
        poi = kspace.points()[0:-1]
    else:
        poi = []

    runcode = ''
    for mode in modes:
        if mode == 'te':
            outputfunc = ' '.join(defaults.output_funcs_te)
        else:
            outputfunc = ' '.join(defaults.output_funcs_tm)
        runcode += (
            '(run-%s %s)\n' % (
                mode, defaults.default_band_func(poi, outputfunc)
            ) +
            '(print-dos 0 1.2 121)\n\n')

    jobname = 'TriHoles2D_{0}_r{1:03.0f}'.format(
                    mat.name, radius * 1000)

    sim = Simulation(
        jobname=jobname + job_name_suffix,
        geometry=geom,
        kspace=kspace,
        numbands=numbands,
        resolution=resolution,
        mesh_size=mesh_size,
        initcode=defaults.default_initcode +
                 '(set! default-material {0})'.format(str(mat)),
        postcode='',
        runcode=runcode,
        work_in_subfolder=path.join(
            containing_folder, jobname + job_name_suffix),
        clear_subfolder=runmode.startswith('s') or runmode.startswith('c'))

    draw_bands_title = ('2D hex. PhC; {0}, radius={1:0.3f}'.format(
                            mat.name, geom.objects[0].radius) + 
                        bands_title_appendix)

    return do_runmode(
        sim, runmode, num_processors, draw_bands_title,
        plot_crop_y=True, # automatic cropping
        convert_field_patterns=convert_field_patterns,
        # don't add gamma point a second time (index 3):
        field_pattern_plot_k_slice=(0,2),
        x_axis_hint=[defaults.default_x_axis_hint, kspace][kspace.has_labels()]
    )


def TriHolesSlab3D(
        material, radius, thickness, numbands=8, k_interpolation=11, 
        resolution=32, mesh_size=7, supercell_z=6,
        runmode='sim', num_processors=2,
        save_field_patterns=True, convert_field_patterns=True,
        job_name_suffix='', bands_title_appendix='',
        modes=('zeven', 'zodd')):
    """Create a 3D MPB Simulation of a slab with a triangular lattice of
    holes.

    :param material: can be a string (e.g. SiN,
    4H-SiC-anisotropic_c_in_z; defined in data.py) or just the epsilon
    value (float)
    :param radius: the radius of holes in units of the lattice constant
    :param thickness: slab thickness in units of the lattice constant
    :param numbands: number of bands to calculate
    :param k_interpolation: number of the k-vectors between every two of
    the used high symmetry points Gamma, M, K and Gamma again, so the
    total number of simulated k-vectors will be 3*k_interpolation + 4
    :param resolution: described in MPB documentation
    :param mesh_size: described in MPB documentation
    :param supercell_z: the height of the supercell in units of the
    lattice constant
    :param runmode: can be one of the following:
        ''       : just create and return the simulation object
        'ctl'    : create the sim object and save the ctl file
        'sim' (default): run the simulation and do all postprocessing
        'postpc' : do all postprocessing; simulation should have run
                   before!
        'display': display all pngs done during postprocessing. This is
                   the only mode that is interactive.
    :param num_processors: number of processors used during simulation
    :param save_field_patterns: indicates whether field pattern h5 files
    are generated during the simulation (at points of high symmetry)
    :param convert_field_patterns: indicates whether field pattern h5
    files should be converted to png (only when postprocessing)
    :param job_name_suffix: Optionally specify a job_name_suffix
    (appendix to the folder name etc.) which will be appended to the
    jobname created automatically from the most important parameters.
    :param bands_title_appendix: will be added to the title of the bands
    diagram.
    :param modes: a list of modes to run. Possible are 'zeven' and
    'zodd'. Default: both
    :return: the Simulation object

    """
    mat = Dielectric(material)

    geom = Geometry(
        width=1,
        height=1,
        depth=supercell_z,
        triangular=True,
        objects=[
            Block(
                x=0, y=0, z=0,
                material=mat,
                #make it bigger than computational cell, just in case:
                size=(2, 2, thickness)), 
            Rod(
                x=0,
                y=0, 
                material='air',
                radius=radius)])

    kspace = KSpaceTriangular(k_interpolation=k_interpolation)

    # points of interest: (output mode patterns at these points)
    if save_field_patterns:
        poi = kspace.points()[0:-1]
    else:
        poi = []

    runcode = ''
    for mode in modes:
        if mode == 'zeven':
            outputfunc = ' '.join(defaults.output_funcs_te)
        else:
            outputfunc = ' '.join(defaults.output_funcs_tm)
        runcode += (
            '(run-%s %s)\n' % (
                mode, defaults.default_band_func(poi, outputfunc)
            ) +
            '(print-dos 0 1.2 121)\n\n')

    jobname = 'TriHolesSlab_{0}_r{1:03.0f}_t{2:03.0f}'.format(
                    mat.name, radius * 1000, thickness * 1000)

    sim = Simulation(
        jobname=jobname + job_name_suffix,
        geometry=geom,
        kspace=kspace,
        numbands=numbands,
        resolution=resolution,
        mesh_size=mesh_size,
        initcode=defaults.default_initcode,
        postcode='',
        runcode=runcode,
        clear_subfolder=runmode.startswith('s') or runmode.startswith('c'))

    draw_bands_title = ('Hex. PhC slab; '
                        '{0}, thickness={1:0.3f}, radius={2:0.3f}'.format(
                            mat.name, 
                            geom.objects[0].size[2], 
                            geom.objects[1].radius) +
                        bands_title_appendix)
    return do_runmode(
        sim, runmode, num_processors, draw_bands_title,
        plot_crop_y=0.8, # crop at 0.8
        convert_field_patterns=convert_field_patterns,
        # don't add gamma point a second time (index 3):
        field_pattern_plot_k_slice=(0,2),
        x_axis_hint=kspace
    )


def TriHoles2D_yWaveguide(
        material, radius, mode='te', numbands=8, k_steps=17,
        supercell_x=5, resolution=32, mesh_size=7,
        runmode='sim', num_processors=2,
        projected_bands_folder='../projected_bands_repo',
        save_field_patterns_kvecs=list(), save_field_patterns_bandnums=list(),
        convert_field_patterns=False,
        job_name_suffix='', bands_title_appendix=''):
    """Create a 2D MPB Simulation of a triangular lattice of holes, with
    a waveguide along the nearest neighbor direction, i.e. Gamma->K
    direction.

    The simulation is done with a rectangular super cell.

    Before the waveguide simulation, additional simulations of the
    unperturbed structure will be run for projected bands data, if these
    simulations where not run before.

    :param material: can be a string (e.g. SiN,
    4H-SiC-anisotropic_c_in_z; defined in data.py) or just the epsilon
    value (float)
    :param radius: the radius of holes in units of the lattice constant
    :param mode: the mode to run. Possible are 'te' and 'tm'.
    :param numbands: number of bands to calculate
    :param k_steps: number of k_y steps between 0 and 0.5 to simulate
    :param supercell_x: the length of the supercell perpendicular to the
    waveguide, in units of sqrt(3) times the lattice constant. It it is
    not a odd number, one will be added.
    :param resolution: described in MPB documentation
    :param mesh_size: described in MPB documentation
    :param runmode: can be one of the following:
        ''       : just create and return the simulation object
        'ctl'    : create the sim object and save the ctl file
        'sim' (default): run the simulation and do all postprocessing
        'postpc' : do all postprocessing; simulation should have run
                   before!
        'display': display all pngs done during postprocessing. This is
                   the only mode that is interactive.
    :param num_processors: number of processors used during simulation
    :param projected_bands_folder: the path to the folder which will
    contain the simulations of the unperturbed PhC, which is needed for
    the projections along k_x. If the folder contains simulations run
    before, their data will be reused.
    :param save_field_patterns_kvecs: a list of k-vectors (3-tuples),
    which indicates where field pattern h5 files are generated during
    the simulation (only at bands in save_field_patterns_bandnums)
    :param save_field_patterns_bandnums: a list of band numbers (int,
    starting at 1), which indicates where field pattern h5 files are
    generated during the simulation (only at k-vectors in
    save_field_patterns_kvecs)
    :param convert_field_patterns: indicates whether field pattern h5
    files should be converted to png (only when postprocessing)
    :param job_name_suffix: Optionally specify a job_name_suffix
    (appendix to the folder name etc.) which will be appended to the
    jobname created automatically from the most important parameters.
    :param bands_title_appendix: will be added to the title of the bands
    diagram.
    :return: the Simulation object

    """
    mat = Dielectric(material)

    # first, make sure all data for projected bands exist, otherwise
    # start their simulations.

    unperturbed_jobname = 'TriHoles2D_{0}_r{1:03.0f}'.format(
        mat.name, radius * 1000)
    # look here for old simulations, and place new ones there:
    repo = path.abspath(
        path.join(
            path.curdir,
            projected_bands_folder,
            unperturbed_jobname
        )
    )
    # create path if not there yet:
    if not path.exists(path.abspath(repo)):
        makedirs(path.abspath(repo))
    # these ky points will be simulated::
    ky_points = np.linspace(0, 0.5, num=k_steps, endpoint=True)

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

    # This list will be forwarded later to this defect simulation's
    # post-process. It contains the folder paths of unperturbed
    # simulations for each k-vec of this simulation:
    project_bands_list = []

    # now, see if we need to simulate:
    for ky in ky_points:
        jobname_suffix = '_ky{0:06.0f}'.format(ky*1e6)
        jobname = unperturbed_jobname + jobname_suffix
        project_bands_list.append(path.join(repo, jobname))
        range_file_name = path.join(
            repo, jobname, jobname + '_' + mode + '_ranges.csv')
        if not path.isfile(range_file_name):
            # does not exist, so start simulation:
            log.info('unperturbed structure not yet simulated at '
                     'k_y={0}. Running now...'.format(ky))
            kspace = KSpace(
                points_list=[
                    rectBZ_K * ky * 2,
                    rectBZ_K * ky * 2 + triBZ_M
                ],
                k_interpolation=15,)

            sim = TriHoles2D(
                    material=material,
                    radius=radius,
                    custom_k_space=kspace,
                    numbands=defaults.num_projected_bands,
                    resolution=resolution,
                    mesh_size=mesh_size,
                    runmode='sim',
                    num_processors=num_processors,
                    containing_folder=repo,
                    save_field_patterns=False,
                    convert_field_patterns=False,
                    job_name_suffix=jobname_suffix,
                    bands_title_appendix=', at ky={0:0.3f}'.format(ky),
                    modes=[mode]
            )

            if not sim:
                log.error(
                    'an error occurred during simulation of unperturbed '
                    'structure. See the .out file in {0}'.format(path.join(
                        repo, jobname
                    ))
                )
                return

    # make it odd:
    if supercell_x % 2 == 0:
        supercell_x += 1
    # half of the supercell (floored):
    scxh = int(supercell_x / 2)

    # Create geometry and add objects.
    # Note: (0, 0, 0) is the center of the unit cell.
    geom = Geometry(
        width='(* (sqrt 3) %i)' % supercell_x,
        height=1,
        triangular=False,
        objects=([
            # center holes:
            Rod(
                x='(* %i (sqrt 3))' % cx,
                y=0,
                material='air',
                radius=radius)
            for cx in list(range(-scxh, 0)) + list(range(1, scxh + 1))] +

            # perimeter holes:
            [Rod(
                x='(* {0:.1f} (sqrt 3))'.format(cx + 0.5),
                y=0.5,
                material='air',
                radius=radius)
            for cx in range(-scxh, scxh + 1)]
        )
    )

    kspace = KSpace(
        points_list=[(0, 0, 0), (0, 0.5, 0)],
        k_interpolation=k_steps - 2,
    )

    jobname = 'TriHoles2D_W1_{0}_r{1:03.0f}'.format(
                    mat.name, radius * 1000)

    if mode == 'te':
        outputfuncs = defaults.output_funcs_te
    else:
        outputfuncs = defaults.output_funcs_tm

    if save_field_patterns_bandnums and save_field_patterns_kvecs:
        runcode = (
            ';function to determine whether an item x is member of list:\n'
            '(define (member? x list)\n'
            '    (cond (\n'
            '        ;false if the list is empty:\n'
            '        (null? list) #f )\n'
            '        ;true if first item (car) equals x:\n'
            '        ( (eqv? x (car list)) #t )\n'
            '        ;else, drop first item (cdr) and make recursive call:\n'
            '        ( else (member? x (cdr list)) )\n'
            '    ))\n\n' +
            '(define output-bands-list (list {0}))\n\n'.format(' '.join(
                map(str, save_field_patterns_bandnums))) +
            '(define (output-func bnum)\n'
            '    (if (member? bnum output-bands-list)\n'
            '        (begin\n' +
            ''.join(12 * ' ' + '({0} bnum)\n'.format(func)
                    for func in outputfuncs) +
            '        )\n'
            '    ))\n\n'
            '(run-{0} {1})\n'.format(
                mode,
                defaults.default_band_func(
                    save_field_patterns_kvecs, 'output-func')) +
            '(print-dos 0 1.2 121)\n\n'
        )
    else:
        runcode = ('(run-{0} {1})\n'.format(
                mode,
                defaults.default_band_func([], None)
            ) +
            '(print-dos 0 1.2 121)\n\n')

    sim = Simulation(
        jobname=jobname + job_name_suffix,
        geometry=geom,
        kspace=kspace,
        numbands=numbands,
        resolution=resolution,
        mesh_size=mesh_size,
        initcode=defaults.default_initcode +
                 '(set! default-material {0})'.format(str(mat)),
        postcode='',
        runcode=runcode,
        clear_subfolder=runmode.startswith('s') or runmode.startswith('c'))

    draw_bands_title = (
        '2D hex. PhC W1; {0}, radius={1:0.3f},'.format(
            mat.name, geom.objects[0].radius) +
        bands_title_appendix)

    return do_runmode(
        sim, runmode, num_processors, draw_bands_title,
        plot_crop_y=False, # no cropping
        convert_field_patterns=convert_field_patterns,
        # don't add gamma point a second time (index 3):
        field_pattern_plot_k_slice=(0,2),
        x_axis_hint=5,
        project_bands_list=project_bands_list
    )

