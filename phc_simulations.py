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

from simulation import Simulation
from geometry import Geometry
from kspace import KSpaceTriangular
from objects import Dielectric, Rod, Block
import defaults
from utility import do_runmode

def TriHoles2D(
        material, radius, numbands=8, k_interpolation=11, 
        resolution=32, mesh_size=7, runmode='sim',
        num_processors=2, convert_field_patterns=True, 
        job_name_suffix='', bands_title_appendix=''):
    """Create a 2D MPB Simulation of a triangular lattice of holes.
    material can be a string (e.g. SiN, 4H-SiC-anisotropic_c_in_z;
    defined in data.py) or just the epsilon value (float),
    radius is the radius of holes in units of the lattice constant,
    numbands is number of bands to calculate.
    k_interpolation is number of the k-vectors between every two of
    the used high symmetry points Gamma, M, K and Gamma again, so the
    total number of simulated k-vectors will be 3*k_interpolation + 4.
    resolution and mesh_size are as described in MPB documentation.
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

    kspace = KSpaceTriangular(k_interpolation=k_interpolation)

    # points of interest: (output mode patterns at these points)
    poi = kspace.points()[0:-1]

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
        runcode='(run-tm %s)\n' % defaults.default_band_func_tm(poi) +
                '(print-dos 0 1.2 121)\n\n' +
                '(run-te %s)\n' % defaults.default_band_func_te(poi) +
                '(print-dos 0 1.2 121)\n\n',
        clear_subfolder=runmode.startswith('s') or runmode.startswith('c'))

    draw_bands_title = ('2D hex. PhC; {0}, radius={1:0.3f}'.format(
                            mat.name, geom.objects[0].radius) + 
                        bands_title_appendix)

    return do_runmode(sim, runmode, num_processors, draw_bands_title,
                      plot_crop_y=True, # automatic cropping
                      convert_field_patterns=convert_field_patterns,
                      # don't add gamma point a second time (index 3):
                      field_pattern_plot_k_slice=(0,2),
                      x_axis_hint=kspace
                     )


def TriHolesSlab3D(
        material, radius, thickness, numbands=8, k_interpolation=11, 
        resolution=32, mesh_size=7, supercell_z=6, runmode='sim',
        num_processors=2, convert_field_patterns=True, 
        job_name_suffix='', bands_title_appendix=''):
    """Create a 3D MPB Simulation of a slab with a triangular lattice of 
    holes. material can be a string (e.g. SiN, 4H-SiC-anisotropic_c_in_z;
    defined in data.py) or just the epsilon value (float),
    radius is the radius of holes, thickness is the slab thickness, both in
    units of the lattice constant, numbands is number of bands to calculate.
    k_interpolation is number of the k-vectors between every two of
    the used high symmetry points Gamma, M, K and Gamma again, so the
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
    poi = kspace.points()[0:-1]

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
        runcode='(run-zodd %s)\n' % defaults.default_band_func_tm(poi) +
                '(print-dos 0 1.2 121)\n\n' +
                '(run-zeven %s)\n' % defaults.default_band_func_te(poi) +
                '(print-dos 0 1.2 121)\n\n',
        clear_subfolder=runmode.startswith('s') or runmode.startswith('c'))

    draw_bands_title = ('Hex. PhC slab; '
                        '{0}, thickness={1:0.3f}, radius={2:0.3f}'.format(
                            mat.name, 
                            geom.objects[0].size[2], 
                            geom.objects[1].radius) +
                        bands_title_appendix)

    return do_runmode(sim, runmode, num_processors, draw_bands_title,
                      plot_crop_y=0.8, # crop at 0.8
                      convert_field_patterns=convert_field_patterns,
                      # don't add gamma point a second time (index 3):
                      field_pattern_plot_k_slice=(0,2),
                      x_axis_hint=kspace
                     )