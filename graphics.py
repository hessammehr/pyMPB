    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
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
import matplotlib.pyplot as plt
#from pylab import figure,show,linspace,savefig,text,griddata,plot,contour,clf,clabel,colorbar
from matplotlib.mlab import griddata
from matplotlib.patches import Ellipse
from matplotlib.text import Text
import numpy as np
from numpy import loadtxt
from utility import max_epsilon, get_gap_bands
from defaults import fig_size, contour_lines, contour_filled, contour_plain,\
    colorbar_style, default_x_axis_hint
from kspace import KSpace
from bandplotter import BandPlotter
import axis_formatter
import objects
import log

def draw_geometry(
        geometry,jobname,format='pdf', display=True, block_when_showing=True,
        anisotropic_component=0):
    """#FIXME: needs to be adapted to work with anisotropic material
    
    Draw and save the geometry. Save as file jobname+.+format. 
    Only show figure if display=True (default). When showing, block script if 
    block_when_showing=True (default). In case of anisotropic material, only 
    draw anisotropic_component (default 0).
    """
    global maxeps
    # I commented clf(), because it just opens an unnecessary new empty
    # figure window, the same that figure() does below.
    #plt.clf()
    maxeps = max_epsilon(geometry, anisotropic_component)
    drawing_dict = {objects.Rod : draw_rod}
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111,aspect='equal')
    X = geometry.width/2
    Y = geometry.height/2
    ax.set_xlim(-X,X)
    ax.set_ylim(-Y,Y)
    
    graphic_objs = (
        [drawing_dict[obj.__class__](index, obj, anisotropic_component) 
            for index,obj in enumerate(geometry.objects)])
    for graph in graphic_objs:
        for elem in graph:
            ax.add_artist(elem)
    
    plt.savefig(jobname+'.'+format,format=format,transparent=True)
    if display:
        plt.show(block=block_when_showing)
    else:
        plt.close(fig)

def draw_rod(index, rod, anisotropic_component=0):
    """#FIXME needs to be adapted to work with anisotropic material
    
    """
    if isinstance(rod.material.epsilon, (list, tuple)):
        eps = rod.material.epsilon[anisotropic_component]
        epsstr = ', '.join(['{0:.2f}'.format(f) for f in rod.material.epsilon])
    else:
        eps = rod.material.epsilon
        epsstr = '{0:.2f}'.format(eps)
        
    return (
        Ellipse(
            (rod.x,rod.y),
            rod.radius*2,
            rod.radius*2,
            alpha=eps / maxeps,
            lw=2.0,
            ls='dashed',
            ec='darkblue',
            fc='lightblue'),
        Text(
            rod.x,
            rod.y,
            '#%s\n$\epsilon=%s$'% (index,epsstr),
            ha='center',
            va='center',
            family='sans-serif'))

def draw_bandstructure_2D(
        jobname, mode, kspace, band, ext='.csv', format='pdf', filled=True,
        levels=15, lines=False, labeled=False, legend=False):
    """Draw 2D band contour map of one band."""
    #clf()
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111, aspect='equal')
    x,y,z = loadtxt(
        "{0}_{1}{2}".format(jobname, mode, ext),
        delimiter=', ',
        skiprows=1,
        usecols=(1, 2, 4 + band),
        unpack=True)
    if hasattr(kspace, 'x_steps') and hasattr(kspace, 'y_steps'):
        # KSpace was created by KSpaceRectangularGrid
        xi = np.linspace(-0.5, 0.5, kspace.x_steps)
        yi = np.linspace(-0.5, 0.5, kspace.y_steps)
        zi = griddata(x, y, z, xi, yi, interp='linear')
        if filled:
            cs = ax.contourf(xi, yi, zi, levels, **contour_filled)
            legend and plt.colorbar(cs, **colorbar_style)
            cs = lines and ax.contour(xi, yi, zi, levels, **contour_lines)
            labeled and lines and plt.clabel(cs, fontsize=8, inline=1)
        else:
            cs = ax.contour(xi, yi, zi, levels, **contour_plain)
            legend and plt.colorbar(cs, **colorbar_style)
            labeled and plt.clabel(cs, fontsize=8, inline=1)
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 0.5)
    else:
        plt.plot(x, y, z)
    plt.savefig(jobname,format=format,transparent=True)
    plt.show()

def draw_bands(
        jobname, modes, x_axis_hint=default_x_axis_hint,
        custom_plotter=None, title='', crop_y=True, light_cone=False):
    """Plot dispersion relation of all bands calculated along all k
    vectors.

    The band data is loaded from previously saved .csv files.
    (filenames: [*jobname* + '_' + *mode* + 'freqs.csv' for mode in
    modes])

    *x_axis_hint* gives a hint on which kind of ticks and labels should
    be shown on the x-axis and provides the data needed. *x_axis_hint*
    can be one of the following:
        -- integer number:
            The axis' labels will be the 3D k-vectors. The number
            denotes the number of major ticks and labels distributed on
            the axis.
        -- [integer, format-string]:
            Same as above, but the labels are formatted with the
            format-string - this gives the possibility to only show one
            of the three vector components, e.g. the string "{2}" to
            only show the k-vector's z-component. The axis title will be
            inferred from the format-string.
        -- KSpace object:
            This must be a KSpace object created with point_labels.
            These labels usually denote the high symmetry or crititical
            points, and they will be shown on the axis.
        -- CustomAxisFormatter object:
            This gives the possibility to completely customize the
            x-axis' tick positions, tick labels and axis label. If the
            CustomAxisFormatter's hover data have not been set, it will
            be set here with the k-vectors read from the .csv file.

    If you want to add the graph to an existing figure, supply a
    BandPlotter with *custom_plotter*, otherwise (default:
    custom_plotter=None) a new BandPlotter is created and returned.

    *title* is the subplot's title.

    If *crop_y* is true (default), the y-axis (frequency) will be
    limited so that only frequency values are shown where all bands are
    known. Alternatively, a numeric value of *crop_y* denotes the upper
    frequency value where the plot will be cropped.

    If *light_cone*, add a light cone and crop the bandgaps at the light
    line.

    """

    # TODO add argument color_by_parity=True, read parity data from files and
    # color the plot lines accordingly.


    if custom_plotter is None:
        plotter = BandPlotter()
    else:
        plotter = custom_plotter
        plotter.next_plot()

    x_axis_formatter = None
    if isinstance(x_axis_hint, axis_formatter.CustomAxisFormatter):
        # use the supplied CustomAxisFormatter:
        x_axis_formatter = x_axis_hint
    elif isinstance(x_axis_hint, KSpace):
        # make a KSpaceAxisFormatter instance from kspace object:
        x_axis_formatter = axis_formatter.KSpaceAxisFormatter(x_axis_hint)
    elif isinstance(x_axis_hint, int):
        # make a standard KVectorAxisFormatter with supplied number of ticks:
        x_axis_formatter = axis_formatter.KVectorAxisFormatter(x_axis_hint)
    else:
        try:
            # is this a sequence?
            hintlen = len(x_axis_hint)
        except TypeError:
            # no sequence
            hintlen = 0
            try:
                # is this a number?
                num = int(x_axis_hint)
            except ValueError:
                # no number
                num = 0
        if hintlen > 1 and (isinstance(x_axis_hint[0], int) and
                hasattr(x_axis_hint[1], 'format')):
            # Supplied a list with at least an int and format_str.
            # Use all items in list as KVectorAxisFormatter arguments:
            x_axis_formatter = axis_formatter.KVectorAxisFormatter(
                *x_axis_hint)
        elif num > 0:
            # make a standard KVectorAxisFormatter with supplied number
            # of ticks:
            x_axis_formatter = axis_formatter.KVectorAxisFormatter(num)
    if x_axis_formatter is None:
        log.warning('draw_bands: Did not understand x_axis_hint, '
            'using default.')
        x_axis_formatter = axis_formatter.KVectorAxisFormatter(
            default_x_axis_hint)

    for i, mode in enumerate(modes):
        fname = '{0}_{1}freqs.csv'.format(jobname, mode)
        data = loadtxt(fname, delimiter=',', skiprows=1)

        # add hover data:
        if x_axis_formatter._hover_func_is_default:
            x_axis_formatter.set_hover_data(data[:, 1:4])

        ''' work in progress... #TODO color by parity
        if color_by_parity:
            # try to load parity files:
            fname = '{0}_{1}yparity.csv'.format(jobname, mode)
            try:
                yparities = loadtxt(fname, delimiter=',')
            except IOError:
                yparities = None
            #print(fname, yparities)
            fname = '{0}_{1}zparity.csv'.format(jobname, mode)
            try:
                zparities = loadtxt(fname, delimiter=',')
            except IOError:
                zparities = None
            #print(fname, zparities)
            # cannot simply color the individual points in a plt.plot.
            # 2 Alternatives:
            # - superimpose a scatter plot (where each point can be
            #   individually colored)
            # - or split band data into 5 categories:
            # (oddz, evenz) x (oddy, eveny) + (everything else)
            # each category has in general multiple sets of connected data
            # points, draw each set with plt.plot in the category's color

        ...end : work in progress'''

        plotter.plot_bands(
            data[:, 5:], data[:, 1:5],
            formatstr='o-',
            x_axis_formatter=x_axis_formatter,
            label=mode.upper(), crop_y=crop_y)
        if light_cone:
            gapbands = get_gap_bands(data[:, 5:], light_line=data[:, 4])
        else:
            gapbands = get_gap_bands(data[:, 5:])
        for band in gapbands:
            plotter.add_band_gap_rectangle(
                band[1], band[2],
                light_line=data[:,4] if light_cone else None)

    if light_cone:
        plotter.add_light_cone()

    plotter.set_plot_title(title)
    plotter.add_legend() 

    return plotter
    
def draw_dos(jobname, modes, custom_plotter=None, title=''):
    """Draw density of states with matplotlib in one subplot. """
    if custom_plotter is None:
        plotter = BandPlotter()
        callnextplot = False
    else:
        plotter = custom_plotter
        callnextplot = True
    
    for i, mode in enumerate(modes):
        fname = jobname + '_dos_' + mode + '.csv' if mode else \
                jobname + '_dos.csv'
        try:
            freqs, dos = loadtxt(fname, delimiter=',', unpack=True)
        except IOError:
            log.error("in graphics.draw_dos: "
                "File not found: {0}\n".format(fname) + 
                "Did you save DOS data in the simulation?")
            return plotter
        if callnextplot:
            callnextplot=False
            plotter.next_plot()
        
        plotter.plot_dos(dos, freqs)
            
    plotter.set_plot_title(title)
        
    return plotter
