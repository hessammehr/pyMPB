    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the imKvecFormatterplied warranty of
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
    colorbar_style
from bandplotter import BandPlotter
import objects

def draw_geometry(
        geometry,jobname,format='pdf', display=True, block_when_showing=True,
        anisotropic_component=0):
# TODO
    """WARNING: 
    needs to be adapted to work with anisotropic material
    
    Draw and save the geometry. Save as file jobname+.+format. 
    Only show figure if display=True (default). When showing, block script if 
    block_when_showing=True (default). In case of anisotropic material, only 
    draw anisotropic_component (default 0).
    """
    global maxeps
    # I commented clf(), because it just opens an unneccessary new empty 
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
# TODO
    """WARNING: 
    needs to be adapted to work with anisotropic material
    
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

def draw_bandstructure(
        jobname, mode, kspace, band, ext='.csv', format='pdf', filled=True,
        levels=15, lines=False, labeled=False, legend=False):
    """
    Draw 2D band contour map of one band"""
    #clf()
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111,aspect='equal')
    x,y,z = loadtxt(
        "{0}_{1}{2}".format(jobname, mode, ext),
        delimiter=', ',skiprows=1,usecols=\
    (1,2,4+band),unpack=True)
    if kspace.dimensions == 1:
        plt.plot(x,y,z)
    elif kspace.dimensions ==2:
        xi = np.linspace(-0.5,0.5,kspace.x_res)
        yi = np.linspace(-0.5,0.5,kspace.y_res)
        zi = griddata(x,y,z,xi,yi)
        if filled:
            cs = ax.contourf(xi,yi,zi,levels,**contour_filled)
            legend and plt.colorbar(cs,**colorbar_style)
            cs = lines and ax.contour(xi,yi,zi,levels,**contour_lines)
            labeled and lines and plt.clabel(cs,fontsize=8,inline=1)
        else:
            cs = ax.contour(xi,yi,zi,levels,**contour_plain)
            legend and plt.colorbar(cs,**colorbar_style)
            labeled and plt.clabel(cs,fontsize=8,inline=1)    
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,0.5)
    plt.savefig(jobname,format=format,transparent=True)
    plt.show()

def draw_bands(
        jobname, modes, custom_plotter=None, title='', crop_y=True, 
        light_cone=False):
    """Draw all bands calculated along all k vecs in a band diagram.
    Draw data with matplotlib in one subplot.
    If crop_y is true (default), the y-axis (frequency) will be limited 
    so that only frequency values are shown where all bands are known.
    Alternatively, a numeric value of crop_y denotes the upper frequency
    value where the plot will be cropped.
    
    """
    if custom_plotter is None:
        plotter = BandPlotter()
    else:
        plotter = custom_plotter
        plotter.next_plot()
    
    for i, mode in enumerate(modes):
        fname = jobname + '_' + mode + '.csv' if mode else jobname + '.csv'
        data = loadtxt(fname, delimiter=',', skiprows=1)
        
        plotter.plot_bands(data[:, 1:],
            formatstr = 'o-', x_data_index = -1,
            label=mode.upper(), crop_y=crop_y)
        if light_cone:
            gapbands = get_gap_bands(data[:, 5:], light_line=data[:, 4])
        else:
            gapbands = get_gap_bands(data[:, 5:])
        for band in gapbands:
            plotter.add_band_gap_rectangle_manual(
                band[1], band[2],
                light_line=data[:,4] if light_cone else None)   
        
    if light_cone:
        plotter.add_light_cone()
            
    plotter.set_plot_title(title)
    plotter.add_legend()
        
    return plotter
    
def draw_dos(jobname, modes, custom_plotter=None, title=''):
    """Draw density of states.
    Draw dos-data with matplotlib in one subplot.
    """
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
            print ("File not found: {0}\n".format(fname) + 
                   "Did you save DOS data in the simulation?")
            return plotter
        if callnextplot:
            callnextplot=False
            plotter.next_plot()
        
        plotter.plot_dos(dos, freqs)
            
    plotter.set_plot_title(title)
        
    return plotter
