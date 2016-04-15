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
from math import sqrt, pi, sin, cos, ceil
from geometry import Geometry
from objects import Rod
from copy import copy
from os import path
from glob import glob1
import re
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from fractions import Fraction

import log

def occupancy_radius(occupancy,n,cell_area = 1.0):
    return sqrt(cell_area*occupancy/n/pi)

def wheel(width,height,n,occupancy,separation,material,priority='None'):
    '''Produces a geometry consisting of n rods of specified material
    with specified cell occupancy in a cell of default area = 1.
    The distance between the center of the cell and the body of each
    rod is adjustable using separation.'''
    distance = lambda point1,point2 : sqrt((point1[0]-point2[0])**2+\
    (point1[1]-point2[1])**2)
    r = occupancy_radius(occupancy,n,width*height)
    if n==1 : return Geometry(width,height,[Rod(0,0,material,r)])
    R = r + separation
    wheel_point = lambda N : (R*cos(N*2*pi/n),R*sin(N*2*pi/n))
    wheel_priority = {\
    'None':lambda R,r,n:(R,r,n)\
    ,'Occupancy':lambda R,r,n:(distance(wheel_point(0),wheel_point(1))>2*r\
    and (R,r,n) or (r/sin(2*pi/(2*n)),r,n))\
    ,'Distance':lambda R,r,n:(distance(wheel_point(0),wheel_point(1))>2*r\
    and (R,r,n) or (R,R*sin(2*pi/(2*n)),n))\
    }
    R,r,n = wheel_priority[priority](R,r,n)
    return Geometry(width,height,[Rod(*wheel_point(N),material=\
    copy(material),radius=r) for N in range(n)])

def max_epsilon(geometry, anisotropic_component=0):
    return max(
        obj.material.epsilon[anisotropic_component] 
        if isinstance(obj.material.epsilon, (list, tuple)) 
        else obj.material.epsilon 
        for obj in geometry.objects)

def get_intersection_freq(freq_left1, freq_right1, freq_left2, freq_right2):
    """Based on two lines, line1 from (n, freq_left1) to (n+1, freq_right1) and 
    line2 from (n, freq_left2) to (n+1, freq_right2), return the frequency 
    (y-value) where they intersect. It is not checked whether they intersect,
    so please take care of that.

    """
    '''y1 = m * x1 + x0
    fr1 = m * 1 + fl1; m = fr1-fl1 /
    fl2 = m * 1 + fr2; m = fl2-fr2 \
    fi = (fr1-fl1)*xi + fl1
    fi = (fl2-fr2)*(1-xi) + fr2 = (fr2-fl2)*xi + fl2
    xi = (fi-fl1)/(fr1-fl1)
    xi = (fi-fl2)/(fr2-fl2)
    (fi-fl1)/(fr1-fl1) = (fi-fl2)/(fr2-fl2)
    (fr2-fl2)*(fi-fl1) = (fr1-fl1)*(fi-fl2)    
    (fr2-fl2)*fi- (fr2-fl2)*fl1 = (fr1-fl1)*fi-(fr1-fl1)*fl2
    (fr2-fl2-fr1+fl1)*fi = (fr2-fl2)*fl1 - (fr1-fl1)*fl2
    fi = ((fr2-fl2)*fl1 - (fr1-fl1)*fl2) / (fr2-fl2-fr1+fl1)
    fi = (fr2*fl1 - fr1*fl2) / (fr2-fl2-fr1+fl1)
    '''
    return ((freq_right2*freq_left1 - freq_right1*freq_left2) / 
            (freq_right2-freq_left2-freq_right1+freq_left1))

def get_intersection_knum(freq_left, freq_right, freq_intersection):
    """Based on two lines, line1 from (0, freq_left) to (1, freq_right) and 
    a horizontal line at freq_intersection, return the knum (x-value) where 
    they intersect. It is not checked whether they intersect, so please take 
    care of that.
    If you use this to get the intersection between to other consecutive knums,
    don't forget to add the left knum to the result.

    """
    #freq_intersection = (freq_right - freq_left) * xi + freq_left
    return (freq_intersection - freq_left) / (freq_right - freq_left)

def get_gap_bands(
        banddata, threshold=5e-4, light_line=None):
    """Calculate the band gaps from the banddata.

    Return the band number after which a gap occurs (first band is band 1),
    the highest frequency of the band just below the gap, the lowest frequency 
    of the band just above the gap and the normalized gap width.

    *banddata* must have shape: (number_of_k_vecs, number_of_bands).
    Gaps smaller than *threshold* will not be counted as gap.
    If *light_line* is given (list of frequency values, one for each k-vector),
    band frequencies higher than the light line frequencies will be ignored.

    """

    bands = []
    # the minimum frequency of each band:
    minfreqs = banddata.min(axis=0)
    # the maximum frequency of each band:
    maxfreqs = banddata.max(axis=0)
    if light_line is not None:
        # A light line is provided. Now, minfreqs and maxfreqs must be the
        # minimum (maximum) of frequencies which are under the light line.
        # So we just go through all k-vecs of each band and throw away all
        # frequencies above the light line, but also linearly interpolate
        # between consecutive k-vecs where one's frequency is above and the
        # other's is below the light line, and include the frequency at the
        # intersection with the light line.

        for bandnum in range(banddata.shape[1]):
            freqs = []
            above = banddata[0, bandnum] > light_line[0]
            for ik, freq in enumerate(banddata[:, bandnum]):
                if freq <= light_line[ik]:
                    if above:
                        # the previous freq was above, so get intersection:
                        freqs.append(get_intersection_freq(
                            freq_left1=light_line[ik-1],
                            freq_right1=light_line[ik],
                            freq_left2=banddata[ik-1, bandnum],
                            freq_right2=banddata[ik, bandnum]))
                        above = False
                    # keep the frequency:
                    freqs.append(freq)
                elif not above:
                    # above light line, so throw away freq, but more important:
                    # we crossed the light line, so add intersection:
                    freqs.append(get_intersection_freq(
                        freq_left1=light_line[ik-1],
                        freq_right1=light_line[ik],
                        freq_left2=banddata[ik-1, bandnum],
                        freq_right2=banddata[ik, bandnum]))
                    above = True

            if len(freqs) > 0:
                maxfreqs[bandnum] = max(freqs)
                minfreqs[bandnum] = min(freqs)
            else:
                maxfreqs[bandnum] = -1
                minfreqs[bandnum] = -1

    for i in range(len(minfreqs) - 1):
        if (minfreqs[i + 1] - maxfreqs[i]) > threshold:
            # the bands are counted from 1:
            bandnum = i + 1
            lo = maxfreqs[i]
            hi = minfreqs[i+1]
            width = 2 * (hi - lo) / (hi + lo)
            bands.append((bandnum, lo, hi, width))
    return bands
    
def strip_format_spec(format_str):
    """Remove all format-specifications from the format-string.

    Changes all format replacement fields in *format_str* containing
    format-specs, e.g.: {0:3.{1}f}, to a field without format-spec, e.g. {0}.

    This is useful if one needs to supply strings to a format_str intended for
    accepting numbers only, e.g. 'x = {0:.3f}'.format('1/3') will not work, but
    strip_format_spec('x = {0:.3f}').format('1/3') will return 'x = 1/3'.

    """
    # Build the regular expression pattern for finding format specifications:

    # re that matches only non-escaped opening curly brackets:
    a = r'(?<!\\){'
    # re that matches anything not containing curly brackets or a colon:
    b = '[^{}:]*'
    # re that matches anything not containing curly brackets:
    c = '[^{}]*'
    # re that is intended to match nested pairs of curly brackets in the
    # format specification, i.e. it matches anything without curly brackets
    # enclosed inside a pair of curly brackets, with an optional suffix of
    # any sign but curly brackets:
    d = ''.join(['(?:{', c, '}', c, ')*'])
    # The final re pattern matches format replacement fields containing
    # format-specs, e.g.: {0:3.{1}f}, and saves the argument number in a
    # group, e.g. in the example '0'.
    re_pattern = ''.join([a, '(', b, '):', c, d, '}'])

    # replace double curly brackets {{ with escaped curly brackets /{:
    fstr = re.sub('{{', r'\\{', format_str)

    # replace {X:...} with {X}, so that format_str accepts strings:
    fstr = re.sub(re_pattern, r'{\1}', fstr)

    # revoke the escaped curly brackets replacement:
    fstr = re.sub(r'\\{', '{{', fstr)

    return fstr


def distribute_pattern_images(
        imgfolder, dstfile, borderpixel=5, only_k_slice=None, 
        title='', show=False):
    """ Read all pngs (from MPB simulation) from imgfolder and distribute
    them according to bandnumber and k vector number in the file(s) dstfile.
    The filenames must be in a format like h.k55.b06.z.zeven.r.png, where
    the order does not matter, but k## and b## must be given.
    borderpixel is the number of pixels that the border around
    the images will take up. (between r and i parts; border between bands and 
    kvecs will take up 3*borderpixel)
    If only_k_slice is None (default) all found images at all k vec numbers
    will be added. Specify a tuple (from, to) to only include those (indices
    into the list of found k-vecs, inclusive).
    """
    if not path.isdir(imgfolder):
        return 0
    # make list of all field pattern h5 files:
    filenames = glob1(imgfolder, "*[edh].*.png")
    if not filenames:
        return 0
    
    log.info("saving field patterns to file(s): %s" % dstfile)
    
    # all files should have same field component, direction and mode,
    # so just read them representatively from first filename:
    first = filenames[0]
    parts = first.split('.')
    reparts = parts[:]
    klabel = [] # the parts of label format of x-axis in plot
    # find where knums and bandnums are in filenames:
    # NOTE: there might be a sophisticated way to do this very nice and short
    # with re.sub (incl. the counting of the number length), but no time to 
    # find out right now.
    btest = re.compile('^b(\d+)$')
    ktest = re.compile('^k(\d+)$')
    irtest = re.compile('^[ri]$')
    tests = [btest, ktest, irtest]
    for i, part in enumerate(parts):
        results = [t.match(part) for t in tests]
        if results[0]:
            # found part with band number
            parts[i] = 'b{bandnum:0%id}' % len(results[0].groups()[0])
            reparts[i] = 'b(?P<bandnum>\d+)'
        elif results[1]:
            # found part with kvec number
            parts[i] = 'k{knum:0%id}' % len(results[1].groups()[0])
            reparts[i] = 'k(?P<knum>\d+)'
            klabel.append(parts[i])
        elif results[2]:
            # found part with r or i (real or imaginary part)
            parts[i] = '{ri}'
            reparts[i] = '(?P<ri>[ri])'
            klabel.append(parts[i])
        
    fname_base = path.join(imgfolder, '.'.join(parts))
    retest = re.compile('.'.join(reparts))
    
    # now get the numbers of bands and kvecs by looking at all files:
    knums = set()
    bnums = set()
    ri = set()    

    for fname in filenames:
        m = retest.match(fname)
        if not m:
            # no match found, maybe some wrong stray file. Warn and ignore:
            log.warning('Warning in distribute_pattern_images: ' \
                  'found non-matching file:', fname, 'in', imgfolder)
            continue
        d = m.groupdict()
        knums.add(int(d['knum']))
        bnums.add(int(d['bandnum']))
        ri.add(d['ri'])
        
    # change sets to sorted lists:
    knums = sorted(knums)
    if only_k_slice is not None:
        knums = knums[only_k_slice[0]:only_k_slice[1] + 1]
    bnums = sorted(bnums)
    ri = sorted(ri, reverse=True) # reverse, because I want real part first
    
    nc = len(ri)
    
    # read img size from first, all images should be the same!   
    img=mpimg.imread(path.join(imgfolder, first))
    imgsize = (img.shape[1], img.shape[0])
    img_aspect = imgsize[1] / imgsize[0]

    # calc pixelsize in data units:
    # I want to force the individual pngs into areas of 1x1 data units,
    # including a half-border around the pngs. That way, the pngs will be
    # placed at integer values of the axes.
    # Because the pngs are generally not rectangular, we will have a different
    # pixelsize in x and y.

    # border belonging to one png in x: (0.5 + 1.5) * bordersize:
    pixelsize_x = 1.0 / (imgsize[0] + 2 * borderpixel)
    # border belonging to one png in y: (1.5 + 1.5) * bordersize:
    pixelsize_y = 1.0 / (imgsize[1] + 3 * borderpixel)
    #print 'pixel sizes (in data space)', pixelsize_x, pixelsize_y
    
    # the aspect ratio for the subplot so the pixels turn out rectangular:
    ax_aspect = pixelsize_x / pixelsize_y
    
    # calc extents of the indivdual pngs:
    # These values are given in data units. They denote the distance between
    # an integer value of an axis (near png center) to where the boundaries
    # of the pngs will be placed in data space, thereby stretching/shrinking
    # the pngs and leaving an empty border between adjacent pngs. 
    ext_thin_border_x = 0.5 - 0.5 * pixelsize_x * borderpixel
    ext_thick_border_x = 0.5 - 1.5 * pixelsize_x * borderpixel
    ext_border_y = 0.5 - 1.5 * pixelsize_y * borderpixel
    
    # size in data units:
    w_dataunits = len(knums) * nc
    h_dataunits = len(bnums)
    
    # now we have all data, so start plotting:
    fig = plt.figure(num=imgfolder, 
                     figsize=(2 * w_dataunits, 2 * h_dataunits))
                     # note: do not play with dpi here, it does not change 
                     # the point size for fonts, so the graphics sizes change
                     # while label sizes stay constant!
    ax = fig.add_subplot(111, axisbg='0.5', aspect=ax_aspect)

    for ib, bandnum in enumerate(bnums):
        for ik, knum in enumerate(knums):
            for ic, comp in enumerate(ri):
                fname = fname_base.format(
                    bandnum=bandnum, knum=knum, ri=comp)
                if not path.isfile(fname):
                    log.warning( 'Warning in distribute_pattern_images: ' \
                          'could not find file:', fname)
                    continue
                x0 = ik * nc + ic
                xl = x0 - ext_thin_border_x if ic else x0 - ext_thick_border_x
                xr = x0 + ext_thick_border_x if ic else x0 + ext_thin_border_x
                y0 = 1 + ib
                img=mpimg.imread(fname)
                ax.imshow(
                    img, 
                    origin='upper', 
                    extent=(xl, xr, y0 - ext_border_y, y0 + ext_border_y),
                    interpolation='none')
                    
    # set aspect; must be done after ax.imshow, as the latter changes it:               
    ax.set_aspect(ax_aspect)
    
    # set ticks, labels etc.:
    kform = '.'.join(klabel)
    xticks = [kform.format(knum=k, ri=c) for k in knums for c in ri]
    ax.set_xticks(range(len(xticks)))
    ax.set_xticklabels(xticks, rotation=45)
    ax.tick_params(which='both', direction='out', length=2)
    ax.set_xlabel('Wave vector index', size='x-large')
    ax.set_ylabel('Band number', size='x-large')
    if title:
        ax.set_title(title, size='x-large')
        
    # choose proper data region:
    # ax.autoscale_view(tight=True)  
    ax.set_xlim(-0.5, w_dataunits - 0.5)
    ax.set_ylim(0.5, 0.5 + h_dataunits)
    
    # width of single png in data units:
    w_png_dataunits = imgsize[0] * pixelsize_x 
    h_png_dataunits = imgsize[1] * pixelsize_y 
    #print 'size of png in data units:', w_png_dataunits, h_png_dataunits
    
    # read here about transformations: 
    # http://matplotlib.org/users/transforms_tutorial.html
    
    # transform sizes in (axis') data units to pixels:
    w_png_currentpixels, h_png_currentpixels = (
        ax.transData.transform([w_png_dataunits, h_png_dataunits]) - 
        ax.transData.transform((0,0)))
    #print 'size of png transformed to pixel:', w_png_currentpixels, \
    #                                           h_png_currentpixels
    
    # transformation to transform pixel sizes to figure units:
    # i.e., what percentage of the whole figure takes up a single png?
    pixel_to_figcoords = fig.transFigure.inverted()    
    w_png_figunits, h_png_figunits = pixel_to_figcoords.transform(
        (w_png_currentpixels, h_png_currentpixels))
    #print 'size of png in figure units:', w_png_figunits, h_png_figunits    
    
    # how many pixels should the whole figure contain, so that the individual
    # pngs have their original resolution?
    w_fig_pixel = imgsize[0] / w_png_figunits
    h_fig_pixel = imgsize[1] / h_png_figunits
    # print 'goal size of whole figure in pixel:', w_fig_pixel, h_fig_pixel
    
    # current size of figure in inches:    
    wfig, hfig = fig.get_size_inches()

    # figure dpi, so that resulting pixel size is same than original images:
    wdpi = w_fig_pixel / wfig
    hdpi = h_fig_pixel / hfig
    #print 'dpi to save figure:', wdpi, hdpi
    # note: I used here that 1 figure unit is 1 inch, but this is not correct
    # since I needed to keep the aspect ratio of the individual images while
    # I forced them (with extent parameter in imshow) on areas of approx. 1x1
    # figure units. Matplotlib then scales everything so it fits, leading to 
    # different x and y DPIs. I believe the biggest DPI is correct, because
    # the smaller axis (smaller figure width in figure units than is actually
    # returned by get_size_inches) will just be padded in the image, leading
    # to a smaller DPI as calculated above.
    dpi = max(wdpi, hdpi); 
    
    if isinstance(dstfile, (list, tuple)):
        for dstfname in dstfile:
            fig.savefig(dstfname, dpi=dpi, bbox_inches='tight', pad_inches=0)
    else:
        fig.savefig(dstfile, dpi=dpi, bbox_inches='tight', pad_inches=0)
    
    if show:
        if show == 'block':
            plt.show(block=True)
        else:
            plt.show(block=False)
    else:
        del fig

def do_runmode(
        sim, runmode, num_processors, bands_plot_title, plot_crop_y,
        convert_field_patterns, field_pattern_plot_k_slice, x_axis_hint):
    """Start a job on the sim object, according to runmode.

    Keyword arguments:

    runmode -- can be one of the following:
        ''       : just create and return the simulation object
        'ctl'    : just write the ctl file to disk
        'sim'    : run the simulation and do all postprocessing
        'postpc' : do all postprocessing; simulation should have run before!
        'display': display all pngs done during postprocessing. This is the
                   only mode that is interactive.

    num_processors   -- the number of processors used for the simulation.
    bands_plot_title -- the title of the band diagrams made in post_processing.
    plot_crop_y      -- the band diagrams are automatically cropped before the
                        last band if plot_crop_y is True, alternatively use
                        plot_crop_y to specify the max. y-value where the plot
                        will be cropped.
    convert_field_patterns --
                        indicates whether field pattern h5 files should be
                        converted to png (only when postprocessing). If this
                        is true, a diagram with all patterns will be created
                        with field patterns for all bands and for the k-vectors
                        included in field_pattern_plot_k_slice
    field_pattern_plot_k_slice --
                        Which k-vecs to include in the field pattern diagram.
                        This slice is a tuple with starting and ending
                        (inclusive) index of the k-vectors where the patterns
                        were exported during simulation, e.g. (0, 2) for the
                        first, second and third exported k-vector.
    x_axis_hint      -- gives a hint on which kind of ticks and labels should
                        be shown on the x-axis of the band diagram and provides
                        the data needed.
                        Can be one of the following:
                        - integer number:
                            The axis' labels will be the 3D k-vectors. The
                            number denotes the number of major ticks and labels
                            distributed on the axis.
                        - list([integer, format-string]):
                            Same as above, but the labels are formatted with
                            the format-string - this gives the possibility to
                            only show one of the three vector components, e.g.
                            the string "{2}" to only show the k-vector's
                            z-component. The axis title will be inferred from
                            the format-string.
                        - KSpace object:
                            This must be a KSpace object created with
                            point_labels. These labels usually denote the high
                            symmetry or crititical points, and they will be
                            shown on the axis.
                        - CustomAxisFormatter object:
                            This gives the possibility to completely customize
                            the x-axis' tick positions, tick labels and axis
                            label. If the CustomAxisFormatter's hover data have
                            not been set, it will be set here with the
                            k-vectors read from the simulation results.
    """
    if not isinstance(runmode, str):
        return sim

    if runmode.startswith('c'): # create ctl file
        sim.write_ctl_file(sim.workingdir)
    elif runmode.startswith('s'): # run simulation
        error = sim.run_simulation(num_processors=num_processors)
        if error:
            return False
        # now continue with postprocessing:
        runmode='postpc'
    if runmode.startswith('p'): # postprocess
        # create csv files of data and pngs:
        sim.post_process(convert_field_patterns=convert_field_patterns)
        # save band diagram as pdf&png:
        sim.draw_bands(
            title=bands_plot_title, crop_y=plot_crop_y,
            x_axis_hint=x_axis_hint)
        # save mode patterns to pdf&png:
        if convert_field_patterns:
            sim.draw_field_patterns(
                title=bands_plot_title,
                only_k_slice=field_pattern_plot_k_slice)
    elif runmode.startswith('d'): # display pngs
        # display png of epsilon:
        sim.display_epsilon()
        # save and show mode patterns in pdf:
        if convert_field_patterns:
            sim.draw_field_patterns(
                title=bands_plot_title,
                only_k_slice=field_pattern_plot_k_slice,
                show=True)
        # show band diagram:
        sim.draw_bands(
            title=bands_plot_title, show=True, crop_y=plot_crop_y,
            x_axis_hint=x_axis_hint, save=False)
    return sim