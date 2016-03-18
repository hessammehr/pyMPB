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
    
    """
    #freq_intersection = (freq_right - freq_left) * xi + freq_left
    return (freq_intersection - freq_left) / (freq_right - freq_left)

def get_gap_bands(
        banddata, threshold=5e-4, light_line=None):
    """Return the band number after which a gap occurs (first band is band 1), 
    the highest frequency of the band just below the gap, the lowest frequency 
    of the band just above the gap and the normalized gap width.
    banddata must have shape (number_of_k_vecs, number_of_bands).
    If light_line is given (list of frequency values, one for each k-vector),
    band frequencies higher than the light line frequencies will be ignored.
    
    """
    bands = []
    # the minimum frequency of each band:
    minfreqs = banddata.min(axis=0)
    # the maximum frequency of each band:
    maxfreqs = banddata.max(axis=0)
    if light_line is not None:
        numk = banddata.shape[0]
        for bandnum in range(banddata.shape[1]):
            freqs = []
            knum = 0
            # increase knum until we cross the light line:
            while knum < numk and banddata[knum, bandnum] > light_line[knum]:
                knum += 1
            if knum == numk:
                # this band is completely above the light line
                maxfreqs[bandnum] = -1
                minfreqs[bandnum] = -1
                continue
            # now, we are below light line:
            # first, add intersection to list:
            if knum > 0:
                freqs.append(get_intersection_freq(
                    freq_left1=light_line[knum-1], 
                    freq_right1=light_line[knum], 
                    freq_left2=banddata[knum-1, bandnum], 
                    freq_right2=banddata[knum, bandnum]))
            # now, add freqs of knum and following to list, until we cross
            # the light line again:
            while knum < numk and banddata[knum, bandnum] <= light_line[knum]:
                freqs.append(banddata[knum, bandnum])
                knum += 1
            # we are above the light line again, so add intersection:
            if knum < numk: # only if we did not leave brillioun zone
                freqs.append(get_intersection_freq(
                    freq_left1=light_line[knum-1], 
                    freq_right1=light_line[knum], 
                    freq_left2=banddata[knum-1, bandnum], 
                    freq_right2=banddata[knum, bandnum]))
            
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
    assert(path.isdir(imgfolder))
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
        
    