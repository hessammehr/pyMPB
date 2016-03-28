    #Copyright 2014-2016 Juergen Probst
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

from __future__ import division, print_function

import matplotlib as mpl
from matplotlib import pyplot as plt
if mpl.__version__ >= '1.4':
    from matplotlib.axes._base import _process_plot_format
else:
    from matplotlib.axes import _process_plot_format
#from matplotlib.widgets import Slider, RadioButtons
import numpy as np
from fractions import Fraction
from itertools import cycle
from utility import get_intersection_knum
import log

"""last edit 2016-01-24:
- added add_band_gap_rectangle
- BandPlotter can be created without knowing the number of plots beforehand
- added automatic color cycle, so multiple calls to plot_bands lead to
  multiple colors, one color for each call to plot_bands
- add_band_gap_rectangle uses last plotted color if no color is given
- possibility to automatically crop y before hightest band
- added reuse_prev_fig (default: True) to constructor parameters.
- moved crop_y from __init__ to plot_bands
- possiblity to specify an upper y-value where to crop manually
- added add_band_gap_manual and stop_at_light_line parameter
- removed self._extra_y_padding, too chaotic
"""

g_fundlength = 270 # fundamental lengthscale (i.e. a) in nm
g_speed_of_light = 299792458 # m/s

class KvecFormatter(mpl.ticker.Formatter):
    """A formatter for the x-axis' ticklabels. Must be initialized with a list
    of ticklabels (usually high-symmetry points like Gamma, K etc.) and the
    k vector data of all k-points. During mouseover, this data will be shown
    in status bar of plot.

    """
    def __init__(self, labels, kvecdata):
        """
        *labels* is a sequence of strings.
        *kvecdata* is a n x 3 array, where n is the number of k-points

        """
        self.labels = labels
        self.kvecdata = kvecdata

    def __call__(self, x, pos=None):
        'Return the format for tick val *x* at position *pos*'
        # pos is not None for the ticks labels,
        # pos is None for coords output during mouseover
        if pos is None or pos >= len(self.labels):
            return self._get_kvec_from_index(int(x + 0.5))
        else:
            return self.labels[pos]

    def _get_kvec_from_index(self, index):
        return self.kvecdata[index]

class BandPlotter:
    _kz = 0.25
    _directions = {
        (0, 0, 0): r'$\Gamma$',
        (0, 0.5, 0): 'M',
        (-0.333333, 0.333333, 0): 'K',
        (0, 0, _kz): r'$\Gamma$_z',
        (0, 0.5, _kz): 'M_z',
        (-0.333333, 0.333333, _kz): 'K_z'}

    def __init__(
            self, fundlength=g_fundlength, figure_size=(12, 9), numrows=1,
            figure_name='bands'):
        """ Prepare a new figure for plotting photonic bands.
        fundlength is the real length of th fundamental length unit a in the
        MPB simulation, usually the lattice constant;
        figure_size is the figure size in inches (x,y-tuple);
        numrows is the number of subplot rows in the figure.
        Provide a unique figure_name to create a new figure, or reuse the
        same figure_name for multiple plots to reuse figures.

        """
        self._fig = plt.figure(figure_name, figsize = figure_size)
        self._fig.clf()
        self._fig.canvas.mpl_connect('pick_event', self._onpick)
        self._numplots = 0
        self._numrows = max(numrows, 1)
        self._axes = []
        # fundamental lengthscale (i.e. a) in nm (only used in _onpick):
        self._fundlength = fundlength
        self._crop_y = False
        self._crop_y_val = None
        self.next_plot()

    def set_num_rows(self, numrows):
        self._numrows = max(numrows, 1)
        self._distribute_subplots()

    def _get_direction_string(self, vec3):
        t = tuple(vec3)
        if t in self._directions:
            return self._directions[t]
        else:
            return ''

    def _set_x_ticks_labels(self, banddata, num_labels=-1):
        """ if num_labels==-1, plot high symmetry points (gamma etc.)"""
        labels = []
        ticks = []
        num_minors = len(banddata)
        if num_labels == -1:
            for i, kdata in enumerate(banddata):
                dstr = self._get_direction_string(kdata[0:3])
                if dstr:
                    labels.append(dstr)
                    ticks.append(i)
        elif num_labels > 0:

            num = len(banddata) - 1
            step = np.floor(num / (num_labels-1))
            ticks = np.arange(0, num+1, step)
            for i in ticks:
                f = Fraction(banddata[i][0]).limit_denominator(100000)
                if f.denominator > 1000:
                    kxstr = '{0:.10f}'.format(banddata[i][0])
                else:
                    kxstr = str(f)
                f = Fraction(banddata[i][1]).limit_denominator(100000)
                if f.denominator > 1000:
                    kystr = '{0:.10f}'.format(banddata[i][0])
                else:
                    kystr = str(f)
                labels.append(r'$\binom{%s}{%s}$' % (kxstr, kystr))
                #labels.append(r'$\left(\stackrel[3]{{{1}}}{{ \stackrel{{{1}}}{{{0}}} }}\right)$'.format(i, i+1))
                #labels.append('{1}\n{1}\n{0}'.format(i, i+1))
        xaxis = self._ax.get_xaxis()
        xaxis.set_label_text(r'Wave vector', size='x-large')
        xaxis.set_minor_locator(
            mpl.ticker.LinearLocator(num_minors))
        xaxis.set_ticks(ticks)
        xaxis.set_major_formatter(KvecFormatter(labels, banddata[:, 0:3]))

    def plot_bands(
            self, banddata, x_data_index = -1, formatstr = '',
            crop_y=True, picker = 3, label = None, **keyargs):
        """
        Plot bands. plt.show() must be called to actually show the figure
        afterwards.
        The first four columns of banddata must correspond to kx-, ky-,
        kz-components and magnitude of k, respectively. The remaining columns
        will be the plotted bands.
        If x_data_index == -1 (default), plots against the indices of banddata.
        High symmetry points will be labeled.
        x_data_index of 0, 1, 2 and 3 corresponds to kx-, ky-, kz-components
        and magnitude of k, respectively.
        If crop_y is true (default), the y-axis (frequency) will be limited
        so that only frequency values are shown where all bands are known.
        Alternatively, a numeric value of crop_y denotes the upper frequency
        value where the plot will be cropped.

        """
        self._crop_y = bool(crop_y)
        if self._crop_y:
            if crop_y is True:
                self._crop_y_val = banddata[:,-1].min()
            else:
                self._crop_y_val = crop_y

        return self.plot_general(
            banddata[:, 4:], x_data_index, banddata[:, :4], formatstr,
            picker, label, **keyargs)

    def plot_general(
            self, data, x_data_index = -1, k_data = None, formatstr = '',
            picker = 3, label = None, **keyargs):
        """
        Plot data. plt.show() must be called to actually show the figure
        afterwards.
        if x_data_index == -1 (default), plots against the indices of data.
        In this case, high symmetry points will only be labeled if k_data is
        given, containing columns corresponding to kx-, ky-, kz-components
        and magnitude of k, respectively.
        x_data_index of 0, 1, 2 and 3 corresponds to kx-, ky-, kz-components
        and magnitude of k, respectively, which must be given in k_data.
        k_data and data must have same amount of rows.

        """
        if len(data) == 0:
            return
        if x_data_index < 0:
            self._x_data = np.arange(len(data))
            if k_data is not None:
                self._set_x_ticks_labels(
                    k_data,
                    num_labels=(-x_data_index if x_data_index < -1 else -1))
        else:
            if k_data is None:
                raise ValueError('please give k_data')
            self._x_data = k_data[:, x_data_index]
            if x_data_index < 3:
                self._ax.set_xlabel(
                    r'Wave vector $k_{0}a/2\pi$'.format(
                        ['x', 'y', 'z'][x_data_index]),
                    size='x-large')

        # keep reference to last data for add_light_cone and fill_between_bands:
        self._last_data = data
        self._last_kdata = k_data

        if self._crop_y_val is not None:
            if self._maxy == -float('inf'):
                # first call top plot_general in this subplot:
                self._maxy = self._crop_y_val
            else:
                # not first call, maybe prev. maximum was less:
                self._maxy = min(self._maxy, self._crop_y_val)
        else:
            # automatic cropping to maximum data:
            self._maxy = max(self._maxy, data.max())

        self._miny = min(self._miny, data.min())

        if (not 'color' in keyargs and not 'c' in keyargs and
                _process_plot_format(formatstr)[-1] is None):
            # if no color for the plot is given, use automatic coloring:
            keyargs['color'] = next(self._colors)

        if 'c' in keyargs:
            self._last_color = keyargs['c']
        else:
            self._last_color = keyargs['color']

        self._ax.plot(
            self._x_data, data, formatstr, label=label, **keyargs)

        if picker:
            # add invisible dots for picker
            # (if picker added to plots above, the lines connecting the dots
            # will fire picker event as well) - this is ok if no dots are shown
            # Also, combine data so it is in one single dataset. This way,
            # the event will fire only once, even when multiple dots coincide:
            # (but then with multiple indexes)
            xnum, bands = data.shape
            newdata = np.zeros((2, xnum * bands))
            for i, x in enumerate(self._x_data):
                for j, y in enumerate(data[i, :]):
                    newdata[:, i + xnum*j] = [x, y]
            frmt = 'o'
            if not ('o' in formatstr or '.' in formatstr):
                frmt += '-'
            self._ax.plot(
                newdata[0], newdata[1], frmt, picker=picker, alpha = 0,
                label=label)


        # matplotlib sometimes adds padding; remove it:
        self._ax.set_xlim(min(self._x_data), max(self._x_data))
        if self._crop_y_val is not None:
            self._ax.set_ylim(self._miny, self._maxy)

    def plot_dos(self, dos, freqs):
        self._ax.plot(dos, freqs)
        if self._crop_y_val is not None:
            self._ax.set_ylim(0, self._crop_y_val)
        self._ax.set_xlabel('DOS', size='x-large')
        self._ax.set_ylabel('')

    def _onpick(self, event):
        try:
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            xaxisformatter = event.mouseevent.inaxes.xaxis.major.formatter
        except AttributeError:
            return
        print(thisline.get_label() + ' mode(s): ', end='')
        if isinstance(xaxisformatter, KvecFormatter):
            for i in ind:
                print('k_index={0:.0f}, freq={1}; '.format(
                    xdata[i] + 1, ydata[i]), end='')
        else:
            for i in ind:
                omega = ydata[i] * g_speed_of_light * 2.0 * np.pi / self._fundlength
                wlen = self._fundlength / ydata[i]
                print('k={0}, freq={1}, omega={2} Hz, wlen={3} nm; '.format(
                        xdata[i], ydata[i], omega, wlen), end = '')

        #plt.ion()
        #fig = plt.figure('mode pattern', figsize=(6, 2))
        #ax = fig.add_subplot(121) #put mode pattern image here
        #ax = fig.add_subplot(122)

        # fmeep = w a / 2 pi c
        # w = fmeep 2 pi c / a
        # wlen = c / freq = c 2 pi / omega = c 2 pi / (2 pi c fmeep / a) = a / fmeep
        print()

    def add_light_cone(
            self, color = 'gray', alpha=0.5):
        if self._last_kdata is None:
            raise ValueError(
                'cannot add light cone: '
                'k_data not given in last plot_general()')
        fillto = max(self._last_kdata[:, 3].max() * 1.1, self._maxy)
        self._ax.fill_between(
            self._x_data, self._last_kdata[:, 3],
            fillto,
            color=color, alpha=alpha)
        self._ax.plot(self._x_data, self._last_kdata[:, 3], color=color)
        # matplotlib sometimes adds padding; remove it:
        self._ax.set_ylim(self._miny, self._maxy)

    def add_band_gap_rectangle(
            self, after_band_num, color=None, alpha=0.5,
            stop_at_light_line=False):
        if self._last_data is None:
            raise ValueError(
                'cannot add band gap rectangle: '
                'data not given in last plot_general()')
        if color is None:
            #use color of last plotted data:
            color = self._last_color
        y_bottom = max(self._last_data[:,after_band_num - 1])
        y_top = min(self._last_data[:,after_band_num])
        return self.add_band_gap_rectangle_manual(
            y_bottom, y_top, color, alpha,
            self._last_kdata if stop_at_light_line[:, 3] else None)

    def add_band_gap_rectangle_manual(
            self, from_freq, to_freq, color=None, alpha=0.35,
            light_line=None):
        if from_freq < 0 or to_freq <= 0:
            return
        if color is None:
            #use color of last plotted data:
            color = self._last_color
        h = to_freq - from_freq
        center = (to_freq + from_freq) / 2
        gapsize = h / center
        gaptext='gap size: {0:.2f}%'.format(gapsize * 100)
        if light_line is None:
            x_left = self._x_data[0]
            x_right = self._x_data[-1]
            w = x_right - x_left
            self._ax.add_patch(
                mpl.patches.Rectangle((x_left, from_freq), width=w, height=h,
                                      color=color, alpha=alpha))
            self._ax.add_artist(
                mpl.text.Text(text=gaptext, x=(x_left + x_right)/2, y=center,
                              horizontalalignment='center',
                              verticalalignment='center'))
        else:
            pts = []
            # find out at which knum light_line crosses from_freq:
            knum = 0
            while light_line[knum] <= from_freq:
                knum += 1
            # now the light line crossed the lower frequency
            lol = get_intersection_knum(
                light_line[knum-1], light_line[knum], from_freq) + knum - 1
            pts.append((lol, from_freq))
            # now add points of light line until we cross upper frequency;
            # or lower freq, then light line does not cross upper freq:
            while light_line[knum] < to_freq and light_line[knum] > from_freq:
                pts.append((knum, light_line[knum]))
                knum += 1
            if (light_line[knum] > from_freq):
                # the light line crossed the upper frequency
                hil = get_intersection_knum(
                    light_line[knum-1], light_line[knum], to_freq) + knum - 1
                pts.append((hil, to_freq))
                # continue until we enter bandgap again:
                while light_line[knum] >= to_freq:
                    knum += 1
                # now the light line crossed the upper frequency again
                hir = get_intersection_knum(
                    light_line[knum-1], light_line[knum], to_freq) + knum - 1
                pts.append((hir, to_freq))
                # now add points of light line until we leave band gap:
                while light_line[knum] > from_freq:
                    pts.append((knum, light_line[knum]))
                    knum += 1
                # now the light line crossed the lower frequency
                lor = get_intersection_knum(
                    light_line[knum-1], light_line[knum], from_freq) + knum - 1
                pts.append((lor, from_freq))
            else:
                # the light line crossed the lower frequency before the upper
                hil = lol
                hir = lor = get_intersection_knum(
                    light_line[knum-1], light_line[knum], from_freq) + knum - 1
                pts.append((hil, from_freq))

            self._ax.add_patch(
                mpl.patches.Polygon(pts, color=color, alpha=alpha))
            self._ax.add_artist(
                mpl.text.Text(text=gaptext,
                              x=(hil+lol+hir+lor)/4, y=center,
                              horizontalalignment='center',
                              verticalalignment='center'))

    def fill_between_bands(
            self, bandfrom, bandto, color = '#7f7fff', alpha=0.5):
        self._ax.fill_between(
            self._x_data, self._last_data[:, bandfrom - 1],
            self._last_data[:, bandto - 1],
            color=color, alpha=alpha)
        # matplotlib sometimes adds padding; remove it:
        self._ax.set_ylim(self._miny, self._maxy)

    def _distribute_subplots(self):
        rows = max(self._numrows, 1)
        if self._numrows == 0:
            numcols = 1
            rows = self._numplots
        else:
            numcols = int(np.ceil(self._numplots / rows))
            rows = min(rows, self._numplots)

        for i, ax in enumerate(self._axes):
            ax.change_geometry(rows, numcols, i + 1)

        self._fig.tight_layout()

    def next_plot(self):
        self._numplots += 1
        # add subplot somewhere where no other is yet. Final position
        # will be set in _distribute_subplots
        self._ax = self._fig.add_subplot(1, self._numplots, self._numplots)
        self._ax.set_xlabel(r'Wave vector $ka/2\pi$', size='x-large')
        self._ax.set_ylabel(r'Frequency $\omega a/2\pi c$', size='x-large')
        self._ax.grid(True)
        self._axes.append(self._ax)
        self._maxy = -float('inf')
        self._miny = float('inf')
        # start new plot at blue again:
        self._colors = cycle('bgrcmky')
        self._distribute_subplots()

    def add_legend(self, loc = 'upper left'):
        handles, labels = self._ax.get_legend_handles_labels()
        filteredlabels = []
        filteredhandles = []
        # filter labels, so each one occurs only once
        # (otherwise legend is too crowded)
        for i, label in enumerate(labels):
            if not label in filteredlabels:
                filteredlabels.append(label)
                filteredhandles.append(handles[i])
        plt.legend(filteredhandles, filteredlabels, loc=loc)

    def set_plot_title(self, title):
        self._ax.set_title(title, size='x-large')
        self._fig.tight_layout()

    def savefig(self, *args, **kwargs):
        return self._fig.savefig(*args, **kwargs)

    def show(self, block=True):
        plt.show(block=block)


def main():
    from MPBOutputReader import MPBOutputReader
##    reader = MPBOutputReader(g_filename5)
##    print reader.get_band_types()
##    print reader.get_data(reader.get_band_types()[0]).shape
##    reader.print_raw_lines()
##    return

    g_filename = '../sims/pympbtest/TriHoles2D_test.out'
    g_filename2 = '../sims/151116_TriPhC/2D_1BZ/out'

    g_fnames = [g_filename, g_filename2]
    g_xdataindices = [-1, -1]
    numplots = 2

    global plotter
    plotter = BandPlotter()
    for i, fname in enumerate(g_fnames[:numplots]):
        reader = MPBOutputReader(fname)

        if reader.has_multiple_bandtypes():
            plotter.plot_bands(reader.get_data('tm'),
                formatstr = 'bo-', x_data_index = g_xdataindices[i],
                label='TM')
            plotter.add_band_gap_rectangle([3, 3][i], color='blue')
            plotter.plot_bands(reader.get_data('te'),
                formatstr = 'ro-', x_data_index = g_xdataindices[i],
                label='TE')
            plotter.add_band_gap_rectangle([1,1][i], color='red')
##            plotter.fill_between_bands(3, 12, alpha = 0.7)
        else:
            plotter.plot_bands(reader.get_data(),
                x_data_index = g_xdataindices[i],
                label=reader.get_band_types()[0])
            plotter.add_band_gap_rectangle(1)


        #plotter.add_light_cone(alpha = 0.7)
##        if i == 1:
##            plotter._ax.set_xlim(0.4, 0.5)
##            plotter._ax.set_ylim(0.3, plotter._maxy)


        plotter.set_plot_title(
            ['this folder',
            'trapezoidal SiN holes r0.375 2D'][i])
        plotter.add_legend()
        if i < numplots - 1:
            plotter.next_plot()
    plt.show()


if __name__ == '__main__':
    main()
##    W1SiN()


    # test:
##    reader = MPBOutputReader('./sims/140319_W1slab/fullfine/out')
##    reader.save_as_CSV('testout.csv', bandtype = 'zevenyodd')
##    reader.save_band_to_CSV('W1_TE_yodd_band13_fine.csv', 13, 'zevenyodd')

