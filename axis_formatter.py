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


from __future__ import division, print_function

from matplotlib import ticker as mticker
import numpy as np
from fractions import Fraction
import re
import defaults
from utility import strip_format_spec
import log


def infer_k_axis_label_from_format_string(format_str):
    """Given a format_str that is intended to format tick labels, infer
    an axis label for the k-axis.

    :param format_str:
        the format string that is intended to be used to create tick
        labels by format_str.format(*k_vector), as used in
        KVectorAxisFormatter.

    :return:
        String with an axis label for a k-vector axis.

        If it does not succeed, returns an empty string.

    """
    # Edit format_str, which usually accepts floats, so that it accepts
    # strings.
    fstr = strip_format_spec(format_str)

    try:
        kvec = fstr.format('k_x', 'k_y', 'k_z', '|k|')
    except IndexError:
        # format_str contains more indices than 0-3:
        log.warning('KVectorAxisFormatter: Could not infer axis label '
                    'from format_str. Please supply axis_label.')
        return ''

    # simplify:
    kvec = re.sub(r'\(k_x, k_y, k_z\)', r'\\vec{k}', kvec)

    if '$' in kvec:
        return defaults.default_x_axis_label.format(kvec)
    else:
        # add latex math mode if not contained in format_str yet:
        return defaults.default_x_axis_label.format('$' + kvec + '$')


class CustomAxisFormatter(mticker.Formatter):
    def __init__(self, ticks=list(), labels=list(), hover_data=None,
                 axis_label=''):
        """A formatter to set custom ticks on a matplotlib.axis.

        Keyword arguments:
            ticks      : a sequence with the major tick positions
                         (first kvec is at position 0, the next at 1 and
                         so on),
            labels     : a sequence of strings for the major tick labels
                         (must have same length than ticks).
            hover_data : If provided, the corresponding item in this sequence
                         will be shown in the plot's status bar if the mouse
                         hovers over the plot. Must have the same number of
                         entries than total number of indices in x.
                         Useful e.g. if this is a list of all k-vectors.
                         If set to None (default), it will just show the
                         x-position.
                         Alternatively, this can also be a callable function
                         which accepts one argument, the x index, and returns
                         the data to be shown (with a __str__ method).
            axis_label : the label printed underneath the x-axis (a string).
        """

        # make sure both ticks and labels have the same length:
        numticks = len(ticks)
        numlabels = len(labels)
        if numticks > numlabels:
            labels.append([''] * (numticks - numlabels))
        elif numlabels > numticks:
            labels = labels[:numticks]

        self._ticks = ticks
        self._labels = labels
        self._hover_data = None
        self._hover_func = lambda x: x
        self._hover_func_is_default = True
        self.set_hover_data(hover_data)
        self._axis_label = axis_label

    def __call__(self, x, tickindex=None):
        """Return the label for the *tickindex*th tick, located on the axes at
        position *x*.

        Matplotlib uses this for two purposes:
        1.: to retrieve the strings for the tick labels and
        2.: to retrieve a string for describing the x-position of the mouse
            (in the status bar) while hovering over the plot.

        """
        if tickindex is None or tickindex >= len(self._labels):
            # case 2 described above or did not get all the labels in __init__:
            return str(self._hover_func(x))
        else:
            # case 1, mpl is building the ticks labels:
            return self._labels[tickindex]

    def _get_hover_data_from_position(self, x):
        # TODO: here, we could return the linear interpolation between
        # the individual data points, if they are numbers or vectors of
        # numbers, because x is a continuous float, but for now,
        # rounding suffices:
        if x >= 0 and int(x + 0.5) < len(self._hover_data):
            return self._hover_data[int(x + 0.5)]
        else:
            return x

    def set_hover_data(self, hover_data):
        """Set the data that will be shown when the mouse hovers over
        the plot.

        :param hover_data:
            can be a sequence of any objects with a __str__ method, it
            just needs the have the same amount of entries than the
            total number of indices (integers) on the x-axis. Then the
            corresponding item will be shown in the plot's status bar
            when the mouse hovers over the plot. This is useful e.g. if
            this is a list of all k-vectors of the simulation.

            Alternatively, this can also be a callable function which
            accepts one argument, the x index, and returns the data to
            be shown.

            This data can be unset by providing None as argument. In
            that case, the status bar will just show the x-position.

        """
        self._hover_func_is_default = False
        if hover_data is None:
            self._hover_func = lambda x: x
            self._hover_func_is_default = True
        elif callable(hover_data):
            self._hover_func = hover_data
        else:
            self._hover_data = hover_data[:]
            self._hover_func = self._get_hover_data_from_position

    def get_longest_label_length(self):
        """Return the length of the longest string in list of axis labels."""
        if len(self._labels) > 0:
            return max([len(label) for label in self._labels])
        else:
            return 0

    def apply_to_axis(self, axis, **kwargs):
        """Set the tick positions and labels on an axis using this formatter.

        :param axis:
            the matplotlib.axis object this formatter will be applied
            to,
        :param kwargs:
            Any remaining keyword arguments will be forwarded to the
            matplotlib.text.Text objects making up the major labels, so
            they can be formatted.
        :return: None

        """
        # set the minor ticks, one at each simulated k-vector:
        axis.set_minor_locator(mticker.IndexLocator(1, 0))

        # set the tick positions:
        axis.set_ticks(self._ticks)

        # set the tick label strings and label formatting:
        axis.set_ticklabels(self._labels, **kwargs)

        # The formatter in set_major_formatter provides the data that is shown
        # in the status bar of the plot while the mouse is moving in the plot:
        # (This overrides the label strings set in xaxis.set_ticklabels,
        # but set_ticklabels is still necessary to set the formatting;
        # Also, set_ticklabels must be called before, because set_ticklabels
        # overrides the major_formatter.)
        axis.set_major_formatter(self)

        # set the axis label:
        axis.set_label_text(self._axis_label, size='x-large')


class KVectorAxisFormatter(CustomAxisFormatter):
    def __init__(
            self, num_ticks,
            format_str=defaults.default_kvecformatter_format_str,
            hover_data=None, axis_label='',
            fractions=defaults.ticks_fractions):
        """A formatter to set ticks labeled with k-vectors on a
        matplotlib.axis.

        :param num_ticks:
            The number of major ticks to place on the axis
        :param format_str:
            The ticks will be labeled with the format_str, which will be
            formatted with format_str.format(*vector), where 'vector' is
            a sequence taken from hover_data at the tick position index.
        :param hover_data:
            This should be a sequence of all k-vectors included in the
            diagramm, i.e. a sequence of sequences. The tick labels will
            be generated from this data. Therefore each k-vector should
            have at least the same number of entries than needed for
            format_str.

            If hover_data is left at the default value (i.e. None), it
            must be set later with set_hover_data, otherwise no labels
            will be shown.

            The corresponding item in this sequence will also be shown
            in the plot's status bar if the mouse hovers over the plot.
        :param axis_label:
            the label printed underneath the x-axis (a string)
        :param fractions:
            If True, the formatter will try to convert the vector's
            decimal components to simple fractions in the label.
            (default: False)
        :return:
            KVectorAxisFormatter object

        """
        self._num_ticks = num_ticks
        self._format_str = format_str
        self._fractions = fractions
        if not axis_label:
            axis_label = infer_k_axis_label_from_format_string(format_str)

        CustomAxisFormatter.__init__(
            self,
            ticks=[],
            labels=[],
            hover_data=hover_data,
            axis_label=axis_label)

    def _make_fraction_str(self, floatnum):
        """Try to make a fraction string of floatnum.

        :param floatnum:
            can be a number, a string with number or a sequence of these.
        :return:
            If the resulting fraction's denominator is greater than
            defaults.tick_max_denominator, or if the fraction could not
            be created because floatnum was no number, it will just
            return the original floatnum.
            Otherwise, it returns a string with floatnum written as
            fraction or, if floatnum was a sequence, a list of fraction
            strings,

        """
        try:
            l = len(floatnum)
            if hasattr(floatnum, 'isalnum'):
                # it's a string
                l = 0
        except TypeError:
            # no sequence
            l = 0

        if l:
            # it's a sequence:
            return [self._make_fraction_str(comp) for comp in floatnum]
        else:
            # it's a single entry:
            try:
                # Limit to rather high denominator just to remove inaccuracies
                # due to floating point error:
                f = Fraction(floatnum).limit_denominator(100000)
                # But don't make a tick label with such a high denominator.
                # Only return fraction if not too high:
                if f.denominator <= defaults.tick_max_denominator:
                    return str(f)
                else:
                    return floatnum
            except ValueError:
                # could not make fraction: probably bad/unknown string
                return floatnum

    def set_hover_data(self, hover_data):
        """Set the data that will be shown when the mouse hovers over
        the plot and used to create the tick labels.

        :param hover_data:
            must be a sequence of all k-vectors of the simulation, i.e.
            with the same amount of entries than the total number of
            indices (integers) on the x-axis.
        :return: None

        """
        # set the hover data using the parent method:
        CustomAxisFormatter.set_hover_data(self, hover_data)

        if hover_data is None:
            return

        # Now we can update the ticks and labels:
        self._ticks = []
        self._labels = []
        try:
            # here, hover_data can't be a callable:
            axis_length = len(hover_data) - 1
        except TypeError:
            log.error('KVectorAxisFormatter: Could not set ticks, '
                      'hover_data must be a sequence.')
            return
        if axis_length:
            step = max(1, np.floor(axis_length / (self._num_ticks - 1)))
            self._ticks = np.arange(0, axis_length + 1, step)

        vecs = [self._get_hover_data_from_position(x) for x in self._ticks]
        if self._fractions:
            vecs = self._make_fraction_str(vecs)
        for vec in vecs:
            try:
                lbl = self._format_str.format(*vec)
            except IndexError:
                log.warning(
                    'KVectorAxisFormatter: format_str "{0}" does not '
                    'match hover_data: {1}'.format(self._format_str, vec))
                lbl = ''
            except ValueError:
                # vec could be a list of strings, rather than a list of
                # numbers, which is normal if self._fractions is True.
                # But if format_str is intended for numbers, we need to
                # strip the format-spec from it:
                lbl = strip_format_spec(self._format_str).format(*vec)
            self._labels.append(lbl)


class KSpaceAxisFormatter(CustomAxisFormatter):

    _greek_symmetry_point_names = ['Gamma', 'Delta', 'Lambda', 'Sigma']
    _symmetry_point_to_latex = dict(
        [(g, r'$\{0}$'.format(g)) for g in _greek_symmetry_point_names])

    def __init__(self, kspace):
        """A formatter used by matplotlib for the k-vector-axis' ticklabels.

        plot high symmetry points (gamma etc.)

        Must be initialized with a list of ticklabels (with the same number of
        entries as number of ticks set with set_ticks before) and a list with
        all k-vectors (usually with more entries than ticklabels). During
        mouseover, the k-vectors will be shown in the status bar of the plot.

        *labels* is a sequence of strings.
        *kvecdata* is a n x 3 array, where n is the number of k-points
        *kvecdata* can also be None, in which case the status bar will just
        show the x-position during hovering.

        """
        if kspace.has_labels():
            ticks = np.arange(
                0, kspace.count_interpolated(), step=kspace.k_interpolation+1
            )
            labels = [
                self._symmetry_point_to_latex.get(l, '${0}$'.format(l))
                for l in kspace.labels()
            ]
        else:
            log.warning('KSpaceAxisFormatter: KSpace object has no labels! '
                        'k-vec axis will have no ticks.')
            ticks = []
            labels = []
        CustomAxisFormatter.__init__(
            self,
            ticks=ticks,
            labels=labels,
            axis_label=defaults.default_kspace_axis_label)
