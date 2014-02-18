# This file is part of Python GRN implementation.
# Architype is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# Architype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License 
# along with GRN.  If not, see <http://www.gnu.org/licenses/>.
# Author Jonathan Byrne 2014

"""A module for plotting results"""

import pylab, sys, os
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

FILETYPE = '.pdf'

def plot_3d(results_list, title, colors=None):
    """show all results in parallel"""
    x_range = range(len(results_list[0]))
    fig = plt.figure()
    axe = Axes3D(fig)
    plt.title(title)
    for idx, result in enumerate(results_list):
        if colors == None:
            axe.plot(x_range, result, idx)
        else:
            axe.plot(x_range, result, idx, color=colors[idx])
    plt.show()

def plot_2d(results_list, title, colors=None, limits=None):
    """multiple runs single graph"""
    pylab.clf()
    pylab.figure().autofmt_xdate()
    x_range = range(len(results_list[0]))
    for idx, result in enumerate(results_list):
        if colors == None:
            pylab.plot(x_range, result, linewidth=0.2)
        else:
            pylab.plot(x_range, result, colors[idx], linewidth=0.2)

    # always set graphs to same axis
    if limits != None:
        pylab.ylim(limits[0],limits[1])

    pylab.title(title)
    if not os.path.exists('./graphs'): os.makedirs('./graphs')
    title = 'graphs/' +title + FILETYPE
    pylab.savefig(title)

def boxplot_data(results_list, title):
    pylab.clf()
    pylab.figure(1)
    result_cols = []
    for i in range(len(results_list[0])):
        res = [result[i] for result in results_list]
        result_cols.append(res)
    pylab.boxplot(result_cols)
    pylab.figure(1).autofmt_xdate()
    title = title + '_boxplot'
    pylab.title(title)
    if not os.path.exists('./graphs'): os.makedirs('./graphs')
    filename  = 'graphs/' + title + FILETYPE
    pylab.savefig(filename)

def plot_ave(results_list, title):
    """ show average with error bars"""
    pylab.clf()
    pylab.figure().autofmt_xdate()

    x_range = range(len(results_list[0]))
    err_x, err_y, std_list = [], [], []

    for i in x_range:
        if i % 10 == 0:
            #get average for each generation
            column = []
            for result in results_list:
                column.append(result[i])
            average = np.average(column)

            std_dev = np.std(column)
            err_x.append(i)
            err_y.append(average)
            std_list.append(std_dev)
    pylab.errorbar(err_x, err_y, yerr=std_list)
    title += '_average'
    pylab.title(title)
    if not os.path.exists('./graphs'): os.makedirs('./graphs')
    filename  = 'graphs/' + title + FILETYPE
    pylab.savefig(filename)
