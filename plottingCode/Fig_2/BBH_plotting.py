# -*- coding: utf-8 -*-

# Primary authors:
# Ben.Farr@ligo.org, Vivien.Raymond@ligo.org, Will.Farr@ligo.org

# Modifications:
# Christopher.Berry@ligo.org
# 10/2018 Michael Pürrer Michael.Puerrer@ligo.org

"""
Various functions taken from the O1 BBH analysis
See http://trac.ligo.caltech.edu/cbc/browser/cbc/searches/O1/papers/BBH.
Some of these are used in the O1 & O2 catalog event plots notebook.
"""

import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patheffects as PathEffects
from matplotlib.ticker import Formatter

from six import string_types
from scipy.spatial import Delaunay
from scipy.stats import gaussian_kde as kde
import pandas as pd

from bounded_2d_kde import Bounded_2d_kde

# Set the random seed for reproducibility
random_seed = 1234

# https://matplotlib.org/users/customizing.html
# Make the font match the document
rc_params = {'backend': 'pdf',
             'axes.labelsize': 10,
             'axes.titlesize': 10,
             'font.size': 10,
             'legend.fontsize': 10,
             'xtick.labelsize': 10,
             'ytick.labelsize': 10,
             #'font.family': 'sans-serif',
             'font.family': 'serif',
             'font.sans-serif': ['Bitstream Vera Sans'],
             'font.serif': ['Times New Roman'],
             'text.usetex':True
            }

# Figure out figure size to avoid rescaling in the document
column_width = 246.0
inches_per_pt = 1.0/72.27
fig_width = column_width * inches_per_pt
rc_params['figure.figsize'] = (fig_width, fig_width/1.6)

plt.rcParams.update(rc_params)

# Path effect for white outlines
white_outline = [PathEffects.withStroke(linewidth=2., foreground="w")]

#------------------------------------------------------------------------------
# Various functions taken from O1 BBH analysis
#------------------------------------------------------------------------------
def estimate_2d_post(params, post, data2=None,
                     xlow=None, xhigh=None, ylow=None, yhigh=None, transform=None,
                     gridsize=500, **kwargs):
    x = post[params[0]]
    y = post[params[1]]

    if transform is None:
        transform = lambda x: x

    deltax = x.max() - x.min()
    deltay = y.max() - y.min()
    x_pts = np.linspace(x.min() - .1*deltax, x.max() + .1*deltax, gridsize)
    y_pts = np.linspace(y.min() - .1*deltay, y.max() + .1*deltay, gridsize)

    xx, yy = np.meshgrid(x_pts, y_pts)

    positions = np.column_stack([xx.ravel(), yy.ravel()])

    # Calculate the KDE
    pts = np.array([x, y]).T
    selected_indices = np.random.choice(len(pts), len(pts)//2, replace=False)
    kde_sel = np.zeros(len(pts), dtype=bool)
    kde_sel[selected_indices] = True
    kde_pts = transform(pts[kde_sel])
    untransformed_den_pts = pts[~kde_sel]
    den_pts = transform(untransformed_den_pts)

    Nden = den_pts.shape[0]

    post_kde=Bounded_2d_kde(kde_pts, xlow=xlow, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
    den = post_kde(den_pts)
    inds = np.argsort(den)[::-1]
    den = den[inds]

    z = np.reshape(post_kde(transform(positions)), xx.shape)

    return {'xx':xx, 'yy':yy, 'z':z, 'kde':den, 'kde_sel':kde_sel}

#------------------------------------------------------------------------------
def plot_bounded_2d_kde(data, data2=None, contour_data=None, levels=None,
                        shade=False, gridsize=500, cut=3, legend=True,
                        shade_lowest=True, ax=None, verbose=False, **kwargs):

    if ax is None:
        ax = plt.gca()

    data = data.astype(np.float64)
    if data2 is not None:
        data2 = data2.astype(np.float64)

    bivariate = False
    if isinstance(data, np.ndarray) and np.ndim(data) > 1:
        bivariate = True
        x, y = data.T
    elif isinstance(data, pd.DataFrame) and np.ndim(data) > 1:
        bivariate = True
        x = data.iloc[:, 0].values
        y = data.iloc[:, 1].values
    elif data2 is not None:
        bivariate = True
        x = data
        y = data2

    if bivariate==False:
        raise TypeError("Bounded 2D KDE are only available for"
                        "bivariate distributions.")

    xx = contour_data['xx']
    yy = contour_data['yy']
    z = contour_data['z']
    den = contour_data['kde']
    kde_sel = contour_data['kde_sel']
    Nden = len(den)

    # Black (thin) contours with while outlines by default
    kwargs['linewidths'] = kwargs.get('linewidths', 1.)
    kwargs['linestyles'] = kwargs.get('linestyles', 1.)

    # Plot the contours
    n_levels = kwargs.pop("n_levels", 10)
    cmap = kwargs.get("cmap", None)

    if cmap is None:
        kwargs['colors'] = kwargs.get('colors', 'k')
    if isinstance(cmap, string_types):
        if cmap.endswith("_d"):
            pal = ["#333333"]
            pal.extend(sns.palettes.color_palette(cmap.replace("_d", "_r"), 2))
            cmap = sns.palettes.blend_palette(pal, as_cmap=True)
        else:
            cmap = plt.cm.get_cmap(cmap)

    kwargs["cmap"] = cmap
    contour_func = ax.contourf if shade else ax.contour
    if levels:
        pts=np.array([x, y]).T
        den_pts = pts[~kde_sel]
        zvalues = np.empty(len(levels))
        for i, level in enumerate(levels):
            ilevel = int(np.ceil(Nden*level))
            ilevel = min(ilevel, Nden-1)
            zvalues[i] = den[ilevel]
        zvalues.sort()

        cset = contour_func(xx, yy, z, zvalues, **kwargs)

        for i, coll in enumerate(cset.collections):
            level = coll.get_paths()[0]
            contour_hull = Delaunay(level.vertices)
            if verbose:
                print("{}% of points found within {}% contour".
                format(int(100*np.count_nonzero(contour_hull.find_simplex(den_pts) > -1)/Nden),
                       int(100*levels[::-1][i])))

        # Add white outlines
        #if kwargs['colors'] == 'k':
        #    plt.setp(cset.collections,
        #             path_effects=[PathEffects.withStroke(linewidth=1.5, foreground="w")])
        #fmt = {}
        #strs = ['${}\%$'.format(int(100*level)) for level in levels[::-1]]
        #for l, s in zip(cset.levels, strs):
        #    fmt[l] = s

        #plt.clabel(cset, cset.levels, fmt=fmt, fontsize=11, manual=False, **kwargs)
        #plt.setp(cset.labelTexts, color='k', path_effects=[PathEffects.withStroke(linewidth=1.5, foreground="w")])
    else:
        cset = contour_func(xx, yy, z, n_levels, **kwargs)

    if shade and not shade_lowest:
        cset.collections[0].set_alpha(0)

    # Avoid gaps between patches when saving as PDF
    if shade:
        for c in cset.collections:
            c.set_edgecolor("face")

    kwargs["n_levels"] = n_levels

    # Label the axes
    if hasattr(x, "name") and legend:
        ax.set_xlabel(x.name)
    if hasattr(y, "name") and legend:
        ax.set_ylabel(y.name)

    return ax

#------------------------------------------------------------------------------
def plot_2d(params, post, contour_data, labels=None,
            cmap='Purples', gridsize=500, levels=[.5, .9], **kwargs):
    """ **kwargs are passed to plot_bounded_2d_kde """
    if labels is None:
        labels = [None for p in params]

    post_series = [pd.Series(post[p], name=label) for p, label in zip(params, labels)]

    #ax = plot_bounded_2d_kde(post_series[0], contour_data, cmap=cmap, shade=True, shade_lowest=False, n_levels=30, **kwargs)
    ax = plot_bounded_2d_kde(post_series[0], post_series[1], contour_data=contour_data, 
                             levels=levels, gridsize=gridsize, **kwargs)

    return ax

#------------------------------------------------------------------------------
def pickle_contour_data(event, post_name, contour_dir, cdata):
    outfile = os.path.join(contour_dir, '{}_{}_contour_data.pkl'.format(event, post_name))
    with open(outfile, 'wb') as outp:
        pickle.dump(cdata, outp)

#------------------------------------------------------------------------------
def unpickle_contour_data(event, post_name, contour_dir):
    infile = os.path.join(contour_dir, '{}_{}_contour_data.pkl'.format(event, post_name))
    with open(infile, 'rb') as inp:
        cdata = pickle.load(inp)
    return cdata

#------------------------------------------------------------------------------
def ms2q(pts):
    """Transformation function from component masses to chirp mass and mass ratio"""
    pts = np.atleast_2d(pts)

    m1 = pts[:, 0]
    m2 = pts[:, 1]
    assert (m1 > 0).all() and (m2 > 0).all()
    mc = np.power(m1 * m2, 3./5.) * np.power(m1 + m2, -1./5.)
    if not np.isfinite(mc).all():
        idx_bad = np.where(np.logical_not(np.isfinite(mc)))
        print('ms2q(): infinite or nan value encountered in chirp mass')
        print('idx:', idx_bad)
        print('mc:', mc[idx_bad])
        print('m1:', m1[idx_bad])
        print('m2:', m2[idx_bad])
    q = m2/m1
    return np.column_stack([mc, q])

#------------------------------------------------------------------------------
class DegreeFormatter(Formatter):
    def __call__(self, x, pos=None):
        return r"${:3.0f}^{{\circ}}$".format(x)


