"""
Functions designed for plotting output files
"""
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.table import Table


def plot_spec(col, spectra, names, log_x=False, log_y=False, scale_to=None, lim_x=False):
    """
    Plots an array of python spectra imported as Tables.

    Use as:
        plot_spec('A40P0.50', [spec_r01, spec_r10], ['1x', '10x'])
    """
    fig, ax = plt.subplots()
    ax.set_xlabel(r"Wavelength (\AA)")

    if scale_to:
        ax.set_ylabel(r"$F_{\lambda}/F_{"+str(scale_to)+r"\AA}$")
    else:
        ax.set_ylabel(r"Flux")
        scale_factor = 1

    maxY = -1
    minY = 99999999

    for spectrum, name in zip(spectra, names):
        if scale_to:
            # Annoyingly, the file is sorted by frequency not wavelength
            # so we need to flip it to get searchsorted to run
            index_scale = np.searchsorted(spectrum['Lambda'][::-1], scale_to)
            scale_factor = spectrum[col][::-1][index_scale]

        if lim_x:
            minY_curr = np.amin(spectrum[col][(spectrum['Lambda'] > lim_x[0]) & (spectrum['Lambda'] < lim_x[1])]) / scale_factor
            if minY_curr < minY:
                minY = minY_curr
            maxY_curr = np.amax(spectrum[col][(spectrum['Lambda'] > lim_x[0]) & (spectrum['Lambda'] < lim_x[1])]) / scale_factor
            if maxY_curr > maxY:
                maxY = maxY_curr

        if log_x and log_y:
            ax.loglog(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)
        elif log_y:
            ax.semilogy(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)
        elif log_x:
            ax.semilogx(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)
        else:
            ax.plot(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)

    if lim_x:
        ax.set_xlim(lim_x[0], lim_x[1])
        ax.set_ylim(minY, maxY)

    ax.legend()
    return fig


def load_grid(filename):
    """
    Loads a pair of grid files from a root name.

    Use as:
        x_r10, z_r10 = load_grid('r10/')
    """
    x = np.loadtxt(filename+'grid_x.txt')
    z = np.loadtxt(filename+'grid_z.txt')
    x[0] = x[1] * 0.1
    z[0] = z[1] * 0.1
    x[-1] = x[-2] * 10
    z[-1] = z[-2] * 10
    return [x, z]


def plot_dat(table, grid_x, grid_z, title, label, volume=True):
    """
    Plots a given .dat file

    Use as:
        plot_dat(table_h1_r01, x_r01, z_r01, 'H-I, radius 1x', 'Log ion fraction', volume=False)
    """
    fig, ax = plt.subplots()
    fig.suptitle(title)
    ax.set_xlabel('Log X (cm)')
    ax.set_ylabel('Log Z (cm)')

    size = (len(grid_x)-1, len(grid_z)-1)
    data = np.reshape(table['var'], size)
    if volume:
        for xi in range(size[0]):
            for zi in range(size[1]):
                # We need to correct the per-area emission to a per-volume emission
                area = 2. * 2. * np.pi * (grid_x[xi+1]-grid_x[xi]) * (grid_z[zi+1]-grid_z[zi])
                area += 2. * np.pi * grid_x[xi+1] * (grid_z[zi+1]-grid_z[zi])
                area += 2. * np.pi * grid_x[xi]   * (grid_z[zi+1]-grid_z[zi])
                if area > 0:
                    data[xi, zi] = data[xi, zi] / area

    im = ax.pcolormesh(np.log10(grid_x), np.log10(grid_z), np.ma.log10(data.T))
    ax.set_xlim(14.5, 18)
    ax.set_ylim(13, 17)
    cbar = fig.colorbar(im, ax=ax).set_label(label)
    return fig


def plot_dat_many(tables, grids_x, grids_z, xlims, zlims, titles, title, label,
                  shared_y=False, shared_cbar=False, volume=True, log=True):
    """
    Plot many dat files on a single plot.

    Use:
        plot_dat_many([table_h1_r01, table_h1_r10, table_h1_r30],
                      [x_r01, x_r10, x_r30], [z_r01, z_r10, z_r30],
                      xlims=[(14.5, 17.5), (15.5, 17.5), (16, 17.5)],
                      zlims=[(13, 17), (13, 17), (13, 17)],
                      titles=['1x Radius', '10x Radius', '30x Radius'],
                      title='H-I ion fraction', label='Log ion fraction',
                      shared_y=True, volume=False, shared_cbar=True)
    """
    if shared_y:
        fig, axes = plt.subplots(1, len(tables), sharey='row')
    else:
        fig, axes = plt.subplots(1, len(tables))

    fig.suptitle(title, y=1.03)
    axes[0].set_ylabel('Log Z (cm)')
    images = []

    if shared_cbar:
        vmin = float('inf')
        vmax = -float('inf')

    for ax, table, grid_x, grid_z, xlim, zlim, title in zip(axes, tables, grids_x, grids_z, xlims, zlims, titles):
        ax.set_title(title)
        ax.set_xlabel('Log X (cm)')

        size = (len(grid_x)-1, len(grid_z)-1)
        data = np.reshape(table['var'], size)
        inwind = np.reshape(table['inwind'], size)

        if volume:
            for xi in range(size[0]):
                for zi in range(size[1]):
                    # We need to correct the per-area emission to a per-volume emission
                    area = 2. * 2. * np.pi * (grid_x[xi+1]-grid_x[xi]) * (grid_z[zi+1]-grid_z[zi])
                    area += 2. * np.pi * grid_x[xi+1] * (grid_z[zi+1]-grid_z[zi])
                    area += 2. * np.pi * grid_x[xi]   * (grid_z[zi+1]-grid_z[zi])
                    if area > 0:
                        data[xi, zi] = data[xi, zi] / area

        if shared_cbar:
            if data.max() > vmax:
                vmax = data.max()
            if data.min() < vmin:
                vmin = data.min()

        if log:
            image = ax.pcolormesh(np.log10(grid_x), np.log10(grid_z),
                                  np.ma.log10(data.T))
        else:
            image = ax.pcolormesh(np.log10(grid_x), np.log10(grid_z),
                                  np.ma.masked_where(inwind.T < 0, data.T))

        images.append(image)
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(zlim[0], zlim[1])

    if not shared_cbar:
        for image, ax in zip(images, axes):
            cbar = fig.colorbar(image, ax=ax).set_label(label)
    else:
        for image, ax in zip(images, axes):
            image.clim = [vmin, vmax]

        cbar = fig.colorbar(images[-1], ax=axes[-1]).set_label(label)

    fig.tight_layout()
    if shared_y:
        fig.subplots_adjust(wspace=0)
    return fig
