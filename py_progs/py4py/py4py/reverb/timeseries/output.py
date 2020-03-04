"""
Timeseries output module

This module contains the output code for the timeseries analysis. This is intended to output to the formats
the CARAMEL and MEMEcho team analyse & plot.
"""
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as ani
from matplotlib.figure import figaspect

from astropy import units as u
from astropy.table import Table
from astropy.units import Quantity

import shutil

from py4py.physics import doppler_shift_wave
from py4py.reverb import TransferFunction

from typing import List, Optional, Tuple


def write_caramel_data(lightcurve: Table, spectra: Table, spectra_times: Table, suffix: str):
    """
    Given a lightcurve, series of spectra and time outputs to CARAMEL format, and then
    compresses the output into a ZIP file.

    Arguments:
        lightcurve (Table): Continuum values and times, in seconds and real units.
        spectra (Table): Table of wavelengths and spectra. Continuum-subtracted.
        spectra_times (Table): Table of spectrum times.
        suffix (str): Suffix appended to filename. Intended to sort outputs as e.g. caramel/qso/caramel_lightcurve_qso.

    Outputs:
        caramel/{suffix}/caramel_lightcurve_{suffix}.txt: Continuum lightcurve.
        caramel/{suffix}/caramel_spectra_{suffix}.txt: Values of the mock observations.
        caramel/{suffix}/caramel_spectra_times_{suffix}.txt: Times of the mock observations.
        caramel/caramel_{suffix}.zip: Zip of all the outputs.
    """
    # Lightcurve file format:
    # Time (s)
    # Value (rescaled from 1-100)
    # Error (rescaled as above)

    time = np.array(lightcurve['time'].quantity.to(u.s))
    # Rescale value and error from 1-100
    value = np.array((lightcurve['value'] - np.amin(lightcurve['value']))
                     / (np.amax(lightcurve['value']) - np.amin(lightcurve['value'])))
    value = (value * 9) + 1
    error = np.array(lightcurve['error']
                     / (np.amax(lightcurve['value']) - np.amin(lightcurve['value'])))
    error = (error * 9)

    np.savetxt(
        os.path.join(
            'caramel', suffix, 'caramel_lightcurve_{}.txt'.format(suffix)
        ),
        np.column_stack((time, value, error))
    )

    # Spectrum times file format:
    # Time (s)
    # Dummy value (0)
    # Dummy value (0)
    times = np.array(spectra_times['time'].to(u.s))
    dummy = [0] * len(times)

    np.savetxt(
        os.path.join(
            'caramel', suffix, 'caramel_spectra_times_{}.txt'.format(suffix)
        ),
        np.column_stack((times, dummy, dummy))
     )

    # Spectra file format:
    # First row is # INT, where INT is number of wavelength bins
    # Second row is central pixel wavelength in angstroms
    # Third row is rescaled flux from 1-100 for spectrum 1
    # Fourth row is rescaled flux error for spectrum 1
    to_save = [spectra['wave']]
    spec_min = 999e99
    spec_max = -999e99
    for column in spectra.colnames[1:]:
        if np.amax(spectra[column]) > spec_max:
            spec_max = np.amax(spectra[column])
        if np.amin(spectra[column]) < spec_min:
            spec_min = np.amin(spectra[column])

    for column in spectra.colnames[5:]:
        value = (np.array(spectra[column]) - spec_min) / (spec_max - spec_min)
        value = (value * 9) + 1
        error = np.array(spectra['error']) / (spec_max - spec_min)
        error = (error * 9)
        to_save.append(value)
        to_save.append(error)

    np.savetxt(
        os.path.join(
            'caramel', suffix, 'caramel_spectra_{}.txt'.format(suffix)
        ),
        to_save, header='{}'.format(len(spectra))
    )
    shutil.make_archive(
        'caramel_{}'.format(suffix), 'zip',
        os.path.join('caramel', suffix)
    )
    return


def write_memecho_data(lightcurve: Table, spectra: Table, spectra_times: Table, suffix: str):
    """
    Given a lightcurve, series of spectra and time outputs to MEMECHO format, outputs them
    to files and zips them up for distribution.

    Arguments:
        lightcurve (Table): Continuum values and times, in seconds and real units.
        spectra (Table): Table of wavelengths and spectra. Not continuum-subtracted.
        spectra_times (Table): Table of spectrum times.
        suffix (str): Suffix appended to file name.

    Outputs:
        memecho/{suffix}/prepspec_{suffix}_{timestep}.txt: The mock observations for each timestep
        memecho/{suffix}/prepspec_times.txt: The times of each mock observation, and spectra file associated with each.
        memecho/{suffix}/memecho_lightcurve.txt: The driving continuum lightcurve.
        memecho/memecho_{suffix}.zip: Zip of all the outputs
    """
    names = []
    for index, column in enumerate(spectra.colnames[5:]):
        np.savetxt(
            os.path.join(
                'memecho', suffix, 'prepspec_{}_{:03d}.txt'.format(suffix, index)
            ),
            np.column_stack((spectra['wave'], spectra[column], spectra['error']))
        )
        names.append(
            os.path.join('memecho', suffix, 'prepspec_{}_{:03d}.txt'.format(suffix, index))
        )

    with open(os.path.join('memecho', suffix, 'prepspec_times_{}.txt'.format(suffix)), 'w') as file:
        for name, time in zip(names, spectra_times['time']):
            file.write("{} {}\n".format(name, time))

    with open(os.path.join('memecho', suffix, 'memecho_lightcurve_{}.txt'.format(suffix)), 'w') as file:
        for time, value, error in zip(lightcurve['time'], lightcurve['value'], lightcurve['error']):
            file.write("{} {} {}\n".format(time, value, error))

    shutil.make_archive(
        'memecho_{}'.format(suffix), 'zip',
        os.path.join('memecho', suffix)
    )


def trailed_spectrogram(spectra: Table, lightcurve: Table, spectra_times: Table, filename: str):
    """
    Generate a trailed spectrogram of both the time series of spectra and difference relative to the mean,
    with the continuum as an adjacent line plot.

    Arguments:
        spectra (Table): Spectra (starting at column 3).
        spectra_times (Table): Times to plot the TS for.
        lightcurve (Table): The continuum lightcurve.
        filename (String): File to write to.

    Outputs:
        {filename}.eps: Time series output.
    """

    # We want a pcolour plot for the time series, with an adjacent
    # fig, (ax_ts, ax_c) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]} , sharey=True)

    w, h = figaspect(1.25/1)
    fig, ((ax_spec, ax_none, ax_none2), (ax_ts, ax_c, ax_cb), (ax_ts2, ax_c2, ax_cb2)) = \
        plt.subplots(
            3, 3, sharex='col', sharey='row',
            gridspec_kw={'width_ratios': [3, 1, 1], 'height_ratios': [1, 3, 3]}, figsize=(w, h)
        )

    ax_none.axis('off')
    ax_none2.axis('off')
    ax_cb.axis('off')
    ax_cb2.axis('off')

    pcolor_data = np.zeros([len(spectra_times), len(spectra)])
    for i, column in enumerate(spectra.colnames[5:]):
        pcolor_data[i, :] = spectra[column]

    # We want to set the upper and lower bounds for each timestep. Use the midpoint between steps,
    # and assume start and end timesteps symmetric.
    time_bounds = np.zeros(len(spectra_times)+1)
    for i in range(0, len(spectra_times)-1):
        time_bounds[i+1] = (spectra_times['time'][i] + spectra_times['time'][i+1]) / 2
    time_bounds[0] = (spectra_times['time'][0] - (spectra_times['time'][1] - spectra_times['time'][0]) / 2)
    time_bounds[-1] = (spectra_times['time'][-1] + (spectra_times['time'][-1] - spectra_times['time'][-2]) / 2)

    # Now we do the pcolour plot, with the *true* lightcurve along the side. Maybe we truncate or remove this...
    pcol = ax_ts.pcolor(
        spectra.meta["bounds"], time_bounds, pcolor_data/np.amax(pcolor_data)
    )
    delta = (pcolor_data - spectra['value']) / np.amax(pcolor_data)
    cb_max = np.amax(np.abs(delta))

    pcol2 = ax_ts2.pcolor(
        spectra.meta["bounds"], time_bounds, delta, vmin=-cb_max, vmax=cb_max, cmap='RdBu_r'
    )

    for ax in (ax_ts, ax_ts2):
        ax.set_xlabel(r'Wavelength ($\AA$)')
        ax.set_ylabel("Time (MJD)")
        ax.set_ylim(time_bounds[0], time_bounds[-1])
        ax.xaxis.set_tick_params(rotation=0, pad=1)
        ax.yaxis.set_tick_params(rotation=45, labelsize=8)
        ax.yaxis.tick_left()
        ax.set_xlim([6300, 6850])

    ax_spec.set_xlim([6300, 6850])
    ax_spec.set_xlabel("λ (Å)")
    ax_spec.set_ylabel(r'$L/L_{\rm max}$')
    ax_spec.plot(spectra['wave'], spectra['value']/np.amax(spectra['value']))

    for ax in (ax_c, ax_c2):
        ax.invert_xaxis()
        ax.plot(
            100.0*(lightcurve['value']-np.mean(lightcurve['value']))/np.mean(lightcurve['value']),
            lightcurve['time'], '-', c='m'
        )
        ax.xaxis.set_tick_params(rotation=0, pad=12)

    ax_c2.set_xlabel(r"ΔC (%)")

    tf_wave = 6562.8
    ax_ts.axvline(tf_wave, color='red')
    ax_ts2.axvline(tf_wave, color='red')
    ax_spec.axvline(tf_wave, color='red')

    ax_vel = ax_spec.twiny()
    ax_vel.set_xlim(ax_spec.get_xlim())
    ax_vel.set_xticks([
        doppler_shift_wave(tf_wave, -1e7),
        tf_wave,
        doppler_shift_wave(tf_wave, +1e7)
    ])
    ax_vel.set_xticklabels([
        r"10", r"0", r"10"
    ])
    ax_vel.set_xlabel(r'Velocity (10$^{3}$ km s$^{-1}$)')

    cbar1 = fig.colorbar(pcol, ax=ax_cb, orientation="vertical", fraction=1)
    cbar1.set_label(r"$L/L_{max}$")
    cbar1.ax.tick_params(labelsize=8)
    cbar2 = fig.colorbar(pcol2, ax=ax_cb2, orientation="vertical", fraction=1)
    cbar2.set_label(r'$\Delta L/L_{\rm max}$')
    cbar2.ax.tick_params(labelsize=8)
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig("{}.eps".format(filename), bbox_inches='tight')
    plt.close(fig)


next_column = None  # Global variable required for plotting write_animation with matplotlib


def write_animation(
        spectra: Table, lightcurve: Table, spectra_times: Table, times: Table, filename: str, is_reversed: bool = False
):
    """
    Given a lightcurve and table containing a time series of spectra,
    generate an animation that shows how the output spectrum changes over time.

    Arguments:
        spectra (Table): Spectra (starting at column 3).
        spectra_times (Table): Times to plot the spectra values for
        times (Table): High time-resolution interpolated continuum lightcurve.
        lightcurve (Table): The continuum lightcurve.
        filename (str): File to write to
        is_reversed (bool): Whether newer points should be overlaid by older ones

    Outputs:
        {filename}.mp4
    """
    global next_column
    figure, (axis, ax_dc) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    figure.subplots_adjust(hspace=.3, wspace=0)

    # Find the maxima and minima over all the spectra to scale the plot appropriately
    ymin = +np.inf
    ymax = -np.inf
    for column in spectra.colnames[1:]:
        if np.amax(spectra[column]) > ymax:
            ymax = np.amax(spectra[column])
        if np.amin(spectra[column]) < ymin:
            ymin = np.amin(spectra[column])

    # For the upper (spectra) axis, label and set limits
    axis.set_xlabel("λ ({})".format(spectra['wave'].meta['name']))
    axis.set_ylabel("L ({})".format(spectra['value'].meta['name']))
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([spectra['wave'][0], spectra['wave'][-1]])

    # For the lower (continuum) axis, label and set limits
    y_dc_min = np.amin(times['dC'])
    y_dc_max = np.amax(times['dC'])
    ax_dc.set_xlabel("T ({})".format(times['time'].meta['name']))
    ax_dc.set_ylabel("ΔC ({})".format(times['dC'].meta['name']))
    ax_dc.set_xlim([lightcurve['time'][0], lightcurve['time'][-1]])
    ax_dc.set_ylim([y_dc_min, y_dc_max])
    # Put a line for the mean continuum on time-continuum
    line_mean = ax_dc.plot((times['time'][0], times['time'][-1]), (0.0, 0.0), '-', c='k', zorder=-len(times))
    # Set up a second axis to show % difference
    ax_dc_percent = ax_dc.twinx()
    ax_dc_percent.set_ylim([y_dc_min * 100.0 / lightcurve.meta['mean'].value,
                            y_dc_max * 100.0 / lightcurve.meta['mean'].value])
    ax_dc_percent.set_ylabel("ΔC (%)")
    ax_dc_percent.get_yaxis().get_major_formatter().set_useOffset(False)

    # Plot the 'original' spectrum in black
    line_start = axis.plot(spectra['wave'], spectra['value'], '-', c='k', zorder=0, label='Baseline')
    line_min = axis.plot(spectra['wave'], spectra['value_min'], '--', c='k', zorder=0, label='Minimum')
    line_max = axis.plot(spectra['wave'], spectra['value_max'], '-.', c='k', zorder=0, label='Maximum')

    # Set up colour palette with time
    # line_colours = [plt.cm.jet(x) for x in np.linspace(0.0, 1.0, len(spectra_times))]

    # Artists is the list of objects redrawn each write_animation (so all plots)
    artists = [line_start]
    # We evaluate the spectra one-by-one for timestamps, starting with 1st (i.e. spectra column 2)
    next_column = 5

    def zorder(step: int, steps: int, reverse: bool = False) -> int:
        """
        Function for returning the z-order for each observation, used when rendering the plot.
        Do we want new points over old, or old points over new?

        Arguments:
            step (int): Current timestep.
            steps (int): Number of timesteps.
            reverse (bool): Whether newer points should get overlaid by older ones.

        Returns:
            int: If reverse, a number that goes down with timestep.
                If not, a number that goes up with timestep (but remains below 0).
        """
        # Straightforward- do we want new points over old, or old over new?
        # Latter is useful for showing how a step change propagates
        if reverse:
            return -step
        else:
            return step-steps

    def update_figure(step: int):
        """
        Called each timestep to update the plot object.

        Arguments:
            step (int): The current timestep.
        """
        global next_column

        time_colour = 'grey'
        if spectra_times['time'][0] <= times['time'][step] <= spectra_times['time'][-1]:
            # If this timestep is in the range of spectra, assign it a colour
            time_colour = plt.cm.jet(
                (times['time'][step]-spectra_times['time'][0]) / (spectra_times['time'][-1] - spectra_times['time'][0])
            )

        if times['time'][step] > spectra.columns[next_column].meta['time']:
            # If we've gone past the current column,
            line_new = axis.plot(spectra['wave'], spectra.columns[next_column], '-', markersize=0.5,
                                 c=time_colour, zorder=zorder(step, len(times), is_reversed))
            times_for_line = np.array([times['time'][step], times['time'][step]])
            values_for_line = np.array([y_dc_min, y_dc_max])
            line_vertical = ax_dc.plot(times_for_line, values_for_line, ':', linewidth=1, c=time_colour, zorder=0)
            artists.append(line_vertical)
            artists.append(line_new)
            next_column += 1

        # Plot the continuum brightness for this step in the appropriate colour
        point_new = ax_dc.plot(
            times['time'][step], times['dC'][step], 'o', c=time_colour, zorder=zorder(step, len(times), is_reversed)
        )
        artists.append(point_new)
        return artists,

    # Generate write_animation and save to file
    animation = ani.FuncAnimation(figure, update_figure, frames=len(times))
    animation.save("{}.mp4".format(filename), fps=24)
    plt.close(figure)


def rescaled_rfs(
        tfs: List[TransferFunction], rescale_max_time: Quantity, figure_max_time: Quantity, keplerian: dict = None
):
    """
    Outputs response functions for rescaled versions of the input.
    Different mass SMBHs scale straightforwardly; accretion disks are generated at locations with the same
    doppler shifts, so the only thing you need to do to scale the mass is to rescale the time delays.

    Arguments:
        tfs (List[TransferFunction]): The transfer functions to plot the response functions for.
        rescale_max_time (Quantity): The new 'maximum time' the TF should extend to.
        figure_max_time (Quantity):
            The maximum time that should be shown on the plot (may be lower than rescale_max_time).
        keplerian (dict): Dictionary containing Keplerian disk profile as used by plot.

    Outputs:
        {tf.name}_resp.eps
    """
    for tf in tfs:
        delay_bins = (tf.delay_bins() * u.s).to(u.day)
        print('rescaled_rfs: Rescale factor is:', (rescale_max_time / delay_bins[-1]).value)
        tf._bins_delay *= (rescale_max_time / delay_bins[-1]).value
        keplerian['rescale'] = (rescale_max_time / delay_bins[-1]).value
        print('rescaled_rfs: Rescale raw peak is:', tf.delay_peak(response=True, days=True))
        print('rescaled_rfs: Rescale centroid is:', tf.delay(days=True))
        tf.plot(response_map=True, name='resp', max_delay=figure_max_time,
                keplerian=keplerian)


def plot_spectra_rms(spectra: List[(Table, Table)], filenames: List[str]):
    """
    Given a list of timeseries of spectra (full and continuum subtracted), produce a trailed
    spectrogram of each, plus the RMS spectra.

    Arguments:
        spectra (List[(Table, Table)]):
            Pairs of tables containing full and continuum subtracted (in that order) timeseries of spectra
        filenames (List[str]): Filenems
            List of filenames for each of the pairs.

    Outputs:
        {filename}.eps for filename in filenames
    """
    for (spec_full, spec_line), filename in zip(spectra, filenames):
        fig, (ax_full, ax_line) = plt.subplots(2, 1, sharex=True)
        plt.tight_layout()
        fig.subplots_adjust(hspace=0, wspace=0)
        fig.legend()

        ax_line.set_ylabel("Spectrum ({})".format(spec_line['value'].meta['name']))
        ax_line.set_xlabel("λ ({})".format(spec_line['wave'].meta['name']))

        full_series = np.array([spec_full[column] for column in spec_full.columns[5:]])
        line_series = np.array([spec_line[column] for column in spec_line.columns[5:]])

        rms_full = np.sqrt(np.sum(np.square(full_series), 1) / len(spec_full))
        rms_line = np.sqrt(np.sum(np.square(line_series), 1) / len(spec_line))

        ax_full.errorbar(
            spec_full['wave'], spec_full['value'], fmt='-', c='r', yerr=spec_full['error'], label='Spectrum'
        )
        ax_full.errorbar(
            spec_full['wave'], rms_full, fmt='-', c='b', label='RMS'
        )
        ax_line.errorbar(
            spec_line['wave'], spec_line['value'], fmt='-', c='r', yerr=spec_line['error'], label='Spectrum'
        )
        ax_line.errorbar(
            spec_line['wave'], rms_line, fmt='-', c='b', label='RMS'
        )

        plt.savefig("{}.eps".format(filename), bbox_inches='tight')
        plt.close(fig)
