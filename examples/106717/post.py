#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Postprocess for the example galaxy.

Note that blobby3d.py is required to be in your path -- I append the relative
path using the sys module in this example. Feel free to solve it in your own
way.

The Blobby3D class is very simple (some would say it's a bad use of classes!).
However, it is useful for organisational purposes of the Blobby3D output. In
this script, I created a Blobby3D object and plotted a handful of sample
attributes.

@author: Mathew Varidel

"""

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join('..', '..', 'python'))
from blobby3d import Blobby3D

from moments import SpectralModel
import b3dplot
import b3dcomp


def plot_maps(sample, figsize=(9.0, 7.0)):
    fig, ax = plt.subplots(1, 4, figsize=figsize)
    b3d.plot_map(
            ax[0], np.ma.masked_invalid(np.log10(b3d.maps[sample, 0, :, :])),
            title=r'log$_{10}$(H$\alpha$)', colorbar=True, cmap='Greys_r')
    b3d.plot_map(
            ax[1], np.ma.masked_invalid(np.log10(b3d.maps[sample, 1, :, :])),
            title=r'log$_{10}$([NII])', colorbar=True, cmap='Greys_r')
    b3d.plot_map(
            ax[2], b3d.maps[sample, 2, :, :],
            title='$v$', colorbar=True, cmap='RdYlBu_r')
    b3d.plot_map(
            ax[3], b3d.maps[sample, 3, :, :],
            title=r'$\sigma_v$', colorbar=True, cmap='YlOrBr')
    fig.tight_layout()

    return fig, ax


def plot_n2_ha(sample, figsize=(9.0, 7.0),
               min_logflux=-np.inf, only_starforming=True):
    fig, ax = plt.subplots(1, 3, figsize=figsize)

    mask = np.log10(b3d.maps[sample, 0, :, :]) < min_logflux
    mask |= b3d.var.sum(axis=2) <= 0.0

    if only_starforming:
        mask |= b3d.maps[sample, 1, :, :]/b3d.maps[sample, 0, :, :] > 1.0

    mask |= ~np.isfinite(mask)

    b3d.plot_map(
            ax[0], np.ma.masked_where(mask, np.log10(b3d.maps[sample, 0, :, :])),
            title=r'log$_{10}$(H$\alpha$)', colorbar=True, cmap='Greys_r')
    b3d.plot_map(
            ax[1], np.ma.masked_where(mask, np.log10(b3d.maps[sample, 1, :, :])),
            title=r'log$_{10}$([NII])', colorbar=True, cmap='Greys_r')

    n2_ha = b3d.maps[sample, 1, :, :]/b3d.maps[sample, 0, :, :]
    n2_ha = np.ma.masked_where(mask, np.log10(n2_ha))
    b3d.plot_map(
            ax[2], n2_ha,
            title=r'log$_{10}$([NII]/H$\alpha$)',
            colorbar=True, cmap='RdYlBu_r',
            clim=(-1.0, 1.0))
    fig.tight_layout()

    return fig, ax


def plot_cube(sample, path=r'cube.pdf'):
    pdf = mpl.backends.backend_pdf.PdfPages(path)
    fig, ax = plt.subplots()
    for i in list(range(b3d.naxis[0]))[:20]:
        for j in list(range(b3d.naxis[1]))[:20]:
            if b3d.data[i, j, :].sum() > 0.0:
                for l in range(len(ax.lines)):
                    ax.lines.pop()
                ax.plot(b3d.precon_cubes[sample, i, j, :],
                        'k', label='Preconvolved')
                ax.plot(b3d.con_cubes[sample, i, j, :],
                        'r', label='Convolved')
                ax.plot(b3d.data[i, j, :],
                        color='0.5', label='Data')
                fig.legend()
                pdf.savefig(fig)
    plt.close()
    pdf.close()


if __name__ == '__main__':
    b3d = Blobby3D(
            samples_path='posterior_sample.txt',
            data_path='data.txt',
            var_path='var.txt',
            metadata_path='metadata.txt',
            nlines=2)

    # choose a sample
    sample = 0

    # Construct comparison plot
    sm = SpectralModel(
            lines=[[6562.81], [6583.1, 6548.1, 0.3333]],
            lsf_fwhm=1.59
            )

    wave = np.linspace(
            b3d.r_lim[0] + 0.5*b3d.dr,
            b3d.r_lim[1] - 0.5*b3d.dr,
            b3d.naxis[2])

    fit, fit_err = sm.fit_cube(
            wave, b3d.con_cubes[sample, :, :, :])

    fit_data, fit_data_err = sm.fit_cube(
            wave, b3d.data)

    fig, ax = b3d.setup_comparison_maps()
    b3d.add_comparison_maps(b3d.maps[sample, :, :], ax, col=0, log_flux=True)
    b3d.add_comparison_maps(fit, ax, col=1, log_flux=True)
    b3d.add_comparison_maps(fit_data, ax, col=2, log_flux=True)

    residuals = fit - fit_data
    residuals[0, :, :] = np.log10(fit[0, :, :]) - np.log10(fit_data[0, :, :])
    residuals[1, :, :] = np.log10(fit[0, :, :]) - np.log10(fit_data[0, :, :])
    b3d.add_comparison_residuals(residuals, ax, col=4)

    b3d.update_comparison_clim(ax[0][:3], ax[0][3], pct=100.0, absolute=False)
    b3d.update_comparison_clim(ax[1][:3], ax[1][3], pct=100.0, absolute=False)
    b3d.update_comparison_clim(ax[2][:3], ax[2][3], pct=100.0, absolute=True)
    b3d.update_comparison_clim(ax[3][:3], ax[3][3], pct=100.0, absolute=False)

    b3d.update_comparison_clim(ax[0][4], ax[0][5], pct=100.0, absolute=True)
    b3d.update_comparison_clim(ax[1][4], ax[1][5], pct=100.0, absolute=True)
    b3d.update_comparison_clim(ax[2][4], ax[2][5], pct=100.0, absolute=True)
    b3d.update_comparison_clim(ax[3][4], ax[3][5], pct=100.0, absolute=True)
