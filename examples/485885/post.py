#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Postprocess for the example galaxy.

The PostBlobby3D class is very simple. However, it is useful for organisational
purposes of the Blobby3D output. In this script, I created a PostBlobby3D
object and plotted a handful of sample attributes.

@author: Mathew Varidel

"""

import numpy as np
import matplotlib.pyplot as plt

import dnest4 as dn4
from pyblobby3d import PostBlobby3D

dn4.postprocess()

post_b3d = PostBlobby3D(
        samples_path='posterior_sample.txt',
        data_path='data.txt',
        var_path='var.txt',
        metadata_path='metadata.txt',
        nlines=2)

# choose a sample
sample = 0

# Plot maps for sample
fig, ax = plt.subplots(1, 4)

ax[0].set_title(r'H$\alpha$ Flux')
ax[0].imshow(
    np.log10(post_b3d.maps[sample, 0]),
    interpolation='nearest', origin='lower')

ax[1].set_title(r'[NII] Flux')
ax[1].imshow(
    np.log10(post_b3d.maps[sample, 1]),
    interpolation='nearest', origin='lower')

ax[2].set_title('V')
ax[2].imshow(post_b3d.maps[sample, 2], interpolation='nearest', origin='lower')

ax[3].set_title('V Disp')
ax[3].imshow(post_b3d.maps[sample, 3], interpolation='nearest', origin='lower')

fig.tight_layout()

# We can also plot the integrated flux across the wavelength axis for a sample
# and compare it to the data. The below does this for H-alpha.
fig, ax = plt.subplots(1, 3)

ax[0].set_title('Preconvolved')
ax[0].imshow(
        np.log10(post_b3d.precon_cubes[sample].sum(axis=2)),
        interpolation='nearest', origin='lower')

ax[1].set_title('Convolved')
ax[1].imshow(
        np.log10(post_b3d.con_cubes[sample].sum(axis=2)),
        interpolation='nearest', origin='lower')

ax[2].set_title('Data')
ax[2].imshow(
        np.log10(post_b3d.data.sum(axis=2)),
        interpolation='nearest', origin='lower')

fig.tight_layout()

# Similarly we can integrate the total cube and look at the flux as
# a function of wavelength. In this case lets compare all convolved samples
# to the data.
fig, ax = plt.subplots()

ax.plot(post_b3d.data.sum(axis=(0, 1)), '--k', label='Data')
ax.plot(post_b3d.con_cubes.sum(axis=(1, 2)).transpose(), '--', color='0.5')

ax.legend()
