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
import matplotlib.pyplot as plt

sys.path.append(os.path.join('..', '..'))
import blobby3d_toolkit as b3dtk

import dnest4 as dn4

dn4.postprocess()

b3d = b3dtk.Blobby3D(
        samples_path='posterior_sample.txt',
        data_path='data.txt',
        var_path='var.txt',
        metadata_path='metadata.txt',
        nlines=2)

# choose a sample
sample = 0

# Plot maps for sample
fig, ax = plt.subplots(1, 4)
map_names = ['FLUåX0', 'FLUX1', 'V', 'VDISP']
for i in range(4):
    ax[i].set_title(map_names[i])
    ax[i].imshow(b3d.maps[sample, i], interpolation='nearest', origin='lower')
fig.tight_layout()

# We can also plot the integrated flux across the wavelength axis for sample
# and compare it to the data
fig, ax = plt.subplots(1, 3)
ax[0].set_title('Preconvolved')
ax[0].imshow(
        b3d.precon_cubes[sample].sum(axis=2),
        interpolation='nearest', origin='lower')

ax[1].set_title('Convolved')
ax[1].imshow(
        b3d.con_cubes[sample].sum(axis=2),
        interpolation='nearest', origin='lower')

ax[2].set_title('Data')
ax[2].imshow(
        b3d.data.sum(axis=2),
        interpolation='nearest', origin='lower')
fig.tight_layout()

# Marginalised samples of flux for each line for sample
b3d.blob_param.loc[sample, ['FLUX0', 'FLUX1']].hist(bins=30)

# Marginalised samples for the mean and standard deviation of the flux for the
# first emission line
b3d.global_param[['FLUX0MU', 'FLUX0SD']].hist()
