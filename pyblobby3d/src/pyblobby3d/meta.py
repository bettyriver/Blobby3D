#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Metadata.

Organisation of the Blobby3D metadata that describes the coordinates.

@author: Mathew Varidel
"""

from pathlib import Path
import numpy as np


class Metadata:

    def __init__(self, metadata_path):
        """Metadata oranisation object.

        Parameters
        ----------
        metadata_path : str or pathlib.Path
            DESCRIPTION.

        Returns
        -------
        None.

        """
        metadata = np.loadtxt(Path(metadata_path))
        self.naxis = metadata[:3].astype(int)
        self.sz = self.naxis.prod()
        self.x_lim = metadata[3:5]
        self.y_lim = metadata[5:7]
        self.r_lim = metadata[7:9]
        self.dx = float(np.diff(self.x_lim)/self.naxis[1])
        self.dy = float(np.diff(self.y_lim)/self.naxis[0])
        self.dr = float(np.diff(self.r_lim)/self.naxis[2])
