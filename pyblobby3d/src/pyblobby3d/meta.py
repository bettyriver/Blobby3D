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

        Returns
        -------
        None.

        Attributes
        ----------
        sz : float
        x_lim : (float, float)
            Left and right most boundaries of the x axis.
        y_lim : float
            Bottom and top most boundaries of the y axis.
        r_lim : float
            Left and right most boundaries of the wavelength axis.
        dx : float
            Width of pixels along the x axis.
        dy : float
            Width of pixels along the y axis.
        dr : float
            Width of pixels along the wavelength axis.

        """
        metadata = np.loadtxt(Path(metadata_path))
        self.naxis = metadata[:3].astype(int)
        self.sz = self.naxis.prod()
        self.x_lim = metadata[3:5]
        self.y_lim = metadata[5:7]
        self.r_lim = metadata[7:9]
        self.dx = np.diff(self.x_lim)[0]/self.naxis[1]
        self.dy = np.diff(self.y_lim)[0]/self.naxis[0]
        self.dr = np.diff(self.r_lim)[0]/self.naxis[2]

    def get_axis_array(self, axis):
        """Get axis array.

        Calculate the x, y or wavelength array.

        Parameters
        ----------
        axis : str
            Axis defind by 'x', 'y', 'r'. Where 'r' is the wavelength.

        Returns
        -------
        axis_array : np.ndarray

        """
        if axis == 'x':
            lim = self.x_lim
            d = self.dx
            n = self.naxis[0]
        elif axis == 'y':
            lim = self.y_lim
            d = self.dy
            n = self.naxis[1]
        elif axis == 'r':
            lim = self.r_lim
            d = self.dr
            n = self.naxis[2]

        return np.linspace(lim[0] + 0.5*d, lim[1] - 0.5*d, n)

