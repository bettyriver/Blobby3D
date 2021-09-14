"""Fitting kinematic moments.

@author: mathewvaridel
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

from .const import PhysicalConstants


class SpectralModel:

    def __init__(self, lines, lsf_fwhm, baseline_order=None, wave_ref=0.0):
        """Spectral model.

        Parameters
        ----------
        lines : list of lists
            DESCRIPTION.
        lsf_fwhm : float
            Full-Width Half Maximum of Line Spread Function.
        baseline_order : int, optional
            Order of baseline polynomial. The default is None.
        wave_ref : float, optional
            Redshift offset in Angstroms. The default is 0.0.

        Returns
        -------
        None.

        """
        self.lines = lines
        self.nlines = len(lines)

        self.lsf_fwhm = lsf_fwhm
        self.lsf_sigma = lsf_fwhm/np.sqrt(8.0*np.log(2.0))

        self.baseline_order = baseline_order

        self.wave_ref = wave_ref

        self.nparam = self.nlines + 2
        if self.baseline_order is not None:
            self.nparam += 1 + self.baseline_order

    def calculate(self, wave, *param):
        """Calculate spectral model.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        param : arra-like
            Model parameters. The format is
            [flux_1, flux_2, .., velocity, velocity dispersion, b0, b1, ...].
            Where flux_i is the flux for the i-th coupled emission line and
            b0, b1, ... are the baseline polynomial coefficients.

        Returns
        -------
        model : np.ndarray

        """
        if self.baseline_order is None:
            model = self._gas_model(wave, param)
        else:
            model = self._gas_model(
                    wave,
                    param[:self.nlines+2])
            model += self._baseline_model(
                    wave,
                    param[-1-self.baseline_order:],
                    )

        return model

    def calculate_cube(self, wave, param):
        """Calculate spectral model for cube.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        param : 3D np.ndarray
            Cube of model parameters. For each spaxel the parameter format is
            [flux_1, flux_2, .., velocity, velocity dispersion]. Where flux_i
            is the flux for the i-th coupled emission line.
        wave_axis : int
            Axis for wavelength can either be 0 or 2. Default is 2.

        Returns
        -------
        model : np.ndarray
            Model cube with axis (i, j, wavelength).

        """
        model_cube = np.zeros((*param.shape[1:], len(wave)))*np.nan
        for i in range(param.shape[1]):
            for j in range(param.shape[2]):
                model_cube[i, j] = self.calculate(wave, *param[:, i, j])

        return model_cube

    def fit_spaxel(self, wave, data, var=None, bounds=None):
        """Fit spaxel using curve_fit.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        data : array-like
        var : array-like, default is None
            Variance array. Default is None, which means curve_fit is not
            supplied a variance array.
        bounds : list of lists
            Bounds supplied to curve_fit for each parameter.
         wave_axis : int
            Axis for wavelength can either be 0 or 2. Default is 2.

        Returns
        -------
        fit : np.ndarray
            Fit parameters returned by curve_fit.
        fit_err : np.ndarray
            Fit errors returned by curve_fit.

        """
        if bounds is None:
            # line flux bounds
            bounds = [
                    [0.0]*self.nlines,
                    [np.inf]*self.nlines,
                    ]

            # v, vdisp bounds
            bounds[0] += [-np.inf, 0.0]
            bounds[1] += [np.inf, np.inf]

            if self.baseline_order is not None:
                bounds[0] += [-np.inf]*(self.baseline_order + 1)
                bounds[1] += [np.inf]*(self.baseline_order + 1)

        if (var is None) & np.any(data != 0.0):
            data_valid = np.isfinite(data)
            data_tmp = data[data_valid]
            w_tmp = wave[data_valid]
            sigma_tmp = None
        elif (var is None) & np.all(data == 0.0):
            popt = np.zeros(self.nparam)*np.nan
            pcov = np.zeros(self.nparam)*np.nan
            return popt, pcov
        elif np.any(var > 0.0):
            data_valid = ((var > 0.0) & np.isfinite(var) & np.isfinite(data))
            data_tmp = data[data_valid]
            w_tmp = wave[data_valid]
            sigma_tmp = np.sqrt(var[data_valid])
        else:
            popt = np.zeros(self.nparam)*np.nan
            pcov = np.zeros(self.nparam)*np.nan
            return popt, pcov

        if len(w_tmp) <= 1:
            popt = np.zeros(self.nparam)*np.nan
            pcov = np.zeros(self.nparam)*np.nan
            return popt, pcov

        guess = self._guess(w_tmp, data_tmp)

        # enforce guess within bounds
        for i in range(len(guess)):
            if guess[i] < bounds[0][i]:
                guess[i] = bounds[0][i]
            elif guess[i] > bounds[1][i]:
                guess[i] = bounds[1][i]
            elif ~np.isfinite(guess[i]):
                guess[i] = 0.5*(bounds[1][i] - bounds[0][i])

        try:
            popt, pcov = curve_fit(
                    self.calculate,
                    w_tmp,
                    data_tmp,
                    sigma=sigma_tmp,
                    bounds=bounds,
                    p0=guess,
                    )

            pcov = pcov.diagonal().copy()

        except RuntimeError:
            # Occurs when curve_fit fails to converge
            popt = np.zeros(guess.size)*np.nan
            pcov = np.zeros(guess.size)*np.nan

        return popt, pcov

    def fit_cube(self, wave, data, var=None, bounds=None, wave_axis=2):
        """Fit cube using curve_fit.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        data : array-like
        var : array-like, default is None
            Variance array. Default is None, which means curve_fit is not
            supplied a variance array.
        bounds : list of lists
            Bounds supplied to curve_fit for each parameter.
         wave_axis : int
            Axis for wavelength can either be 0 or 2. Default is 2.

        Returns
        -------
        fit : np.ndarray
            Fit parameters returned by curve_fit.
        fit_err : np.ndarray
            Fit errors returned by curve_fit.

        """
        if wave_axis == 2:
            shp = (self.nparam, *data.shape[:2])
            fit = np.zeros(shp)*np.nan
            fit_err = np.zeros(shp)*np.nan
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if var is None:
                        var_tmp = None
                    else:
                        var_tmp = var[i, j]

                    fit[:, i, j], fit_err[:, i, j] = self.fit_spaxel(
                            wave, data[i, j], var_tmp, bounds)

        elif wave_axis == 0:
            shp = (self.nparam, *data.shape[1:])
            fit = np.zeros(shp)*np.nan
            fit_err = np.zeros(shp)*np.nan
            for i in range(data.shape[1]):
                for j in range(data.shape[2]):
                    if var is None:
                        var_tmp = None
                    else:
                        var_tmp = var[:, i, j]

                    fit[:, i, j], fit_err[:, i, j] = self.fit_spaxel(
                            wave, data[:, i, j], var_tmp, bounds)

        else:
            raise ValueError('Wave axis needs to be 0 or 2.')

        return fit, fit_err

    def _guess(self, wave, data, lambda_win=10.0):
        """Guess parameters using method of moments.

        This takes a window around each emission line then calculates the
        mean and standard deviation. For coupled lines it just uses the
        first line in the lines array to make a guess.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        data : array-like
        lambda_win : float
            Window to estimate parameters

        Returns
        -------
        guess : np.ndarray

        """
        dwave = wave[1] - wave[0]
        wave_left = wave - 0.5*dwave
        wave_right = wave + 0.5*dwave

        if self.baseline_order is None:
            guess = np.zeros(self.nlines + 2)
        else:
            guess = np.zeros(self.nlines + 2 + self.baseline_order + 1)

        # tmp_flux = np.zeros(self.nlines)*np.nan
        tmp_v = np.zeros(self.nlines)*np.nan
        tmp_vdisp = np.ones(self.nlines)*np.nan

        # guess[-1] =
        # guess[:self.nlines] = tmp_flux
        # guess[self.nlines] = tmp_v
        # guess[self.nlines+1] = np.mean(tmp_vdisp)
        for i, line in enumerate(self.lines):
            win = (
                (wave_right >= line[0] - lambda_win)
                & (wave_left <= line[0] + lambda_win)
                )
            if win.sum() <= 0:
                # if no valid pixels around emission line, then use initial
                # guesses
                continue

            win_data = data[win]
            win_wave = wave[win]

            # Guess flux.
            guess[i] = max(0.0, win_data.sum())

            # Calculate weights by normalising to minimum flux in array. Set
            # all weights to one if all are equal to the minimum weight.
            weights_v = win_data - np.min(win_data)
            if np.all(weights_v == 0.0):
                weights_v = np.ones(weights_v.shape)

            mean_wave = np.average(win_wave, weights=weights_v)
            tmp_v[i] = (mean_wave/line[0] - 1.0)*PhysicalConstants.C

            # TODO : Below is a biased estimator. But this will do for an
            # initial guess. Also, I just throw out all negative values
            # here, which may not be ideal but seems to work just fine.
            weights_vdisp = win_data.copy()
            weights_vdisp[win_data <= 0.0] = 0.0
            var_wave = np.average((win_wave - mean_wave)**2,
                                  weights=weights_vdisp)
            var_wave = var_wave - self.lsf_sigma**2
            tmp_vdisp[i] = max(1e-9, np.sqrt(var_wave))
            tmp_vdisp[i] *= PhysicalConstants.C/line[0]

        guess[self.nlines] = np.average(tmp_v, weights=guess[:self.nlines])
        guess[self.nlines+1] = np.average(tmp_vdisp, weights=guess[:self.nlines])

        return guess

    def _gas_model(self, wave, gas_param):
        """Gas model.

        Parameters
        ----------
        wave : array-like.
            Wavelength array.
        gas_param : array-like.
            DESCRIPTION.

        Returns
        -------
        model : TYPE
            DESCRIPTION.

        """
        model = np.zeros(len(wave))

        rel_lambda = 1.0 + gas_param[self.nlines]/PhysicalConstants.C
        rel_lambda_sigma = gas_param[self.nlines+1]/PhysicalConstants.C

        # add emission line contribution
        for i, line in enumerate(self.lines):
            # model first line
            line_wave = line[0]
            line_flux = gas_param[i]

            lam = rel_lambda*line_wave
            lam_sigma = rel_lambda_sigma*line_wave

            model += self._gas_line_model(wave, line_flux, lam, lam_sigma)

            nclines = len(line)//2
            for i in range(nclines):
                line_wave = line[1+2*i]
                factor = line[2+2*i]
                lam = rel_lambda*line_wave
                lam_sigma = rel_lambda_sigma*line_wave
                model += self._gas_line_model(
                        wave, factor*line_flux, lam, lam_sigma)

        return model

    def _gas_line_model(self, wave, flux, lam, lam_sigma):
        """Gas emission line model.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        flux : float
            Line flux.
        lam : float
            Line wavelength.
        lam_sigma : float
            Line standard deviation.

        Returns
        -------
        gas_line_model : np.ndarray

        """
        dwave = wave[1] - wave[0]
        wave_left = wave - 0.5*dwave
        wave_right = wave + 0.5*dwave

        var = lam_sigma**2 + self.lsf_sigma**2

        cdf_left = 0.5*erf((wave_left - lam)/np.sqrt(2.0*var))
        cdf_right = 0.5*erf((wave_right - lam)/np.sqrt(2.0*var))

        return flux*(cdf_right - cdf_left)

    def _baseline_model(self, wave, baseline_param):
        """Baseline model.

        Basline is assumed to follow a polynomial function of order
        self.baseline_order. if self.baseline_order is None this function is
        not called.

        Parameters
        ----------
        wave : array-like
            Wavelength array.
        baseline_param : array-like
            List of coefficients for polynomial model.

        Returns
        -------
        baseline : np.ndarray

        """
        wave_shft = wave - self.wave_ref
        baseline = np.ones(len(wave))*baseline_param[0]
        for i in range(self.baseline_order):
            baseline += baseline_param[1+i]*wave_shft**(1+i)

        return baseline
