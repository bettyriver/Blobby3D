#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fitting kinematic moments.

@author: mathewvaridel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf

from .const import PhysicalConstants

# constants
C = PhysicalConstants.C


class SpectralModel:

    def __init__(self, lines, lsf_fwhm, baseline_order=None, wave_ref=0.0):
        self.lines = lines
        self.nlines = len(lines)

        self.lsf_fwhm = lsf_fwhm
        self.lsf_sigma = lsf_fwhm/np.sqrt(8.0*np.log(2.0))

        self.baseline_order = baseline_order

        self.wave_ref = wave_ref

        self.nparam = self.nlines + 2
        if self.baseline_order is not None:
            self.nparam += 1 + self.baseline_order

    def calculate(self, wavelength, *param):
        if self.baseline_order is None:
            model = self._gas_model(wavelength, param)
        else:
            model = self._gas_model(
                    wavelength,
                    param[:self.nlines+2])
            model += self._baseline_model(
                    wavelength,
                    param[-1-self.baseline_order:],
                    )

        return model

    def fit_spaxel(self, wavelength, data, var=None, bounds=None):
        if bounds is None:
            # line flux bounds
            bounds = [
                    [np.log(1e-9)]*self.nlines,
                    [np.inf]*self.nlines,
                    ]

            # v, vdisp bounds
            bounds[0] += [-np.inf, np.log(1e-9)]
            bounds[1] += [np.inf, np.inf]

            if self.baseline_order is not None:
                bounds[0] += [-np.inf]*(self.baseline_order + 1)
                bounds[1] += [np.inf]*(self.baseline_order + 1)

        if (var is None) & np.any(data != 0.0):
            data_valid = np.isfinite(data)
            data_tmp = data[data_valid]
            w_tmp = wavelength[data_valid]
            sigma_tmp = None
        elif (var is None) & np.all(data == 0.0):
            popt = np.zeros(self.nparam)*np.nan
            pcov = np.zeros(self.nparam)*np.nan
            return popt, pcov
        elif np.any(var > 0.0):
            data_valid = ((var > 0.0) & np.isfinite(var) & np.isfinite(data))
            data_tmp = data[data_valid]
            w_tmp = wavelength[data_valid]
            sigma_tmp = np.sqrt(var[data_valid])
        else:
            popt = np.zeros(self.nparam)*np.nan
            pcov = np.zeros(self.nparam)*np.nan
            return popt, pcov

        if len(w_tmp) <= 1:
            popt = np.zeros(self.nparam)*np.nan
            pcov = np.zeros(self.nparam)*np.nan
            return popt, pcov

        try:
            guess = self._guess(w_tmp, data_tmp)
        except (ZeroDivisionError, IndexError):
            print(w_tmp, data_tmp)

        # enforce guess within bounds
        for i in range(len(guess)):
            if guess[i] < bounds[0][i]:
                guess[i] = bounds[0][i]
            elif guess[i] > bounds[1][i]:
                guess[i] = bounds[1][i]
            elif ~np.isfinite(guess[i]):
                guess[i] = 0.5*(bounds[1][i] - bounds[0][i])

        try:
            # print('GUESS', guess, bounds)
            # print(w_tmp, data_tmp, sigma_tmp)
            popt, pcov = curve_fit(
                    self.calculate,
                    w_tmp,
                    data_tmp,
                    sigma=sigma_tmp,
                    bounds=bounds,
                    p0=guess,
                    )

            pcov = pcov.diagonal().copy()

            # convert back to linear space
            popt[:self.nlines] = np.exp(popt[:self.nlines])
            popt[self.nlines+1] = np.exp(popt[self.nlines+1])
            pcov[:self.nlines] = np.exp(pcov[:self.nlines])
            pcov[self.nlines+1] = np.exp(pcov[self.nlines+1])
        except RuntimeError:
            # Occurs when curve_fit fails to converge
            popt = np.zeros(guess.size)*np.nan
            pcov = np.zeros(guess.size)*np.nan

        return popt, pcov

    def fit_cube(self, wavelength, data, var=None, wave_axis=2, **fit_kwargs):
        if wave_axis == 2:
            fit = np.zeros((self.nparam, *data.shape[:2]))*np.nan
            fit_err = np.zeros((self.nparam, *data.shape[:2]))*np.nan
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if var is not None:
                        fit[:, i, j], fit_err[:, i, j] = self.fit_spaxel(
                            wavelength, data[i, j, :], var[i, j, :],
                            **fit_kwargs)
                    else:
                        fit[:, i, j], fit_err[:, i, j] = self.fit_spaxel(
                            wavelength, data[i, j, :], **fit_kwargs)
        elif wave_axis == 0:
            fit = np.zeros((self.nparam, *data.shape[1:]))*np.nan
            fit_err = np.zeros((self.nparam, *data.shape[1:]))*np.nan
            for i in range(data.shape[1]):
                for j in range(data.shape[2]):
                    if var is not None:
                        fit[:, i, j], fit_err[:, i, j] = self.fit_spaxel(
                                wavelength, data[:, i, j], var[:, i, j],
                                **fit_kwargs)
                    else:
                        fit[:, i, j], fit_err[:, i, j] = self.fit_spaxel(
                                wavelength, data[:, i, j], **fit_kwargs)
        else:
            raise ValueError('Wave axis needs to be 0 or 2.')

        return fit, fit_err

    def _guess(self, wavelength, data, lambda_win=10.0):
        """Guess parameters using method of moments.

        data : array-like
        lambda_win : floatp
            Window to estimate parameters
        """
        dwave = wavelength[1] - wavelength[0]
        wave_left = wavelength - 0.5*dwave
        wave_right = wavelength + 0.5*dwave

        if self.baseline_order is None:
            guess = np.zeros(self.nlines + 2)
        else:
            guess = np.zeros(self.nlines + 2 + self.baseline_order + 1)

        tmp_flux = np.ones(self.nlines)*1e-9
        tmp_v = np.zeros(self.nlines)
        tmp_vdisp = np.ones(self.nlines)*1e-9*C/self.lines[0]

        guess[:self.nlines] = tmp_flux.copy()
        guess[self.nlines] = np.mean(tmp_v)
        guess[self.nlines+1] = np.mean(tmp_vdisp)
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
            win_wave = wavelength[win]

            tmp_flux[i] = max(1e-9, win_data.sum())
            guess[i] = max(1e-9, win_data.sum())

            weights = win_data - np.min(win_data)
            if np.all(weights == 0.0):
                weights = np.ones(weights.shape)

            mean_wave = np.average(win_wave, weights=weights)
            tmp_v[i] = (mean_wave/line[0] - 1.0)*C

            # TODO : Think below is wrong + it's a biased estimator
            tmp_vdisp[i] = np.sqrt(np.average(
                    (win_wave - mean_wave)**2, weights=weights))
            tmp_vdisp[i] = max(
                    1e-9,
                    np.sqrt(tmp_vdisp[i]**2 - self.lsf_sigma**2)
                    )
            tmp_vdisp[i] *= C/line[0]
            tmp_vdisp[i] = 20.0  # just for testing

        guess[self.nlines] = np.average(tmp_v, weights=tmp_flux)
        guess[self.nlines+1] = np.average(tmp_vdisp, weights=tmp_flux)

        # convert flux and vdisp to log
        guess[:self.nlines] = np.log(guess[:self.nlines])
        guess[self.nlines+1] = np.log(guess[self.nlines+1])

        return guess

    def _gas_model(self, wave, gas_param):
        model = np.zeros(len(wave))

        rel_lambda = 1.0 + gas_param[self.nlines]/C
        rel_lambda_sigma = np.exp(gas_param[self.nlines+1])/C

        # add emission line contribution
        for i, line in enumerate(self.lines):
            # model first line
            line_wave = line[0]
            line_flux = gas_param[i]

            lam = rel_lambda*line_wave
            lam_sigma = rel_lambda_sigma*line_wave

            model += self._gas_line_model(
                    wave, line_flux, lam, lam_sigma)

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
        dwave = wave[1] - wave[0]
        wave_left = wave - 0.5*dwave
        wave_right = wave + 0.5*dwave

        var = lam_sigma**2 + self.lsf_sigma**2

        cdf_left = 0.5*erf((wave_left - lam)/np.sqrt(2.0*var))
        cdf_right = 0.5*erf((wave_right - lam)/np.sqrt(2.0*var))
        return np.exp(flux)*(cdf_right - cdf_left)  # /dx

    def _baseline_model(self, wave, baseline_param):
        wave_shft = wave - self.wave_ref
        baseline = np.ones(len(wave))*baseline_param[0]
        for i in range(self.baseline_order):
            baseline += baseline_param[1+i]*wave_shft**(1+i)
        return baseline


if __name__ == '__main__':
    import os
    from blobby3d import Blobby3D

    root = r'/Users/mathewvaridel/Google Drive/Code/Uni/Blobby3D/Examples/106717'
    b3d = Blobby3D(
        samples_path=os.path.join(root, 'posterior_sample.txt'),
        data_path=os.path.join(root, 'data.txt'),
        var_path=os.path.join(root, 'var.txt'),
        metadata_path=os.path.join(root, 'metadata.txt'),
        nlines=2)

    sm = SpectralModel(
            lines=[[6562.81], [6583.1, 6548.1, 0.3333]],
            lsf_fwhm=1.61
            )

    wave = np.linspace(b3d.r_lim[0] + 0.5*b3d.dr,
                       b3d.r_lim[1] - 0.5*b3d.dr,
                       b3d.naxis[2])

    i, j = 15, 15

    guess = sm._guess(wave, b3d.data[i, j, :])
    guess_model = sm._gas_model(wave, guess)

    fit, fit_err = sm.fit_spaxel(wave, b3d.data[i, j, :])
    fit_model = sm.calculate(wave, *fit)

    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(wave, guess_model)
    plt.plot(wave, fit_model)
    plt.plot(wave, b3d.data[i, j, :])

    fit_cube, fit_cube_err = sm.fit_cube(wave, b3d.data)
