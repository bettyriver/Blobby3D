"""Post analysis for Blobby3D output.

To be used after running Blobby3D to organise the Blobby3D output.


Attributes
----------
Classes: Blobby3D

@author: mathewvaridel

"""

from pathlib import Path
import numpy as np
import pandas as pd

from .meta import Metadata


class PostBlobby3D:

    def __init__(
            self, samples_path, data_path, var_path, metadata_path,
            save_maps=True, save_precon=True, save_con=True,
            nlines=1, nsigmad=2):
        """Blobby3D postprocess object.

        This provides a read and storage object for analysis.

        Parameters
        ----------
        samples_path : str or pathlib object
            Sample path. This can be either DNest4 samples or posterior
            samples.
        data_path : str or pathlib object
            Data cube path. Data cubes are assumed to be in whitespace
            separated text files. The values should be in row-major order.
        var_path : str or pathlib object
            Variance cube path. Variance cubes are assumed to be in whitespace
            separated text files. The values should be in row-major order.
        metadata_path : str or pathlib object
            Metadata file path. This path records the coordinates for your
            cube.
        save_maps : bool, optional
            Maps were saved to sampled by Blobby3D. The default is True.
        save_precon : bool, optional
            Preconvolved cubes were saved to samples by Blobby3D. The default
            is True.
        save_con : bool, optional
            Convolved cubes were saved to samples by Blobby3D. The default is
            True.
        nlines : int, optional
            Number of lines modelled during run. A coupled line count as one.
            Therefore modelling H-alpha and [NII] at 6548.1 A and 6583.1 A
            with a fixed ratio is counted as two lines in total. The default is
            1.
        nsigmad : int, optional
            The degree to which white and shot noise are modelled. The default
            is 2.

        Attributes
        ----------
        nsamples : int
            Number of samples recorded.
        max_blobs : int
            Maximum number of blobs allowed in run.
        metadata : Metadata object
            See Metadata class docstring for details.
        data : np.ndarray
        var : np.ndarray
        maps : np.ndarray
            2D maps of flux per line, velocity, and velocity dispersion for
            each sample..
        precon_cubes : np.ndarray
            Preconvolved model cubes for each sample.
        con_cubes : np.ndarray
            Convolved model cubes for each sample.
        global_param : pd.DataFrame
            Global parameters.
        blob_param :
            Parameters for each blob.

        """
        self._posterior_path = Path(samples_path)
        self._data_path = Path(data_path)
        self._var_path = Path(var_path)
        self._metadata_path = Path(metadata_path)
        self._nlines = nlines
        self._nsigmad = nsigmad

        # import data
        self.metadata = Metadata(self._metadata_path)
        self.data = np.loadtxt(data_path).reshape(self.metadata.naxis)
        self.var = np.loadtxt(var_path).reshape(self.metadata.naxis)

        # posterior samples
        samples = np.atleast_2d(np.loadtxt(samples_path))
        self.nsamples = samples.shape[0]

        if save_maps:
            self.maps = np.zeros((
                    self.nsamples,
                    self._nlines+2,
                    *self.metadata.naxis[:2]))
            map_shp = self.metadata.naxis[:2].prod()
            for s in range(self.nsamples):
                # Flux
                for ln in range(self._nlines):
                    self.maps[s, ln, :, :] = samples[
                            s, ln*map_shp:(1+ln)*map_shp
                            ].reshape(self.metadata.naxis[:2])

                # LoS velocity
                self.maps[s, self._nlines, :, :] = samples[
                        s, self._nlines*map_shp:(1+self._nlines)*map_shp
                        ].reshape(self.metadata.naxis[:2])

                # LoS vdisp
                self.maps[s, 1+self._nlines, :, :] = samples[
                        s, (1+self._nlines)*map_shp:(2+self._nlines)*map_shp
                        ].reshape(self.metadata.naxis[:2])

        if save_con:
            self.precon_cubes = np.zeros((self.nsamples, *self.metadata.naxis))
            st = save_maps*(2+self._nlines)*map_shp
            for s in range(self.nsamples):
                self.precon_cubes[s, :, :, :] = samples[
                        s, st:st+self.metadata.sz].reshape(self.metadata.naxis)

        if save_precon:
            self.con_cubes = np.zeros((self.nsamples, *self.metadata.naxis))
            st = save_maps*(2+self._nlines)*map_shp
            st += save_precon*self.metadata.sz
            for s in range(self.nsamples):
                self.con_cubes[s, :, :, :] = samples[
                        s, st:st+self.metadata.sz].reshape(self.metadata.naxis)

        # blob parameters
        st = save_maps*(2+self._nlines)*map_shp
        st += self.metadata.sz*(save_precon + save_precon)
        self.max_blobs = int(samples[0, st+1])

        global_names = ['WD', 'RADMIN', 'RADMAX', 'QMIN']
        for i in range(self._nlines):
            global_names += ['FLUX%iMU' % (i), 'FLUX%iSD' % (i)]
        global_names += [
                'NUMBLOBS',
                'XC', 'YC',
                'DISKFLUX', 'DISKMU',
                'VSYS', 'VMAX', 'VSLOPE', 'VGAMMA', 'VBETA'
                ]
        global_names += ['VDISP%i' % (i) for i in range(self._nsigmad)]
        global_names += ['INC', 'PA', 'SIGMA0', 'SIGMA1']
        global_param = np.concatenate((
                samples[:, st+2:st+7+2*self._nlines],
                samples[:, -13-self._nsigmad:]), axis=1)
        self.global_param = pd.DataFrame(global_param, columns=global_names)
        self.global_param.index.name = 'SAMPLE'

        if self.max_blobs > 0:
            blob_names = ['RC', 'THETAC', 'W', 'Q', 'PHI']
            blob_names += ['FLUX%i' % (i) for i in range(self._nlines)]
            n_bparam = len(blob_names)
            self.blob_param = np.zeros(
                    (self.nsamples*self.max_blobs, n_bparam))
            st_bparam = st+7+2*self._nlines
            end_bparam = -13-self._nsigmad
            for s in range(self.nsamples):
                row_st = self.max_blobs*s
                row_end = self.max_blobs*(s + 1)
                sblob_param = samples[s, st_bparam:end_bparam]
                sblob_param = sblob_param.reshape(
                        (self.max_blobs, n_bparam), order='F')
                self.blob_param[row_st:row_end, :] = sblob_param

            self.blob_param = pd.DataFrame(self.blob_param, columns=blob_names)

            self.blob_param['SAMPLE'] = np.repeat(
                    np.arange(self.nsamples), self.max_blobs)
            self.blob_param['BLOB'] = np.tile(
                    np.arange(self.max_blobs), self.nsamples)
            self.blob_param.set_index(['SAMPLE', 'BLOB'], inplace=True)
            self.blob_param = self.blob_param[self.blob_param['RC'] > 0.0]

    # def plot_global_marginalised(self, save_file=None):
    #     if isinstance(save_file, str):
    #         pdf_file = mpl.backends.backend_pdf.PdfPages(save_file)

    #     for key in self.global_param.keys():
    #         fig = plt.figure()
    #         ax = plt.gca()
    #         ax.hist(self.global_param[key])
    #         ax.set_title(r'%s' % (str(key)))
    #         if save_file:
    #             pdf_file.savefig(fig)
    #             plt.close()

    #     if isinstance(save_file, str):
    #         pdf_file.close()

    # def plot_map(self, ax, map_2d, **kwargs):
    #     """Plot individual map to a given axes object."""
    #     b3dplot.plot_map(ax, map_2d, **kwargs)

    # def setup_comparison_maps(
    #         self, figsize=(10.0, 10.0), log_flux=True, **kwargs):
    #     """Setup comparsion maps for a given sample."""
    #     fig, ax = b3dcomp.setup_comparison_maps(
    #             comp_shape=(2 + self._nlines, 3),
    #             map_shape=self.metadata.naxis[:2],
    #             figsize=figsize,
    #             **kwargs)

    #     return fig, ax

    # def add_comparison_maps(self, maps, ax, col, log_flux=False):
    #     """Add a column worth of maps to the comparison maps."""
    #     for line in range(self._nlines):
    #         if log_flux:
    #             with np.errstate(invalid='ignore', divide='ignore'):
    #                 flux_map = np.log10(maps[line, :, :])
    #                 flux_map[~np.isfinite(flux_map)] = np.nan
    #         else:
    #             flux_map = maps[line, :, :]
    #         self.plot_map(ax[line][col], flux_map, cmap=b3dplot.cmap.flux)
    #     self.plot_map(
    #             ax[self._nlines][col], maps[self._nlines, :, :],
    #             cmap=b3dplot.cmap.v)
    #     self.plot_map(
    #             ax[self._nlines+1][col], maps[self._nlines+1, :, :],
    #             cmap=b3dplot.cmap.vdisp)

    # def add_comparison_residuals(self, maps, ax, col, log_flux=False):
    #     """Add a column worth of maps to the comparison maps."""
    #     for line in range(self._nlines):
    #         if log_flux:
    #             flux_map = np.log10(maps[line, :, :])
    #         else:
    #             flux_map = maps[line, :, :]
    #         self.plot_map(ax[line][col], flux_map, cmap=b3dplot.cmap.residuals)
    #     self.plot_map(
    #             ax[self._nlines][col], maps[self._nlines, :, :],
    #             cmap=b3dplot.cmap.residuals)
    #     self.plot_map(
    #             ax[self._nlines+1][col], maps[self._nlines+1, :, :],
    #             cmap=b3dplot.cmap.residuals)

    # def update_comparison_clim(
    #         self, ax, cax,
    #         pct=100.0, vlim=None, absolute=False, **cb_kwargs):
    #     """Update each row colour scale for comparison plots."""
    #     if isinstance(ax, (list, tuple)):
    #         data = np.array([a.images[0]._A.data for a in ax])
    #         clim = b3dcomp.map_limits(data, pct, vlim, absolute)
    #         for a in ax:
    #             a.images[0].set_clim(clim)
    #         b3dcomp.colorbar(ax[-1].images[0], cax, clim, **cb_kwargs)
    #     else:
    #         data = ax.images[0]._A.data
    #         clim = b3dcomp.map_limits(data, pct, vlim, absolute)
    #         ax.images[0].set_clim(clim)
    #         b3dcomp.colorbar(ax.images[0], cax, clim, **cb_kwargs)

    # def update_comparison_colorbar(self, ax, **cb_kwargs):
    #     """
    #     Update colorbars in accordance with corresponding maps.

    #     ax : list of matplotlib.axes
    #         The right-most axis is the colorbar. The preceding maps are the
    #         maps with the corresponding clim.
    #     """
    #     clims = [a.images[0].get_clim() for a in ax[:-1]]
    #     for clim in clims:
    #         assert clim == clims[0]

    #     b3dcomp.colorbar(ax[-1], clims[0], **cb_kwargs)

    # def update_comparison_mask(self, mask, ax):
    #     """
    #     Update mask for maps.

    #     Parameters
    #     ----------
    #     ax : matplotlib.axes._subplots.AxesSubplot
    #     """
    #     if isinstance(ax, (list, tuple)):
    #         for i, a in enumerate(ax):
    #             a.images[0]._A = np.ma.masked_where(mask, a.images[0]._A)
    #     else:
    #         ax.images[0]._A = np.ma.masked_where(mask, ax.images[0]._A)
