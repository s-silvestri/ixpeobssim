import os

from astropy.io import fits
import numpy
import scipy

from glob import glob
from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube
from ixpeobssim.core.hist import xHistogram1d
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.deconvolution import circular_kernel
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.units_ import arcmin_to_arcsec
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, save_gcf
from ixpeobssim.utils.os_ import rm

IRF_NAME = 'ixpe:obssim:v11'
X_REF, Y_REF = 83.633, 22.0144
E_MIN, E_MAX = 2., 5.5
NUM_PIXELS = 15
IMG_SIZE = 8. # arcmin
CONV_NSIDE = 15 # pixels

PIXEL_SIZE = arcmin_to_arcsec(IMG_SIZE / NUM_PIXELS)

l2_path = '/work/msh_1552/01001101/event_l2'
work = 'work'
workdir = os.path.join(l2_path, work)


def _l2_file_list():
    return glob(os.path.join(l2_path, '*v01.fits'))

def _map_file_list(type, suffix):
    return glob(f'{l2_path}/*{type}*{suffix}.fits')


def _create_pmaps(nside, emin, emax, img_size=IMG_SIZE, overwrite=True):
    """
    """
    pixel_size = arcmin_to_arcsec(img_size / nside)
    suffix = 'pmap_%d_%d_%d' % (nside, emin, emax)
    kwargs = dict(npix=nside, pixsize=pixel_size, xref=X_REF, yref=Y_REF,
        emin=emin, emax=emax, irfname=IRF_NAME, suffix=suffix)
    pipeline.xpbin(*_l2_file_list(), algorithm='PMAP', overwrite=overwrite, **kwargs)


def plot_convolved_pmap(size):
    #_create_pmaps(200,2,8)
    kernelsize = size
    pmap = xBinnedPolarizationMapCube.from_file_list(_map_file_list('pmap', '200_2_8'))
    kernel = circular_kernel(kernelsize)
    pmap.convolve(kernel)
    pmap.plot_significance()
    pmap.plot_polarization_degree(arrows = None, num_sigma = 5)
    psf = plt.Circle((10,10), size/2., color = 'r', label = 'Convolution kernel', fill = False)
    plt.add_patch(psf)
    plt.savefig(os.path.join(workdir, 'pd_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_significance()
    plt.savefig(os.path.join(workdir, 'sig_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_mdp_map()
    plt.savefig(os.path.join(workdir, 'mdp99_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_polarization_angle(num_sigma = 5)
    plt.savefig(os.path.join(workdir, 'pa_200_2_8_' + str(size)))
    plt.close()
    #_write_pmap(pmap, _file_list('pmap_200_2_8_kernel'))

_create_pmaps(200, 2, 8)
plot_convolved_pmap(15)