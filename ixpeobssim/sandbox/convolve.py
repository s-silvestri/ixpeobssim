import os

from astropy.io import fits
import numpy
import scipy

from astropy.io import fits
from regions import Regions
from astropy.wcs import WCS
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
X_REF, Y_REF = 228.48133, -59.13578
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
    kwargs = dict(npix=nside, pixsize=pixel_size,
        emin=emin, emax=emax, irfname=IRF_NAME, suffix=suffix)
    pipeline.xpbin(*_l2_file_list(), algorithm='PMAP', overwrite=overwrite, **kwargs)


def plot_convolved_pmap(size, pmap):
    kernelsize = size
    #pmap = xBinnedPolarizationMapCube.from_file_list(_map_file_list('pmap', '200_2_8'))
    kernel = circular_kernel(kernelsize)
    pmap.convolve(kernel)
    pmap.plot_significance()
    pmap.plot_polarization_degree(arrows = None, num_sigma = 3)
    #psf = plt.Circle((10,10), size/2., color = 'r', label = 'Convolution kernel', fill = False)
    #plt.add_patch(psf)
    plt.savefig(os.path.join(workdir, 'pd_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_significance()
    plt.savefig(os.path.join(workdir, 'sig_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_mdp_map()
    plt.savefig(os.path.join(workdir, 'mdp99_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_polarization_angle(num_sigma = 3)
    plt.savefig(os.path.join(workdir, 'pa_200_2_8_' + str(size)))
    plt.close()
    pmap.plot_stokes_parameters()
    plt.show()
    #_write_pmap(pmap, _file_list('pmap_200_2_8_kernel'))

def subtract (polcube, dark):
    ''' both dark and polcube must be strings that point to the file dir
    dark is a region, but one day it might be a real dark.
    returns the subtracted pmap
    '''
    pmap = xBinnedPolarizationMapCube.from_file_list(polcube)
    dark_reg = Regions.read(os.path.join(l2_path, dark), format='ds9')
    print (dark_reg[0])
    print (dark_reg[1])
    #input()
    # Get the WCS
    hdul = fits.open(glob(f'{l2_path}/ixpe01001101_det1_evt2_v01_pmap_200_2_8.fits')[0]) 
    w = WCS(hdul[1].header)
    #Create a dark mask for the region   
    dark_mask_1 = dark_reg[0].to_pixel(w).to_mask().to_image(numpy.shape(pmap.SIGNIF[0]))
    dark_mask_2 = dark_reg[1].to_pixel(w).to_mask().to_image(numpy.shape(pmap.SIGNIF[0]))
    dark_mask = dark_mask_1+dark_mask_2
    figname = 'mask'
    plt.figure(figname)
    plt.imshow(dark_mask, origin='lower')
    plt.colorbar()
    plt.savefig(os.path.join(workdir, figname))
    figname = 'masked region'
    plt.figure (figname)
    plt.imshow(pmap.I[0]*dark_mask, origin='lower')
    plt.colorbar()
    plt.savefig(os.path.join(workdir, figname))
    figname = 'masked region (PD)'
    plt.figure (figname)
    plt.imshow(pmap.PD[0]*dark_mask, origin='lower')
    plt.colorbar()
    plt.savefig(os.path.join(workdir, figname))
    plt.show()
    #Create darks
    dark_I = numpy.median((pmap.I[0]*dark_mask)[dark_mask>0])
    dark_Q = numpy.median((pmap.Q[0]*dark_mask)[dark_mask>0])
    dark_U = numpy.median((pmap.U[0]*dark_mask)[dark_mask>0])
    dark_PD = numpy.median((pmap.PD[0]*dark_mask)[dark_mask>0])
    # print some stuff
    print (f'Median dark rate = {dark_I}\n')
    print (f'Median counts= {numpy.median(pmap.I[0])}\n')
    print (f'Median mask PD = {dark_PD}\n')

    #Subtract
    pmap.I[0] -= dark_I
    pmap.Q[0] -= dark_Q
    pmap.U[0] -= dark_U
    return (pmap)

#_create_pmaps(200, 2, 8)
#plot_convolved_pmap(15, xBinnedPolarizationMapCube.from_file_list(_map_file_list('pmap', '200_2_8')))
pmap_dirs = glob(f'{l2_path}/ixpe01001101_det*_pmap_200_2_8.fits')
reg_dir = os.path.join(l2_path, 'background2.reg')
plot_convolved_pmap(15, subtract(pmap_dirs, reg_dir))