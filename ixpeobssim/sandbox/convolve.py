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
from matplotlib.colors import LogNorm
from ixpeobssim.sandbox.csgro.friendFiles import xEventFileFriend
from ixpeobssim.evt.event import xEventFile

IRF_NAME = 'ixpe:obssim:v11'
X_REF, Y_REF = 228.48133, -59.13578
E_MIN, E_MAX = 2., 5.5
NUM_PIXELS = 15
IMG_SIZE = 8. # arcmin
CONV_NSIDE = 15 # pixels

PIXEL_SIZE = arcmin_to_arcsec(IMG_SIZE / NUM_PIXELS)

l2_path = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/observations/msh_1552/01001101/event_l2'
l1_path = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/observations/msh_1552/01001101/event_l1'
work = 'work'
workdir = os.path.join(l2_path, work)



def _l2_file_list():
    return glob(os.path.join(l2_path, '*v01.fits'))
    
def _l1_file_list():
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


def plot_convolved_pmap(size, pmap, sigma = 3, suffix = ''):
    kernelsize = size
    #pmap = xBinnedPolarizationMapCube.from_file_list(_map_file_list('pmap', '200_2_8'))
    kernel = circular_kernel(kernelsize)
    pmap.convolve(kernel)
    pmap.plot_significance()
    pmap.plot_polarization_degree(arrows = None, num_sigma = sigma)
    #psf = plt.Circle((10,10), size/2., color = 'r', label = 'Convolution kernel', fill = False)
    #plt.add_patch(psf)
    plt.savefig(os.path.join(workdir, 'pd_200_2_8_' + str(size) + suffix))
    plt.close()
    pmap.plot_significance()
    plt.savefig(os.path.join(workdir, 'sig_200_2_8_' + str(size) + suffix))
    plt.close()
    pmap.plot_mdp_map()
    plt.savefig(os.path.join(workdir, 'mdp99_200_2_8_' + str(size) + suffix))
    plt.close()
    pmap.plot_polarization_angle(num_sigma = sigma)
    plt.savefig(os.path.join(workdir, 'pa_200_2_8_' + str(size) + suffix))
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

def cleanup(l2, l1, y0 = 0.85):
    ''' Do the filtering of l2 events based on l1 parameters.
    Based on carmelo's. Full paths to the single data file have to be
    specified.
    '''
    mySource = xEventFileFriend(l2, l1)
    tgr_id_l2 = mySource.l2value("TRG_ID")
    tgr_id_l1 = mySource.l1value("TRG_ID")
    #assert numpy.sum(tgr_id_l1 - tgr_id_l2) = 0, 'Error: Ids differ'
    ra, dec  = mySource.sky_position_data()
    energy   = mySource.energy_data()
    num_pix  = mySource.l1value("NUM_PIX")
    trk_size = mySource.l1value("TRK_SIZE")
    trk_bord = mySource.l1value("TRK_BORD")
    evt_fra  = mySource.l1value("EVT_FRA")
    evt_x = mySource.l1value('ABSX')
    evt_y = mySource.l1value('ABSY')
    cut_enerange = numpy.logical_and(energy>2, energy<8)
    x = numpy.linspace(2,8,10)
    y_cut = (x<4)*(0.07*(x-4)) + (0.008*(x-4))*(x>=4) + y0
    cut_efra = evt_fra > ((energy<4)*(0.07*(energy-4)) + (0.008*(energy-4))*(energy>=4) + y0)
    y_test = 0.8*(1-numpy.exp(-(x-0.23)/1.1)) + 0.0055
    # Diagnostic plots
    plt.figure("EvtFravsE")
    plt.hist2d(energy, evt_fra,
           range=((0,10),(0,1.5)), bins=(100, 150), norm=LogNorm())
    plt.plot(x, y_cut , "g") # black h line 
    plt.plot(x, y_test , "b") # black h line 
    plt.xlabel("Energy [keV]")
    plt.ylabel("Fraction of Energy in the main track")
    plt.figure("RaDec_all", tight_layout=True, figsize=(16,5))

    plt.subplot(131)
    plt.title("All evt")
    plt.hist2d(ra[cut_enerange],dec[cut_enerange],
           bins=100, norm=LogNorm())
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.colorbar()

    plt.subplot(132)
    plt.title("Sgn (efra cut) evt")
    plt.hist2d(ra[numpy.logical_and(cut_enerange, cut_efra)],
           dec[numpy.logical_and(cut_enerange, cut_efra)],
           bins=100, norm=LogNorm())
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.colorbar()

    plt.subplot(133)
    plt.title("Bkg (efra cut) evt DU")
    plt.hist2d(ra[numpy.logical_and(cut_enerange, numpy.logical_not(cut_efra))],
           dec[numpy.logical_and(cut_enerange, numpy.logical_not(cut_efra))],
           bins=100, norm=LogNorm())
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.colorbar()

    plt.figure("NUMPIX")
    plt.hist2d(energy[cut_efra], num_pix[cut_efra],
           range=((0,10),(0,1000)), bins=(100, 500), norm=LogNorm())
    plt.plot(x, 130 + 30*(x-2))
    plt.xlabel("Energy [keV]")
    plt.ylabel("NUM_PIX")

    plt.figure("absorption point")
    plt.subplot(221)
    plt.hist2d(evt_x[cut_efra], evt_y[cut_efra], bins=(100, 100))
    plt.xlabel("X absorption point (clean)")
    plt.ylabel("Y absorption point (clean)")
    plt.colorbar()
    plt.subplot(222)
    plt.hist2d(evt_x[~cut_efra], evt_y[~cut_efra], bins=(100, 100))
    plt.xlabel("X absorption point (bkg)")
    plt.ylabel("Y absorption point (bkg)")
    plt.colorbar()
    plt.subplot(223)
    plt.hist(evt_x[~cut_efra], bins=100)
    plt.xlabel('X absorption point')
    plt.ylabel('counts')
    plt.xlim((-7,7))
    plt.subplot(224)
    plt.hist(evt_y[~cut_efra], bins=100)
    plt.xlabel('Y absorption point')
    plt.ylabel('counts')
    plt.xlim((-7,7))
    

    plt.figure("TRKBORD")
    plt.hist2d(energy[cut_efra], trk_bord[cut_efra],
           range=((0,10),(0,50)), bins=(100, 50), norm=LogNorm())
    plt.plot(x, 130 + 30*(x-2))
    plt.xlabel("Energy [keV]")
    plt.ylabel("TRKBORD")
    plt.show()
    print (f'Are you happy for the results? input Yes or new y0 value [{y0}]')
    a = input()
    if a == 'Y':
        return (cut_efra)
    else :
        a = numpy.float64(a)
        cleanup(l2, l1, y0 = a)

def fits_filter (fits_file, mask, suffix = 'clean'):
    event_file = xEventFile(fits_file)
    event_file.filter(mask)
    event_file.write(fits_file + suffix, overwrite = 'True')

DU = 3
# # Signal extraction and display scripts
good_events = cleanup(glob(f'{l2_path}/*det{DU}*v01.fits'), glob(f'{l1_path}/*det{DU}*v01.fits'))

#print (numpy.shape(good_events))
#print (numpy.sum(good_events))
#input()
#fits_filter(glob(f'{l2_path}/*det{DU}*v01.fits')[0], good_events, 'clean')
#_create_pmaps(200, 2, 8)
#plot_convolved_pmap(15, xBinnedPolarizationMapCube.from_file_list(_map_file_list('clean_pmap', '')), suffix = '_signal')
#pmap_dirs = glob(f'{l2_path}/ixpe01001101_det*clean_pmap.fits')
#reg_dir = os.path.join(l2_path, 'background2.reg')
#plot_convolved_pmap(15, subtract(pmap_dirs, reg_dir))

# # Noise extraction and display scritps
#bad_events = ~good_events
#bad_events = ~cleanup(glob(f'{l2_path}/*det{DU}*v01.fits'), glob(f'{l1_path}/*det{DU}*v01.fits'))
#fits_filter(glob(f'{l2_path}/*det{DU}*v01.fits')[0], bad_events, 'noise')
#plot_convolved_pmap(15, xBinnedPolarizationMapCube.from_file_list(_map_file_list('noise_pmap', '')), suffix = '_noise')

