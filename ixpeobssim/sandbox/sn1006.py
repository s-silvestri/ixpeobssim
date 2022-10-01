import os
import numpy
import glob

from astropy.io import fits
from regions import Regions
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from ixpeobssim.utils.misc import pairwise
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube
from ixpeobssim.binning.misc import xBinnedPulseProfile, xEventBinningCMAP


#file_path = '/work/sn1006/acisf04391_repro_evt2.fits'

#hdul = fits.open(file_path)
#events = hdul[1]
#print (events.data['energy'])
#plt.hist (events.data['energy']/1e3, bins=250)
#plt.xlim((2,8))
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
SIMS = numpy.arange (200)
SUBSIMS = numpy.arange(50)
DU_ID = 2
DU_IDS = [1,2,3]

file_path = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/sn1006'

def _file_name (sim_idx, du_id, suffix):
    return f'xpphotonlist_{sim_idx}__du{du_id}_photon_list_ixpesim_recon_simfmt_{suffix}.fits'

#def _file_name (sim_idx, du_id):
#    return f'xpphotonlist_{sim_idx}__du{du_id}_pmap.fits'
    
def _file_list(suffix):
    return [os.path.join(file_path, _file_name(sim_idx, DU_ID, suffix)) for sim_idx in SIMS]

def _file_list_small(suffix):
    return [os.path.join(file_path, _file_name(sim_idx, DU_ID, suffix)) for sim_idx in SUBSIMS]

def _glob_pmap(sim_n, suffix):
    return (glob.glob(f'{file_path}/xpphotonlist_{sim_n}_*{suffix}.fits'))



#pmap = xBinnedPolarizationMapCube.from_file_list(_file_list_small())
#pmap.plot_polarization_degree(arrows = None, num_sigma = 1)
#pmap.plot_significance()
#mask = pmap.MDP_99<0.2
#mask = mask[0]
#pmap.plot_stokes_parameters()
#pmap.plot_mdp_map()
#pmap.plot_polarization_angle()
#print (numpy.sum(pmap.Q)/numpy.sum(pmap.COUNTS))
#print (numpy.sum(pmap.U)/numpy.sum(pmap.COUNTS))

def _shift_wcs(crval, shift):
    wcs_input_dict = {
    'CTYPE1': 'RA---TAN', 
    'CUNIT1': 'deg', 
    'CDELT1': -0.0014583333333333, 
    'CRPIX1': 100.5 - shift, 
    'CRVAL1': crval[0], 
    'NAXIS1': 200,
    'CTYPE2': 'DEC--TAN', 
    'CUNIT2': 'deg', 
    'CDELT2': 0.0014583333333333, 
    'CRPIX2': 100.5 + shift, 
    'CRVAL2': crval[1], 
    'NAXIS2': 200, 
    'CTYPE3': '',
    'CUNIT3': '',
    'CRPIX3': 0.,
    'CRVAL3': 0.,
    'NAXIS3': 1
    }
    return WCS(wcs_input_dict)



pmap = xBinnedPolarizationMapCube.from_file_list(_file_list('pmap'))
# Load the WCS from a file
w = WCS(fits.open(os.path.join(file_path, _file_name(1,1,'pmap')))[1].header)
# region in WCS format
background = Regions.read(os.path.join(file_path, 'Background.reg'), format='ds9')
signal = Regions.read(os.path.join(file_path, 'signal.reg'), format='ds9')
# convert the region in pixel coordinates using the WCS
background = background[1].to_pixel(w)
signal = signal[1].to_pixel(w)
# make a mask for the background
bkg = background.to_mask()
sig = signal.to_mask()
# convert the mask to image shape
shape = numpy.shape(pmap.SIGNIF[0])
bkg_mask = bkg.to_image(shape)
sig_mask = sig.to_image(shape)

#input()
#plt.show()
#plt.figure('Mask')
#plt.imshow(mask, origin = 'lower')
#plt.show()
plt.figure('I map')
plt.imshow(pmap.I[0], origin = 'lower')
plt.colorbar()
plt.figure('Background')
plt.imshow(pmap.I[0]*bkg_mask, origin = 'lower')
plt.colorbar()
plt.figure('Signal')
plt.imshow(pmap.I[0]*sig_mask, origin = 'lower')
plt.colorbar()
plt.show()
print (f'Signal counts: {numpy.sum(pmap.I[0]*sig_mask)}')
print (f'Background counts: {numpy.sum(pmap.I[0]*bkg_mask)}')

reg = Regions.read(os.path.join(file_path, 'Background.reg'), format='ds9')
fluxes = []
for j in range (-10, 30):
    window = _shift_wcs([225.9444625755, -41.725407320703], j)
    moving_banana = reg[1].to_pixel(window).to_mask().to_image(numpy.shape(pmap.SIGNIF[0]))
    counts = (pmap.I[0]*moving_banana)
    fluxes.append(numpy.sum(counts))

plt.plot(fluxes)
plt.xlabel('Total counts')
plt.ylabel('radius')
plt.show()
input()
plt.figure('Inner end')
plt.imshow (pmap.I[0]*reg[1].to_pixel( _shift_wcs([225.9444625755, -41.725407320703], -10)).to_mask().to_image(numpy.shape(pmap.SIGNIF[0])), origin = 'lower')
plt.colorbar()
plt.figure('Outer end')
plt.imshow (pmap.I[0]*reg[1].to_pixel( _shift_wcs([225.9444625755, -41.725407320703], 30)).to_mask().to_image(numpy.shape(pmap.SIGNIF[0])), origin = 'lower')
plt.colorbar()
plt.show()


def statistics():
    smaps = []
    #pmap = xBinnedPolarizationMapCube.from_file_list([_file_list()[0]])
    for j in range (len(_file_list())):
        signif = xBinnedPolarizationMapCube.from_file_list([_file_list()[j]]).SIGNIF[0]*mask
        smaps.append(signif)
    return (smaps)

def statistics_3du():
    smaps = []
    for j in SIMS:
        signif = xBinnedPolarizationMapCube.from_file_list(_glob_pmap(j)).SIGNIF[0]*mask
        smaps.append(signif)
    return (smaps)
'''
smaps = statistics_3du()
smaps = numpy.array(smaps)
#plt.figure('average')
#plt.imshow(numpy.average(smaps, axis=0), origin = 'lower')
#plt.colorbar()
#plt.figure('std')
#plt.imshow(numpy.std((smaps), axis=0), origin = 'lower')
#plt.colorbar()
#plt.figure('maximum')
#plt.imshow(numpy.max((smaps), axis=0), origin = 'lower')
#plt.colorbar()
sigmas = 2
overthresh = 100.*(smaps>sigmas).sum(axis=0)/len(smaps)
plt.figure ('Over 2sigma threshold (percentage)')
plt.imshow(overthresh, origin = 'lower')
plt.colorbar()
sigmas = 2.5
overthresh = 100.*(smaps>sigmas).sum(axis=0)/len(smaps)
plt.figure ('Over 2.5 sigma threshold (percentage)')
plt.imshow(overthresh, origin = 'lower')
plt.colorbar()
plt.show()
input()
plt.close()
bins = numpy.linspace(1,6, num=50)
#for j in range (25):
#    for i in range (25):
#        plt.hist(smaps[:,j,i], bins = bins)
plt.hist(smaps.reshape(88200), bins=bins)
plt.ylabel ('counts')
plt.xlabel('sigma')
plt.show()
'''



