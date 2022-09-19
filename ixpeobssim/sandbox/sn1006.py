import os
import numpy

from astropy.io import fits
import matplotlib.pyplot as plt
from ixpeobssim.utils.misc import pairwise
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube


#file_path = '/work/sn1006/acisf04391_repro_evt2.fits'

#hdul = fits.open(file_path)
#events = hdul[1]
#print (events.data['energy'])
#plt.hist (events.data['energy']/1e3, bins=250)
#plt.xlim((2,8))
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
SIMS = numpy.arange (10)
DU_ID = 3

file_path = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/sn1006'

def _file_name (sim_idx, du_id):
    return f'xpphotonlist_{sim_idx}__du{du_id}_photon_list_ixpesim_recon_simfmt_pmap.fits'
    
def _file_list():
    return [os.path.join(file_path, _file_name(sim_idx, DU_ID)) for  sim_idx in SIMS]


pmap = xBinnedPolarizationMapCube.from_file_list(_file_list())
pmap.plot_polarization_degree(arrows = None, num_sigma = 2.5)
pmap.plot_significance()
pmap.plot_mdp_map()
plt.show()