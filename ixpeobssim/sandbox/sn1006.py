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
SIMS = numpy.arange (50)
DU_ID = 2
DU_IDS = [1,2,3]

file_path = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/sn1006'

def _file_name (sim_idx, du_id):
    return f'xpphotonlist_{sim_idx}__du{du_id}_photon_list_ixpesim_recon_simfmt_pmap.fits'
    
def _file_list():
    return [os.path.join(file_path, _file_name(sim_idx, DU_ID)) for  sim_idx in SIMS]

pmap = xBinnedPolarizationMapCube.from_file_list(_file_list())
#pmap.plot_polarization_degree(arrows = None, num_sigma = 1)
#pmap.plot_significance()
mask = pmap.MDP_99<0.2
mask = mask[0]
pmap.plot_mdp_map()
plt.show()
input()

def statistics():
    smaps = []
    #pmap = xBinnedPolarizationMapCube.from_file_list([_file_list()[0]])
    for j in range (len(_file_list())):
        signif = xBinnedPolarizationMapCube.from_file_list([_file_list()[j]]).SIGNIF[0]*mask
        smaps.append(signif)
    return (smaps)

#def img_mean(maps, param):
#    numpy.average (maps[j].SIGNIF[0])

smaps = statistics()
smaps = numpy.array(smaps)
plt.figure('average')
plt.imshow(numpy.average(smaps, axis=0))
plt.colorbar()
plt.figure('std')
plt.imshow(numpy.std((smaps), axis=0))
plt.colorbar()
plt.figure('maximum')
plt.imshow(numpy.max((smaps), axis=0))
plt.colorbar()
plt.show()
input()
plt.close()
bins = numpy.linspace(1,6, num=50)
#for j in range (25):
#    for i in range (25):
#        plt.hist(smaps[:,j,i], bins = bins)
plt.hist(smaps.reshape(31250), bins=bins)
plt.show()




