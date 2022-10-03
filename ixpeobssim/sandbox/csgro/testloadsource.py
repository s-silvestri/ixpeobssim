import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ixpeobssim.sandbox.csgro.friendFiles import xEventFileFriend
from ixpeobssim.evt.subselect import xEventSelect

# For CasA DU1
#commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CasA/obs01001301/"
#l1path = commonpath+"event_l1/ixpe01001301_det1_evt1_v01.fits"
#l2path = commonpath+"event2_corr/ixpe01001301E_det1_evt2_v02_picorr.fits"

# For CenX3
commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CenX3/"
l1path = []
l2path = []
for DuId in [1]:
    l1path.append(commonpath+"event_l1/ixpe01006501_det%d_evt1_v02.fits" % DuId)
    l2path.append(commonpath+"event_l2/ixpe01006501_det%d_evt2_v02.fits" % DuId)


        
myFiles = xEventFileFriend(l2path, l1path)

tgr_id_l2 = myFiles.l2value("TRG_ID")
tgr_id_l1 = myFiles.l1value("TRG_ID")
print("Trigger Id diff check, must be zero, and the winner is:",
      sum(tgr_id_l1 - tgr_id_l2))

energy   = myFiles.energy_data()
num_clu  = myFiles.l1value("NUM_CLU")
num_pix  = myFiles.l1value("NUM_PIX")
trk_size = myFiles.l1value("TRK_SIZE")
trk_bord = myFiles.l1value("TRK_BORD")
evt_fra  = myFiles.l1value("EVT_FRA")
w_mom    = myFiles.l2value("W_MOM")
ra, dec  = myFiles.sky_position_data()
ene_dens = energy/trk_size
pix_fra  = trk_size/num_pix

cut_enerange = numpy.logical_and(energy>2, energy<8)
my_cut       = evt_fra>0.8


plt.ion()

plt.figure("Energy")
plt.hist(energy, range=(0,20), bins = 100, label="all")
plt.hist(energy[num_clu<=2], range=(0,20), bins = 100, label="nClu=1&2")
plt.hist(energy[num_clu==1], range=(0,20), bins = 100, label="nClu=1")
plt.hist(energy[evt_fra>0.8], range=(0,20), bins = 100, label="EFrac>0.8")
plt.hist(energy[num_clu>2], range=(0,20), bins = 100, label="nClu>2")

plt.legend()


plt.figure("Num Clu")
plt.hist(num_clu, range=(0,20), bins = 20)
plt.hist(num_clu[cut_enerange], range=(0,20), bins = 20)
plt.hist(num_clu[numpy.logical_and(cut_enerange, my_cut)], range=(0,20), bins = 20)

plt.figure("TKR Board")
#plt.hist(trk_bord, range=(0,50), bins = 50)
plt.hist(trk_bord[cut_enerange], range=(0,50), bins = 50, label="Erange")
plt.hist(trk_bord[numpy.logical_and(cut_enerange, num_clu<=2)], range=(0,50), bins = 50, label="Erange&&nClu=1&2")
plt.legend()



plt.figure("EvtFra")
plt.hist(evt_fra, range=(0,1), bins = 51)
plt.hist(evt_fra[cut_enerange], range=(0,1), bins = 51)
plt.figure("EvtFra_vsE")
plt.hist2d(energy[cut_enerange], evt_fra[cut_enerange], range=((2,8),(0,1)), bins=(60, 100), norm=LogNorm())

plt.figure("PixFra_vsE")
plt.hist2d(energy[cut_enerange], pix_fra[cut_enerange], range=((2,8),(0,1)), bins=(60, 100), norm=LogNorm())

#plt.figure("EneDens")
#plt.hist(ene_dens[cut_enerange], range=(0,1), bins=100))
#plt.figure("EneDens_vsE")
#plt.hist2d(energy[cut_enerange], ene_dens[cut_enerange], range=((2,8),(0,1)), bins=(60, 100), norm=LogNorm())


plt.figure("Eve_trksize_all")
plt.hist2d(energy, trk_size, range=((0, 15),(0, 600)), bins=(30, 120), norm=LogNorm())
plt.colorbar()



plt.figure("RaDec")
plt.hist2d(ra[cut_enerange], dec[cut_enerange], bins=200, norm=LogNorm())
plt.colorbar()

#input()

plt.figure("RaDec_clu1")
plt.hist2d(ra[numpy.logical_and(cut_enerange, num_clu==1)],
           dec[numpy.logical_and(cut_enerange, num_clu==1)],
           bins=200, norm=LogNorm())
plt.colorbar()

plt.figure("RaDec_efrac1")
plt.hist2d(ra[numpy.logical_and(cut_enerange, evt_fra<0.7)],
           dec[numpy.logical_and(cut_enerange, evt_fra<0.7)],
           bins=200, norm=LogNorm())
plt.colorbar()

# write mask to out file
if False:
    mask_out = numpy.logical_and(cut_enerange, my_cut )
    with open("mask_nobkg_test.npy", 'wb') as f:
        numpy.save(f, mask_out)
    with open("mask_nobkg_test.npy", 'rb') as f:
        mask_in = numpy.load(f)
            
    print(sum(mask_out!=mask_in))
