import numpy
from ixpeobssim.evt.event import xEventFile, xEventFileFriend
from ixpeobssim.irfgen import NUM_CHANNELS, channel_to_energy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/Crab/"
l1path = [commonpath+"01001099/event_l1/ixpe01001001_det1_evt1_v04.fits",
          commonpath+"01001099/event_l1/ixpe01001002_det1_evt1_v03.fits"]
#l2path = [commonpath+"01001099/event_l2/ixpe01001099_det1_evt2_v02.fits"]
l2path = [commonpath+"01001099/event_l2/ixpe01001099_det1_evt2_v02_livet50.fits"]
# friendsing doesn't work on bary-corr files.
#commonpath+"Lev2Corrected/ixpe01001099_det1_evt2_v01_aspectcorr_picorr_wcscorr_barycorr_tlcorr.fits"]
DuId = 1

myFiles = xEventFileFriend(l2path, l1path)

tgr_id_l2 = myFiles.l2value("TRG_ID")
tgr_id_l1 = myFiles.l1value("TRG_ID")
#print(tgr_id_l1)
print("Trigger Id diff check, must be zero, and the winner is:",
          sum(tgr_id_l1 - tgr_id_l2))


energy   = myFiles.energy_data()
pha      = myFiles.l1value("PHA_EQ")
#num_clu  = myFiles.l1value("NUM_CLU")
#trk_size = myFiles.l1value("TRK_SIZE")
#trk_bord = myFiles.l1value("TRK_BORD")
evt_fra  = myFiles.l1value("EVT_FRA")
livetime = myFiles.l1value("LIVETIME")
ra, dec  = myFiles.sky_position_data()
absx     = myFiles.l1value("ABSX")
absy     = myFiles.l1value("ABSY")

livet_cut = livetime>50
ene_cut = numpy.logical_and(energy>2, energy<8)
#ene_cut = numpy.logical_and(ene_cut, livet_cut)
ra_cut  = numpy.logical_and(ra>83.615, ra<83.645)
dec_cut = numpy.logical_and(dec>22.005, dec<22.03)
roi_cut = numpy.logical_and(ra_cut, dec_cut)
roi_round_cut =  numpy.logical_and(roi_cut, ((ra-83.63304)**2 + (dec-22.01449)**2 < 0.003**2))
_cut = numpy.logical_and(ene_cut, roi_round_cut)

# WRITE MASK ON FILE
if False:
    mask_out = livetime>50
    with open("mask_livet50_du1_Crab.npy", 'wb') as f:
        numpy.save(f, mask_out)
        
# plots various stuff

plt.figure()
#plt.hist2d(energy, pha, bins = 300)
plt.hist(energy/pha, bins = 300)

plt.figure("Livetime")
plt.hist(livetime, bins=100, range=(0, 200))
plt.hist(livetime[ene_cut], bins=100, range=(0, 200))
plt.hist(livetime[_cut], bins=100, range=(0, 200))

plt.figure("Livetime2")
plt.hist2d(livetime, energy, bins = (100, 200), range = ((0, 200), (0,10)), norm=LogNorm())

plt.figure("Spectrum")
plt.subplot(211)
lvtthr = 50
h0 = plt.hist(energy, bins=150, range=(0, 15), label="all")
h1 = plt.hist(energy[livetime>lvtthr], bins=150, range=(0, 15), label="evtlivet>%d" % lvtthr)
plt.legend()
#plt.hist(energy[livetime>50], bins=200, range=(0, 20), label="evtlivet>50")
plt.subplot(212)
hx = 0.5*(h0[1][1:]+h0[1][:-1])
hxe = 0.5*(h0[1][1:]-h0[1][:-1])
plt.errorbar(hx, h1[0]/h0[0], xerr=hxe )
plt.grid()
plt.xlim(0,15)

print("Livetime cut in range 2-8 keV, selecting %d /%d evts." %
      (len(energy[numpy.logical_and(ene_cut, livet_cut)]), len(energy[ene_cut])))

plt.figure("RaDec_All")
plt.hist2d(ra[numpy.logical_and(ene_cut, livet_cut)],
           dec[numpy.logical_and(ene_cut, livet_cut)],
           bins=200, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
ax = plt.gca(); ax.invert_xaxis();

plt.figure("ABSXY_All")
plt.hist2d(absx[ene_cut],
           absy[ene_cut],
           bins=200, range=((-8,8), (-8,8)), norm=LogNorm())
plt.xlabel("ABSX")
plt.ylabel("ABSY")
plt.colorbar()

plt.figure("ABSXY_RegionCut")
plt.hist2d(absx[_cut],
           absy[_cut],
           bins=200, range=((-8,8), (-8,8)), norm=LogNorm())
plt.xlabel("ABSX")
plt.ylabel("ABSY")
plt.colorbar()

plt.show()
