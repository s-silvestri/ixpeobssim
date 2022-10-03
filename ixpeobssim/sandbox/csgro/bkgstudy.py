import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ixpeobssim.sandbox.csgro.friendFiles import xEventFileFriend

def hratio(h0, h1):
    ratio = h0/h1
    ratio_err = ratio*numpy.sqrt(1./h0 + 1./h1)
    return ratio, ratio_err
#
# For CenX3
#
commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CenX3/"
l1path = []
l2path = []
for DuId in [1,2,3]:
    l1path.append(commonpath+"event_l1/ixpe01006501_det%d_evt1_v02.fits" % DuId)
    l2path.append(commonpath+"event_l2/ixpe01006501_det%d_evt2_v02.fits" % DuId)

myCenX3 = xEventFileFriend(l2path, l1path)

tgr_id_l2 = myCenX3.l2value("TRG_ID")
tgr_id_l1 = myCenX3.l1value("TRG_ID")
print("Trigger Id diff check, must be zero, and the winner is:",
      sum(tgr_id_l1 - tgr_id_l2))

ra_CenX3, dec_CenX3  = myCenX3.sky_position_data()
energy_CenX3   = myCenX3.energy_data()
#num_clu_CenX3  = myCenX3.l1value("NUM_CLU")
#num_pix_CenX3  = myCenX3.l1value("NUM_PIX")
#trk_size_CenX3 = myCenX3.l1value("TRK_SIZE")
#trk_bord_CenX3 = myCenX3.l1value("TRK_BORD")
evt_fra_CenX3  = myCenX3.l1value("EVT_FRA")
du_status_CenX3 = myCenX3.l1value("DU_STATUS")
livetime_CenX3 = myCenX3.l1value("LIVETIME")
#w_mom_CenX3    = myCenX3.l2value("W_MOM")

#ene_dens_CenX3 = energy_CenX3/trk_size_CenX3
#pix_fra_CenX3  = trk_size_CenX3/num_pix_CenX3

cut_enerange_CenX3 = numpy.logical_and(energy_CenX3>2, energy_CenX3<8)
#my_cut_CenX3 = evt_fra_CenX3>0.7
my_cut_CenX3 = evt_fra_CenX3>((0.2/6.)*(energy_CenX3-2) + 0.65)
x = numpy.linspace(2,8,10)
y_cut_CenX3 = (0.15/6.)*(x-2) + 0.65

# plots
plt.ion()

plt.figure("EvtFravsE_CenX3")
plt.hist2d(energy_CenX3[cut_enerange_CenX3], evt_fra_CenX3[cut_enerange_CenX3],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CenX3 , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/EvtFravsE_CenX3.png")

plt.figure("EvtFravsE_CenX3_DGN")
plt.hist2d(energy_CenX3[numpy.logical_and(cut_enerange_CenX3, du_status_CenX3==2)],
           evt_fra_CenX3[numpy.logical_and(cut_enerange_CenX3, du_status_CenX3==2)],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CenX3 , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")


plt.figure("Energy_CenX3")
plt.subplot(211)
HEne_CenX3_all = plt.hist(energy_CenX3,
                          range=(0,20), bins = 100, label="all")
HEne_CenX3_cut = plt.hist(energy_CenX3[my_cut_CenX3],
                          range=(0,20), bins = 100, label="w/ cut")
plt.yscale('log')
plt.xlabel("Energy [keV]")
plt.xlim(0,10)
plt.legend()
plt.subplot(212)
HEne_CenX3_ratio, HEne_CenX3_ratio_err  = \
    hratio(HEne_CenX3_cut[0], HEne_CenX3_all[0])
HEne_CenX3_x = numpy.linspace(0+0.1,20-0.1, 100)
plt.errorbar(HEne_CenX3_x, HEne_CenX3_ratio, HEne_CenX3_ratio_err, xerr=0.1,
             fmt=".")
plt.xlabel("Energy [keV]")
plt.ylabel("Evt fractionEnergy")
plt.xlim(0,10)
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/Energy_CenX3.png")

print("Global event fraction after cut: %d/%d = %f" % \
      (len(energy_CenX3[my_cut_CenX3]), len(energy_CenX3),
       len(energy_CenX3[my_cut_CenX3])/len(energy_CenX3) ))

print("Global event fraction in [2,8] kev range after cut: %d/%d = %f" % \
      (len(energy_CenX3[numpy.logical_and(cut_enerange_CenX3, my_cut_CenX3)]),
           len(energy_CenX3[cut_enerange_CenX3]),
           len(energy_CenX3[numpy.logical_and(cut_enerange_CenX3, my_cut_CenX3)])/len(energy_CenX3[cut_enerange_CenX3]) ))



plt.figure("RaDec_SgnCut_CenX3")
plt.hist2d(ra_CenX3[numpy.logical_and(cut_enerange_CenX3, my_cut_CenX3)],
           dec_CenX3[numpy.logical_and(cut_enerange_CenX3, my_cut_CenX3)],
           bins=100, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/RaDec_SgnCut_CenX3.png")

plt.figure("RaDec_BkgCut_CenX3")
plt.hist2d(ra_CenX3[numpy.logical_and(cut_enerange_CenX3,
                                      numpy.logical_not(my_cut_CenX3))],
           dec_CenX3[numpy.logical_and(cut_enerange_CenX3,
                                       numpy.logical_not(my_cut_CenX3))],
           bins=100, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/RaDec_BkgCut_CenX3.png")

plt.figure("Livetime")
plt.hist(livetime_CenX3, bins=100)

#
# For Cas A
#
commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CasA/obs01001301/"
l1path = []
l2path = []
for DuId in [1]:
    l1path.append(commonpath+"event_l1/ixpe01001301_det%d_evt1_v01.fits" % DuId)
    l2path.append(commonpath+"event2_corr/ixpe01001301E_det%d_evt2_v02_picorr.fits" % DuId)

myCasA = xEventFileFriend(l2path, l1path)

tgr_id_l2 = myCasA.l2value("TRG_ID")
tgr_id_l1 = myCasA.l1value("TRG_ID")
print("Trigger Id diff check, must be zero, and the winner is:",
      sum(tgr_id_l1 - tgr_id_l2))

ra_CasA, dec_CasA  = myCasA.sky_position_data()
energy_CasA   = myCasA.energy_data()
#num_clu_CasA  = myCasA.l1value("NUM_CLU")
num_pix_CasA  = myCasA.l1value("NUM_PIX")
trk_size_CasA = myCasA.l1value("TRK_SIZE")
#trk_bord_CasA = myCasA.l1value("TRK_BORD")
evt_fra_CasA  = myCasA.l1value("EVT_FRA")
#w_mom_CasA    = myCasA.l2value("W_MOM")
du_status_CasA = myCasA.l1value("DU_STATUS")
livetime_CasA = myCasA.l1value("LIVETIME")
livetime_all_CasA = myCasA.l1value("LIVETIME", allevts = True)
#time_all_CasA = myCasA.l1value("TIME", allevts = True)
#time_min_CasA = round(numpy.min(time_all_CasA))
#time_max_CasA = round(numpy.max(time_all_CasA))
#time_l2_CasA = myCasA.l1value("TIME")
status2_all_CasA = myCasA.l1value("STATUS2", allevts = True)
evt_fra_all_CasA = myCasA.l1value("EVT_FRA", allevts = True)
pha_all_CasA     = myCasA.l1value("PHA_EQ", allevts = True)

ene_dens_CasA = energy_CasA/trk_size_CasA
pix_fra_CasA  = trk_size_CasA/num_pix_CasA

cut_enerange_CasA = numpy.logical_and(energy_CasA>2, energy_CasA<8)
#my_cut_CasA = evt_fra_CasA>0.7
my_cut_CasA = evt_fra_CasA>((0.2/6.)*(energy_CasA-2) + 0.65)
my_cut_CasA =numpy.logical_and(my_cut_CasA, livetime_CasA>15)
x = numpy.linspace(2,8,10)
y_cut_CasA = (0.15/6.)*(x-2) + 0.65

# WRITE MASK ON FILE
if False:
    mask_out = numpy.logical_and(cut_enerange_CasA, my_cut_CasA )
    with open("mask_nobkg_du3_CasA.npy", 'wb') as f:
        numpy.save(f, mask_out)


# plots
plt.ion()

plt.figure("PixFravsE_CasA")
plt.hist2d(energy_CasA[cut_enerange_CasA], pix_fra_CasA[cut_enerange_CasA],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Pixels in the main track")


plt.figure("EvtFravsE_CasA")
plt.hist2d(energy_CasA[cut_enerange_CasA], evt_fra_CasA[cut_enerange_CasA],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/EvtFravsE_CasA.png")


plt.figure("EvtFravsE_CasA_DGN")
plt.hist2d(energy_CasA[numpy.logical_and(cut_enerange_CasA, du_status_CasA==2)],
           evt_fra_CasA[numpy.logical_and(cut_enerange_CasA, du_status_CasA==2)],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")

plt.figure("EvtFravsE_CasA_Livetime")
plt.hist2d(energy_CasA[numpy.logical_and(cut_enerange_CasA, livetime_CasA<15)],
           evt_fra_CasA[numpy.logical_and(cut_enerange_CasA, livetime_CasA<15)],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")


plt.figure("LivetimevsE_CasA")
plt.hist2d(energy_CasA[cut_enerange_CasA], livetime_CasA[cut_enerange_CasA],
           range=((2,8),(0,100)), bins=(60, 100), norm=LogNorm())
plt.xlabel("Energy [keV]")
plt.ylabel("Livetime")

plt.figure("Energy_CasA")
plt.subplot(211)
HEne_CasA_all = plt.hist(energy_CasA,
                          range=(0,20), bins = 100, label="all")
HEne_CasA_cut = plt.hist(energy_CasA[my_cut_CasA],
                          range=(0,20), bins = 100, label="w/ cut")
plt.yscale('log')
plt.xlabel("Energy [keV]")
plt.xlim(0,10)
plt.legend()
plt.subplot(212)
HEne_CasA_ratio, HEne_CasA_ratio_err  = \
    hratio(HEne_CasA_cut[0], HEne_CasA_all[0])
HEne_CasA_x = numpy.linspace(0+0.1,20-0.1, 100)
plt.errorbar(HEne_CasA_x, HEne_CasA_ratio, HEne_CasA_ratio_err, xerr=0.1,
             fmt=".")
plt.xlabel("Energy [keV]")
plt.ylabel("Evt fractionEnergy")
plt.xlim(0,10)
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/Energy_CasA.png")

print("Global event fraction after cut: %d/%d = %f" % \
      (len(energy_CasA[my_cut_CasA]), len(energy_CasA),
       len(energy_CasA[my_cut_CasA])/len(energy_CasA) ))

print("Global event fraction in [2,8] kev range after cut: %d/%d = %f" % \
      (len(energy_CasA[numpy.logical_and(cut_enerange_CasA, my_cut_CasA)]),
           len(energy_CasA[cut_enerange_CasA]),
           len(energy_CasA[numpy.logical_and(cut_enerange_CasA, my_cut_CasA)])/len(energy_CasA[cut_enerange_CasA]) ))



plt.figure("RaDec_SgnCut_CasA")
plt.hist2d(ra_CasA[numpy.logical_and(cut_enerange_CasA, my_cut_CasA)],
           dec_CasA[numpy.logical_and(cut_enerange_CasA, my_cut_CasA)],
           bins=100, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/RaDec_SgnCut_CasA.png")

plt.figure("RaDec_BkgCut_CasA")
plt.hist2d(ra_CasA[numpy.logical_and(cut_enerange_CasA,
                                      numpy.logical_not(my_cut_CasA))],
           dec_CasA[numpy.logical_and(cut_enerange_CasA,
                                       numpy.logical_not(my_cut_CasA))],
           bins=100, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
#plt.savefig("/home/carmelo/xpe/xpework/CasA/bkgstudy/plots/RaDec_BkgCut_CasA.png")

#plt.figure("TimeHist")
#time_bin_size = 2*60 # seconds
#time_bins = round((time_max_CasA-time_min_CasA)/time_bin_size) +1
#print("{{{{{\nTimeHist with %d bins\n}}}}}" % time_bins)
#plt.hist(time_all_CasA-time_min_CasA, bins=time_bins)


# CalD, CalC 
# Energy calibration by hand
# DU1 peak 4755 & 17750
myenergy_all_CasA = (5.9-1.5)/(17750.-4755.)*(pha_all_CasA - 4755.) + 1.5
plt.figure("EvtFravsPha_CasA_CalD")
plt.hist2d(pha_all_CasA[status2_all_CasA[:,7]],
           evt_fra_all_CasA[status2_all_CasA[:,7]],
           range=((0,40000),(0,1.5)), bins=(80, 150), norm=LogNorm())
plt.xlabel("PHA [adc]")
plt.ylabel("Fraction of Energy in the main track")

plt.figure("EvtFravsPha_CasA_CalC")
plt.hist2d(pha_all_CasA[status2_all_CasA[:,6]],
           evt_fra_all_CasA[status2_all_CasA[:,6]],
           range=((0,40000),(0,1.5)), bins=(80, 150), norm=LogNorm())
plt.xlabel("PHA [adc]")
plt.ylabel("Fraction of Energy in the main track")

plt.figure("EvtFravsE_CasA_CalD")
plt.hist2d(myenergy_all_CasA[numpy.logical_and(status2_all_CasA[:,7], livetime_all_CasA>15)],
           evt_fra_all_CasA[numpy.logical_and(status2_all_CasA[:,7], livetime_all_CasA>15)],
           range=((0,12),(0,1.5)), bins=(80, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "b") # blue h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")

plt.figure("EvtFravsE_CasA_CalC")
plt.hist2d(myenergy_all_CasA[numpy.logical_and(status2_all_CasA[:,6], livetime_all_CasA>15)],
           evt_fra_all_CasA[numpy.logical_and(status2_all_CasA[:,6], livetime_all_CasA>15)],
           range=((0,12),(0,1.5)), bins=(80, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "b") # blue h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")


# & Occultation TBD

#
# For Cas A SIMULATION (No POL 890k secons)
#

#ixpesim --src-photon-list  casa_ph_seed1_du1_photon_list.fits  --dme-pressure 645. --output-file casa_ph_seed1_du1_ixpesim_cfg18_off1_gmap.fits -n 10000000
#--diagnostic-prescale 100
#--zero-sup-threshold 18
#--calib-xpol-gain /run/media/carmelo/TOSHIBAEXT/ixpedata/108/gainmaps/108_0000134_data_gainmap_gpd31.fits 


#commonpath = "/run/media/carmelo/TOSHIBAEXT/obssim_sim/casA_sim_nopol_890k/"
commonpath = "/run/media/carmelo/TOSHIBAEXT/obssim_sim/casA_sim_nopol_890k_config18/"
l2path = []
for DuId in [1]:
    for seedId in [1]:
        #l2path.append(commonpath+"casa_nopol_seed%d_du%d_ixpesim_recon_simfmt.fits" % (seedId, DuId))
        l2path.append(commonpath+"casa_ph_seed%d_du%d_ixpesim_cfg18_off1_gmap_recon_simfmt.fits" % (seedId, DuId))


# Level1 is useless in this case - use None
myCasASim = xEventFileFriend(l2path, filel1=None) 

ra_CasASim, dec_CasASim  = myCasASim.sky_position_data()
energy_CasASim   = myCasASim.energy_data()
evt_fra_CasASim  = myCasASim.l2value("EVT_FRA")
du_status_CasASim  = myCasASim.l2value("DU_STATUS")

cut_enerange_CasASim = numpy.logical_and(energy_CasASim>2, energy_CasASim<8)
#my_cut_CasA = evt_fra_CasA>0.7
my_cut_CasASim = evt_fra_CasASim>((0.2/6.)*(energy_CasASim-2) + 0.65)

# plots
plt.ion()

plt.figure("EvtFravsE_CasASim")
plt.hist2d(energy_CasASim[cut_enerange_CasASim],
           evt_fra_CasASim[cut_enerange_CasASim],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")

plt.figure("EvtFravsE_CasASim_DGN")
plt.hist2d(energy_CasASim[numpy.logical_and(cut_enerange_CasASim, du_status_CasASim==2)],
           evt_fra_CasASim[numpy.logical_and(cut_enerange_CasASim, du_status_CasASim==2)],
           range=((2,8),(0,1.5)), bins=(60, 150), norm=LogNorm())
plt.plot(x, y_cut_CasA , "k") # black h line 
plt.xlabel("Energy [keV]")
plt.ylabel("Fraction of Energy in the main track")

plt.figure("RaDec_SgnCut_CasASim")
plt.hist2d(ra_CasASim[numpy.logical_and(cut_enerange_CasASim, my_cut_CasASim)],
           dec_CasASim[numpy.logical_and(cut_enerange_CasASim, my_cut_CasASim)],
           bins=100, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()

plt.figure("RaDec_BkgCut_CasASim")
plt.hist2d(ra_CasASim[numpy.logical_and(cut_enerange_CasASim,
                                      numpy.logical_not(my_cut_CasASim))],
           dec_CasASim[numpy.logical_and(cut_enerange_CasASim,
                                       numpy.logical_not(my_cut_CasASim))],
           bins=100, norm=LogNorm())
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
