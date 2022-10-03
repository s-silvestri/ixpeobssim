import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ixpeobssim.sandbox.csgro.friendFiles import xEventFileFriend

#
# For Cas A
#

commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CasA/obs01001301/"

myDU1 = xEventFileFriend(commonpath+"event2_corr/ixpe01001301E_det1_evt2_v02_picorr.fits",
                           commonpath+"event_l1/ixpe01001301_det1_evt1_v01.fits")
myDU2 = xEventFileFriend(commonpath+"event2_corr/ixpe01001301E_det2_evt2_v02_picorr.fits",
                           commonpath+"event_l1/ixpe01001301_det2_evt1_v01.fits")
myDU3 = xEventFileFriend(commonpath+"event2_corr/ixpe01001301E_det3_evt2_v02_picorr.fits",
                           commonpath+"event_l1/ixpe01001301_det3_evt1_v01.fits")

#
# For CenX3
#
"""
commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CenX3/"
myDU1 = xEventFileFriend(commonpath+"event_l2/ixpe01006501_det1_evt2_v02.fits",
                           commonpath+"event_l1/ixpe01006501_det1_evt1_v02.fits")
myDU2 = xEventFileFriend(commonpath+"event_l2/ixpe01006501_det2_evt2_v02.fits",
                           commonpath+"event_l1/ixpe01006501_det2_evt1_v02.fits")
myDU3 = xEventFileFriend(commonpath+"event_l2/ixpe01006501_det3_evt2_v02.fits",
                           commonpath+"event_l1/ixpe01006501_det3_evt1_v02.fits")
"""
#
# For 4U
#
"""
commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/4U_0142+65/ix01003201/01003201/"
myDU1 = xEventFileFriend(commonpath+"event_l2/ixpe01003201_det1_evt2_v01.fits",
                           commonpath+"event_l1/ixpe01003201_det1_evt1_v01.fits")
myDU2 = xEventFileFriend(commonpath+"event_l2/ixpe01003201_det2_evt2_v01.fits",
                           commonpath+"event_l1/ixpe01003201_det2_evt1_v01.fits")
myDU3 = xEventFileFriend(commonpath+"event_l2/ixpe01003201_det3_evt2_v01.fits",
                           commonpath+"event_l1/ixpe01003201_det3_evt1_v01.fits")
"""

#
# Cas A sim
#
"""
commonpath = "/run/media/carmelo/TOSHIBAEXT/obssim_sim/casA_sim_nopol_890k/"
myDU1 = xEventFileFriend(commonpath+"casa_nopol_seed1_du1_ixpesim_recon_simfmt.fits",
                         commonpath+"casa_nopol_seed1_du1_ixpesim_recon_simfmt.fits")
myDU2 = xEventFileFriend(commonpath+"casa_nopol_seed1_du2_ixpesim_recon_simfmt.fits",
                         commonpath+"casa_nopol_seed1_du2_ixpesim_recon_simfmt.fits")
myDU3 = xEventFileFriend(commonpath+"casa_nopol_seed1_du3_ixpesim_recon_simfmt.fits",
                         commonpath+"casa_nopol_seed1_du3_ixpesim_recon_simfmt.fits")
"""

ra_DU1, dec_DU1  = myDU1.sky_position_data()
energy_DU1   = myDU1.energy_data()
detphi_DU1 = myDU1.l1value("DETPHI2") +numpy.deg2rad(229)#+numpy.deg2rad(109)
q_DU1      = myDU1.l2value("Q")
u_DU1      = myDU1.l2value("U")
cut_DU1 = energy_DU1>5
phi_DU1 = 0.5*numpy.arctan2(u_DU1,q_DU1)
detq_DU1 = 2. * numpy.cos(2. * detphi_DU1)
detu_DU1 = 2. * numpy.sin(2. * detphi_DU1)

ra_DU2, dec_DU2  = myDU2.sky_position_data()
energy_DU2   = myDU2.energy_data()
detphi_DU2 = myDU2.l1value("DETPHI2") +numpy.deg2rad(109)#+numpy.deg2rad(229)
q_DU2      = myDU2.l2value("Q")
u_DU2      = myDU2.l2value("U")
cut_DU2 = energy_DU2>5
phi_DU2 = 0.5*numpy.arctan2(u_DU2,q_DU2)
detq_DU2 = 2. * numpy.cos(2. * detphi_DU2)
detu_DU2 = 2. * numpy.sin(2. * detphi_DU2)

ra_DU3, dec_DU3  = myDU3.sky_position_data()
energy_DU3   = myDU3.energy_data()
detphi_DU3 = myDU3.l1value("DETPHI2") +numpy.deg2rad(349-180)#+numpy.deg2rad(349-180)
q_DU3      = myDU3.l2value("Q")
u_DU3      = myDU3.l2value("U")
cut_DU3 = energy_DU3>5
phi_DU3 = 0.5*numpy.arctan2(u_DU3,q_DU3)
detq_DU3 = 2. * numpy.cos(2. * detphi_DU3)
detu_DU3 = 2. * numpy.sin(2. * detphi_DU3)

clk_cor_DU1 = numpy.deg2rad(120)
clk_cor_DU2 = 0.
clk_cor_DU3 = numpy.deg2rad(240)
detq_cor_DU1 = 2. * numpy.cos(2. * (detphi_DU1 + clk_cor_DU1))
detu_cor_DU1 = 2. * numpy.sin(2. * (detphi_DU1 + clk_cor_DU1))
detq_cor_DU2 = 2. * numpy.cos(2. * (detphi_DU2 + clk_cor_DU2))
detu_cor_DU2 = 2. * numpy.sin(2. * (detphi_DU2 + clk_cor_DU2))
detq_cor_DU3 = 2. * numpy.cos(2. * (detphi_DU3 + clk_cor_DU3))
detu_cor_DU3 = 2. * numpy.sin(2. * (detphi_DU3 + clk_cor_DU3))



# plots
plt.ion()

plt.figure("Phi_DU", figsize=(10,6), tight_layout=True)
# DU1
plt.subplot(231)
hh = plt.hist2d(detphi_DU1, phi_DU1, bins=100)
plt.title("DU 1")
plt.xlabel("DET Phi")
plt.ylabel("Lv2 Phi")
plt.subplot(234)
hh = plt.hist(numpy.rad2deg(detphi_DU1 - phi_DU1), bins=100)
plt.xlabel("DET Phi - Lv2 Phi [deg]")
# DU2
plt.subplot(232)
hh = plt.hist2d(detphi_DU2, phi_DU2, bins=100)
plt.title("DU 2")
plt.xlabel("DET Phi")
plt.ylabel("Lv2 Phi")
plt.subplot(235)
hh = plt.hist(numpy.rad2deg(detphi_DU2 - phi_DU2), bins=100)
plt.xlabel("DET Phi - Lv2 Phi [deg]")
# DU3
plt.subplot(233)
hh = plt.hist2d(detphi_DU3, phi_DU3, bins=100)
plt.title("DU 3")
plt.xlabel("DET Phi")
plt.ylabel("Lv2 Phi")
plt.subplot(236)
hh = plt.hist(numpy.rad2deg(detphi_DU3 - phi_DU3), bins=100)
plt.xlabel("DET Phi - Lv2 Phi [deg]")

#########3
raise

plt.figure("QU_DU", figsize=(10,6), tight_layout=True)
plt.subplot(231)
hh = plt.hist2d(detq_DU1, q_DU1, bins=100, range=((-5,5), (-5, 5)))
plt.title("DU 1")
plt.xlabel("DET Q")
plt.ylabel("LV2 Q")
plt.subplot(234)
hh = plt.hist2d(detu_DU1, u_DU1, bins=100, range=((-5,5), (-5, 5)))
plt.xlabel("DET U")
plt.ylabel("LV2 U")

#plt.figure("QU_DU_DU2")
plt.subplot(232)
hh = plt.hist2d(detq_DU2, q_DU2, bins=100, range=((-5,5), (-5, 5)))
plt.title("DU 2")
plt.xlabel("DET Q")
plt.ylabel("LV2 Q")
plt.subplot(235)
hh = plt.hist2d(detu_DU2, u_DU2, bins=100, range=((-5,5), (-5, 5)))
plt.xlabel("DET U")
plt.ylabel("LV2 U")

plt.subplot(233)
hh = plt.hist2d(detq_DU3, q_DU3, bins=100, range=((-5,5), (-5, 5)))
plt.title("DU 3")
plt.xlabel("DET Q")
plt.ylabel("LV2 Q")
plt.subplot(236)
hh = plt.hist2d(detu_DU3, u_DU3, bins=100, range=((-5,5), (-5, 5)))
plt.xlabel("DET U")
plt.ylabel("LV2 U")


plt.figure("QU_DU_Corr", figsize=(10,6), tight_layout=True)
# DU1
plt.subplot(231)
hh = plt.hist2d(detq_cor_DU1, q_DU1, bins=100, range=((-5,5), (-5, 5)))
plt.title("DU 1 (phi+%.1f)" % numpy.rad2deg(clk_cor_DU1))
plt.xlabel("DET Q")
plt.ylabel("LV2 Q")
plt.subplot(234)
hh = plt.hist2d(detu_cor_DU1, u_DU1, bins=100, range=((-5,5), (-5, 5)))
plt.xlabel("DET U")
plt.ylabel("LV2 U")
# DU2
plt.subplot(232)
hh = plt.hist2d(detq_cor_DU2, q_DU2, bins=100, range=((-5,5), (-5, 5)))
plt.title("DU 2 (phi+%.1f)" % numpy.rad2deg(clk_cor_DU2))
plt.xlabel("DET Q")
plt.ylabel("LV2 Q")
plt.subplot(235)
hh = plt.hist2d(detu_cor_DU2, u_DU2, bins=100, range=((-5,5), (-5, 5)))
plt.xlabel("DET U")
plt.ylabel("LV2 U")
# DU3
plt.subplot(233)
hh = plt.hist2d(detq_cor_DU3, q_DU3, bins=100, range=((-5,5), (-5, 5)))
plt.title("DU 3 (phi+%.1f)" % numpy.rad2deg(clk_cor_DU3))
plt.xlabel("DET Q")
plt.ylabel("LV2 Q")
plt.subplot(236)
hh = plt.hist2d(detu_cor_DU3, u_DU3, bins=100, range=((-5,5), (-5, 5)))#, norm=LogNorm())
plt.xlabel("DET U")
plt.ylabel("LV2 U")
