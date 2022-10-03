import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ixpeobssim.sandbox.csgro.friendFiles import xEventFileFriend

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
#num_pix_CasA  = myCasA.l1value("NUM_PIX")
#trk_size_CasA = myCasA.l1value("TRK_SIZE")
#trk_bord_CasA = myCasA.l1value("TRK_BORD")
evt_fra_CasA  = myCasA.l1value("EVT_FRA")
#w_mom_CasA    = myCasA.l2value("W_MOM")
du_status_CasA = myCasA.l1value("DU_STATUS")
livetime_CasA = myCasA.l1value("LIVETIME")
time_all_CasA = myCasA.l1value("TIME", allevts = True)
time_min_CasA = round(numpy.min(time_all_CasA))
time_max_CasA = round(numpy.max(time_all_CasA))
time_l2_CasA = myCasA.l1value("TIME")
status2_CasA = myCasA.l1value("STATUS2", allevts = True)
pha_CasA     = myCasA.l1value("PHA_EQ", allevts = True)
status2_l2_CasA = myCasA.l1value("STATUS2")

cut_enerange_CasA = numpy.logical_and(energy_CasA>2, energy_CasA<8)
#my_cut_CasA = evt_fra_CasA>0.7
my_cut_CasA = evt_fra_CasA>((0.2/6.)*(energy_CasA-2) + 0.65)
my_cut_CasA =numpy.logical_and(my_cut_CasA, livetime_CasA>15)
#x = numpy.linspace(2,8,10)
#y_cut_CasA = (0.15/6.)*(x-2) + 0.65
occultation_cut = numpy.logical_or(status2_CasA[:,10], status2_CasA[:,11])

#print("N evt %d" % len(time_all_CasA))
#print("N evt in occultation %d" % sum(occultation_cut))

#occultation_cut = occultation_cut[:100001]

import time
t0 = time.time()

transitions = numpy.diff(occultation_cut)
edges       = numpy.where(transitions==True)[0]
occ_start   = edges[::2] +1
occ_stop    = edges[1::2]
n = len(occ_start)
print("Found %d, %d transitions in %f s" %\
      (n, len(occ_stop), time.time() -t0))

good_occ = 0
occ_duration_list = []
occultation_goodtime_cut = time_all_CasA < 0 # all false
for i in []:#range(n):
    print(i,
          time_all_CasA[occ_start[i]]-time_min_CasA,
          time_all_CasA[occ_stop[i]]-time_min_CasA)
    occ_duration = time_all_CasA[occ_stop[i]] - time_all_CasA[occ_start[i]]
    occ_duration_list.append(occ_duration)
    if occ_duration>900:
        # min 15 minutes of duration
        # take central 5 minutes
        occ_avg_time = 0.5*(time_all_CasA[occ_stop[i]] + time_all_CasA[occ_start[i]])
        occ_min_time = occ_avg_time-2.5*60.
        occ_max_time = occ_avg_time+2.5*60.
        good_occ+=1
        tcut = numpy.logical_and(time_all_CasA > occ_min_time,
                                 time_all_CasA < occ_max_time)
        occultation_goodtime_cut = numpy.logical_or(occultation_goodtime_cut,
                                                    tcut)
print("found %d good occultations" % good_occ)

# plots
plt.ion()

plt.figure()
plt.subplot(211)
plt.hist(pha_CasA[status2_CasA[:,6]], range=(0,40000), bins=100)
plt.title("CalC (CalHI - bit 6) pha distribution")
#plt.xlabel("PHA [adc]")
plt.subplot(212)
plt.hist(pha_CasA[status2_CasA[:,7]], range=(0,40000), bins=100)
plt.title("CalD (CalLO - bit 7) pha distribution")
plt.xlabel("PHA [adc]")

raise

plt.figure()
plt.hist(occ_duration_list, bins=100)

plt.figure("TimeHist")
time_bin_size = 0.5*60 # seconds
time_bins = round((time_max_CasA-time_min_CasA)/time_bin_size) +1
plt.hist(time_all_CasA[occultation_cut]-time_min_CasA, bins=time_bins, label="occultation")
 
 
plt.figure("TimeHist1")
time_bin_size = 2*60 # seconds
time_bins = round((time_max_CasA-time_min_CasA)/time_bin_size) +1
plt.hist(time_all_CasA[occultation_goodtime_cut]-time_min_CasA, bins=time_bins, label="occultation")

"""
plt.figure("TimeHist")
time_bin_size = 2*60 # seconds
time_bins = round((time_max_CasA-time_min_CasA)/time_bin_size) +1
print("{{{{{\nTimeHist with %d bins\n}}}}}" % time_bins)
plt.hist(time_all_CasA-time_min_CasA, bins=time_bins, label="all")
plt.hist(time_all_CasA[status2_CasA[:,6]]-time_min_CasA, bins=time_bins, label="calC")
plt.hist(time_all_CasA[occultation_cut]-time_min_CasA, bins=time_bins, label="occultation")
plt.hist(time_all_CasA[status2_CasA[:,12]]-time_min_CasA, bins=time_bins, label="in SAA")

plt.legend()
plt.xlim(0, 10000)
"""

