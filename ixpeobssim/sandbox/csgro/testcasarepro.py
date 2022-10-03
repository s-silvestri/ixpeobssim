import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ixpeobssim.sandbox.csgro.friendFiles import xEventFileFriend
from ixpeobssim.evt.event import xEventFile

#
# For Cas A
#
commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CasA/obs01001301/"
"""
l1path = []
l2path = []
for DuId in [1]:
    l1path.append(commonpath+"event_l1/ixpe01001301_det%d_evt1_v01.fits" % DuId)
    #l1path.append(commonpath+"event_l2_repro/ixpe01001301_det%d_evt2_v02.fits" % DuId)
    l2path.append(commonpath+"event2_corr/ixpe01001301E_det%d_evt2_v02_picorr.fits" % DuId)

myCasA = xEventFileFriend(l2path, l1path)

tgr_id_l2 = myCasA.l2value("TRG_ID")
tgr_id_l1 = myCasA.l1value("TRG_ID")
print("Trigger Id diff check, must be zero, and the winner is:",
      sum(tgr_id_l1 - tgr_id_l2))
"""

DUID = 2

# UNCORRECTED
"""
filel2    = xEventFile(commonpath+
                       "event_l2/ixpe01001301_det%d_evt2_v01.fits" %DUID)
filel2rep = xEventFile(commonpath+
                       "event_l2_repro/ixpe01001301_det%d_evt2_v02.fits" %DUID)
"""
# PI-CORRECTED FILES

#filel2    = xEventFile(commonpath+
#                       "event_l2/ixpe01001301_det%d_evt2_v01_picorr.fits" %DUID)
#filel2    = xEventFile(commonpath+
#                       "event2_corr/ixpe01001301E_det%d_evt2_v02_picorr.fits" %DUID)
filel2 = xEventFile(commonpath+
                       "event_l2_repro/ixpe01001301_det%d_evt2_v02_select_picorr.fits" %DUID)
filel2rep = xEventFile(commonpath+
                       "event_l2_repro/ixpe01001301_det%d_evt2_v02_select_picorr1.fits" %DUID)


timel2 = filel2.time_data()
timel2r = filel2rep.time_data()
time_min = min(min(timel2), min(timel2r))
time_max = max(max(timel2), max(timel2r))

trgidl2  =  filel2.event_data["TRG_ID"]
trgidl2r =  filel2rep.event_data["TRG_ID"]

pil2  =  filel2.event_data["PI"]
pil2r =  filel2rep.event_data["PI"]

enel2 = filel2.energy_data()
enel2r = filel2rep.energy_data()

q = filel2.event_data["Q"]
qr = filel2rep.event_data["Q"]
u = filel2.event_data["U"]
ur = filel2rep.event_data["U"]

print("DU %d" % DUID)
print("Total numer of events: org %d, repro %d" % (len(timel2), len(timel2r)))

# GTI data
startl2 = filel2.gti_data()['START']
stopl2  = filel2.gti_data()['STOP']
startl2r = filel2rep.gti_data()['START']
stopl2r  = filel2rep.gti_data()['STOP']

gtil2 = stopl2-startl2
gtil2r = stopl2r-startl2r
print("GIT max/min:", max(gtil2r-gtil2), min(gtil2r-gtil2))
print("START max/min:", max(startl2r-startl2), min(startl2r-startl2))
print("STOP max/min:", max(stopl2r-stopl2), min(stopl2r-stopl2))

evalGTImask = False
if evalGTImask:
    print("Eval GTI mask for org files")
    gti_cut_l2 = numpy.array([False]* len(timel2))
    for i in range(len(stopl2)):
        gti_cut_l2 = numpy.logical_or(gti_cut_l2,
                                      numpy.logical_and(timel2>=startl2[i],
                                                        timel2<=stopl2[i]))
        
    print("Mask done, evts in GTI: %d/%d " % (sum(gti_cut_l2), len(gti_cut_l2)))

    print("Eval GTI mask for repro files")
    gti_cut_l2r = numpy.array([False]* len(timel2r))
    for i in range(len(stopl2r)):
        gti_cut_l2r = numpy.logical_or(gti_cut_l2r,
                                       numpy.logical_and(timel2r>=startl2r[i],
                                                         timel2r<=stopl2r[i]))

    print("Mask done, evts in GTI: %d/%d " % (sum(gti_cut_l2r), len(gti_cut_l2r)))

    # WRITE MASK ON FILE
    with open("mask_inGTI_l2_CasA_DU%d.npy" % DUID, 'wb') as f:
        numpy.save(f, gti_cut_l2 )
    with open("mask_inGTI_l2rep_CasA_DU%d.npy" % DUID, 'wb') as f:
        numpy.save(f, gti_cut_l2r )
else:
    # LOAD MASK from file:
    with open("mask_inGTI_l2_CasA_DU%d.npy" % DUID, 'rb') as f:
        gti_cut_l2 = numpy.load(f)
    print("Mask loaded, evts in GTI: %d/%d " % (sum(gti_cut_l2), len(gti_cut_l2)))
    with open("mask_inGTI_l2rep_CasA_DU%d.npy" % DUID, 'rb') as f:
        gti_cut_l2r = numpy.load(f)
    print("Mask loaded, evts in GTI: %d/%d " % (sum(gti_cut_l2r), len(gti_cut_l2r)))    
"""
# how many events in l2rep are in l2?
print("how many events in l2rep are in l2?")
t2r_IN_t2 =numpy.in1d(timel2r, timel2)  
print(sum(t2r_IN_t2), "out of", len(t2r_IN_t2))

# how many events in l2 are in l2rep?
print("how many events in l2 are in l2rep?")
t2_IN_t2r =numpy.in1d(timel2, timel2r)  
print(sum(t2_IN_t2r), "out of", len(t2_IN_t2r))
"""

commonevt = numpy.intersect1d(timel2, timel2r, assume_unique=True)# return_indices=False)[source]
print("Num of common events %d" % len(commonevt))

commonevtInGTI = numpy.intersect1d(timel2[gti_cut_l2], timel2r[gti_cut_l2r], assume_unique=True)# return_indices=False)[source]
print("Num of common eventsIn GTI %d" % len(commonevtInGTI))

mask_l2  = numpy.in1d(timel2,  commonevt)
mask_l2r = numpy.in1d(timel2r, commonevt)
# ONLY EVT IN GTI - NOT IN COMON
mask_nc_l2  = numpy.logical_and(numpy.logical_not(mask_l2),  gti_cut_l2) 
mask_nc_l2r = numpy.logical_and(numpy.logical_not(mask_l2r), gti_cut_l2r) 
mask_l2  = numpy.logical_and(mask_l2,  gti_cut_l2) # ONLY EVT IN GTI
mask_l2r = numpy.logical_and(mask_l2r, gti_cut_l2r) # ONLY EVT IN GTI
print("mask len: ", len(mask_l2), len(mask_l2r))

plt.ion()

# PLOT time of evt not in GIT
time_bin_size = 5*60#
time_bin_min  = round(time_min-1)
time_bin_max  = round(time_max+1)
time_bin_n    = round((time_bin_max -  time_bin_min)/time_bin_size +1)
plt.figure("EVT IN GTI")
plt.subplot(211)
hh = plt.hist(timel2[numpy.logical_not(gti_cut_l2)] - time_bin_min , bins =time_bin_n,  range = (0, time_bin_max-time_bin_min), label='outof GIT original')
plt.xlabel("Time [s]")
plt.legend()
plt.subplot(212)
hh = plt.hist(timel2r[numpy.logical_not(gti_cut_l2r)] - time_bin_min , bins =time_bin_n,  range = (0, time_bin_max-time_bin_min), label='outof GIT repro')
plt.xlabel("Time [s]")
plt.legend()
plt.yscale("log")

plt.figure() # EVT NOT IN COMMON
plt.subplot(211)
hh = plt.hist(timel2[mask_nc_l2] - time_bin_min , bins =time_bin_n,  range = (0, time_bin_max-time_bin_min), label='NOT in common; original')
plt.xlabel("Time [s]")
plt.legend()
plt.subplot(212)
hh = plt.hist(timel2r[mask_nc_l2r] - time_bin_min , bins =time_bin_n,  range = (0, time_bin_max-time_bin_min), label='NOT in common;  repro')
plt.xlabel("Time [s]")
plt.legend()


plt.figure()
hh = plt.hist2d(pil2[mask_l2], pil2r[mask_l2r], bins=400, range=((0,400),(0,400)) , norm=LogNorm())

plt.figure()
plt.subplot(211)
#hh = plt.hist(pil2r[mask_l2r]-pil2[mask_l2], bins=200, range=(-100,100))
hh = plt.hist(pil2r[gti_cut_l2r]-pil2[gti_cut_l2], bins=200, range=(-100,100))
plt.xlabel("PI_reprocess - PI")
plt.yscale('log')
plt.subplot(212)
#hh = plt.hist2d(pil2[mask_l2], pil2r[mask_l2r]-pil2[mask_l2], bins=(400,200), range=((0,400),(-100,100)) , norm=LogNorm())
hh = plt.hist2d(pil2[gti_cut_l2], pil2r[gti_cut_l2r]-pil2[gti_cut_l2], bins=(400,200), range=((0,400),(-100,100)) , norm=LogNorm())
plt.ylim(-5,5)
plt.ylabel("PI_reprocess - PI")
plt.xlabel("PI")

print("Nevt (in common) with different PI %d " %
      sum(pil2r[mask_l2r]-pil2[mask_l2] != 0) )
