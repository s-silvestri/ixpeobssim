# just my test script to make a few plot for calib source degub, after running:
# xpobssim.py  --configfile config/toy_point_source.py --duration 50000  --occult True --onbrdcalib True and several options...

from astropy.io import fits
import numpy
from ixpeobssim.utils.matplotlib_ import plt

#plt.figure("SrcTime", figsize=(12,10), tight_layout=True)
#plt.figure("CalCSpectrum", figsize=(8,10))
#plt.figure("CalCPhi", figsize=(8,10))
#plt.figure("CalCDetXY", figsize=(6,10))
#plt.figure("CalCDeltaT", figsize=(8,10))

DU_LIST = [1]
for DU in DU_LIST:

    file_path = "/home/carmelo/xpe/ixpeobssimdata/testCalC/simpleSimCalCOpt/toy_point_source_du%d.fits" % DU

    hdu = fits.open(file_path)
    gti_data = hdu['GTI'].data
    mc_data  =  hdu['MONTE_CARLO'].data
    evt_data =  hdu['EVENTS'].data
    
    srcid    = mc_data['SRC_ID']
    evttime  = evt_data["TIME"]
    gtistart = gti_data["START"]
    gtistop  = gti_data["STOP"]

    timeOffset = gtistart[0]
    print("LEN of GTI:", gtistop - gtistart)
    print("LEN of NON- GTI:", gtistart[1:] - gtistop[:-1])

    energy =  evt_data["ENERGY"][srcid==101]
    phi =  evt_data["DETPHI"][srcid==101]
    x =  evt_data["DETX"][srcid==101]
    y =  evt_data["DETY"][srcid==101]
    livetime  = evt_data["LIVETIME"]
    deltat = numpy.ediff1d(evttime[srcid==101])
    
    if False:
        plt.figure("SrcTime", figsize=(14, 3*len(DU_LIST)), tight_layout=True)
        plt.subplot(len(DU_LIST),1,DU)
        plt.plot(evttime-timeOffset, srcid, ".")
        
        plt.xlabel("Rel. Time")
        plt.ylabel("Scr Id")
        #plt.grid()
        plt.axis([-500, max(evttime)-timeOffset+500, -10, 120])
    
        for xc in gtistart:
            plt.axvline(x=xc-timeOffset, color='g')
        for xc in gtistop:
            plt.axvline(x=xc-timeOffset, color='r', linestyle='--')

    if True:
        R_c = 0
        N_c = 0
        CalCTI_start = gtistop[:-1]
        CalCTI_stop  = gtistart[1:]
        CalCTI_n = len(CalCTI_start)
        for T0, T1 in zip(CalCTI_start, CalCTI_stop):
            _timeCut = numpy.logical_and(evttime>T0, evttime<T1)
            LT = livetime[_timeCut]
            if len(LT)>0:
                print(T0, T1, len(LT), sum(LT), (1e6)*len(LT)/sum(LT))
                R_c += (1e6)*len(LT)/sum(LT)
                N_c += 1
        plt.figure("CalCDeltaT", figsize=(8,10))
        plt.subplot(len(DU_LIST),1,DU)
        plt.hist(deltat, range=(0, 0.1), bins=500, label = "mean %f - Rate %.1f" % (numpy.mean(deltat), R_c/N_c))
        plt.yscale("log")
        plt.legend()
    if False:
        plt.figure("CalCSpectrum", figsize=(8,10))
        plt.subplot(len(DU_LIST),1,DU)
        plt.hist(energy, range=(2, 10), bins=100, label = "DU %d" % DU)
    if False:
        plt.figure("CalCPhi", figsize=(8,10))
        plt.subplot(len(DU_LIST),1,DU)
        plt.hist(phi, range=(-3.2, 3.2), bins=100, label = "DU %d" % DU)
    if False:
        plt.figure("CalCDetXY", figsize=(8,10))
        plt.subplot(len(DU_LIST),1,DU)
        plt.hist2d(x,y, bins=(100,100), range= ((-8, 8), (-8,8)))

    
plt.show()
