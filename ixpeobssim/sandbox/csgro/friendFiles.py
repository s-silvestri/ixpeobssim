import numpy
from ixpeobssim.evt.event import xEventFile, xEventFileFriend
#from ixpeobssim.irfgen import NUM_CHANNELS, channel_to_energy

## code moved to ixpeobssim.evt.event.py

if __name__ == '__main__':
# For CenX3
    commonpath = "/run/media/carmelo/TOSHIBAEXT/ixpedata/Science/CenX3/"
    l1path = []
    l2path = []
    for DuId in [1,2]:
        l1path.append(commonpath+"event_l1/ixpe01006501_det%d_evt1_v02.fits" % DuId)
        l2path.append(commonpath+"event_l2/ixpe01006501_det%d_evt2_v02.fits" % DuId)

    
    myFiles = xEventFileFriend(l2path, l1path)
    #myFiles = xEventFileFriend('/home/carmelo/xpe/xpedata/Science/smcx1/ixpe01903701_det1_evt2_v09.fits', '/home/carmelo/xpe/xpedata/Science/smcx1/ixpe01903701_det1_evt1_v06.fits')

    t2 = myFiles.l2value("TIME")
    t1 = myFiles.l1value("TIME")
    t1a = myFiles.l1value("TIME", True)

    print(t2)
    print(t1)
    print(t1a)
    

    tgr_id_l2 = myFiles.l2value("TRG_ID")
    tgr_id_l1 = myFiles.l1value("TRG_ID")
    print("Trigger Id diff check, must be zero, and the winner is:",
          sum(tgr_id_l1 - tgr_id_l2))


    energy   = myFiles.energy_data()
    num_clu  = myFiles.l1value("NUM_CLU")
    trk_size = myFiles.l1value("TRK_SIZE")
    trk_bord = myFiles.l1value("TRK_BORD")
    
    import matplotlib.pyplot as plt


    plt.figure("Energy")
    plt.hist(energy, range=(0,20), bins = 100)

    plt.figure("Num Clu")
    plt.hist(num_clu, range=(0,21), bins = 20)

    plt.figure("TKR Board")
    plt.hist(trk_bord, range=(0,51), bins = 50)

    plt.figure()
    plt.hist2d(energy, trk_size)

    plt.show()
    
