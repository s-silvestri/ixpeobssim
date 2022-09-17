from astropy.io import fits
import matplotlib.pyplot as plt


file_path = '/work/sn1006/acisf04391_repro_evt2.fits'

hdul = fits.open(file_path)
events = hdul[1]
print (events.data['energy'])
plt.hist (events.data['energy']/1e3, bins=250)
plt.xlim((2,8))
plt.xscale('log')
plt.yscale('log')
plt.show()