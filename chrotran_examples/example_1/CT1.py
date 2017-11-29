import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from numpy import ma
import h5py

filenames = ['CT1di.h5','CT1xx.h5','CT1xi.h5','CT1dx.h5']


timestep = 100
dt = 5.0
myvars = ['Total_Cr(VI) [M]']
majorFormatter = plt.matplotlib.ticker.FormatStrFormatter("%0.0e")
MW_Cr = 52.0

for i in np.arange(0,4):
	f = h5py.File(filenames[i], 'r')
	keys = f.keys()
	timekeys = [key for key in keys if 'Time:' in key]
	times = [float(timekey.split('Time:  ')[1].split(' d')[0]) for timekey in timekeys]
	st = sorted(times)
	keydict = dict(zip(times,timekeys))

	# get coordinates
	xgrid = f['Coordinates']['X [m]'].value
	ygrid = f['Coordinates']['Y [m]'].value
	x, y = np.meshgrid(xgrid, ygrid)
	# x, y = np.mgrid[min(xgrid):max(xgrid):xgrid[1]-xgrid[0], min(ygrid):max(ygrid):ygrid[1]-ygrid[0]]

	results = dict()
	for myvar in myvars:
		results[myvar] = f[keydict[timestep*dt]][myvar].value
		results[myvar] = results[myvar].reshape(len(xgrid)-1,len(ygrid)-1)

	plt.subplot(2,2,i+1)
	plt.pcolor(x, y, np.transpose(results['Total_Cr(VI) [M]'])*MW_Cr*10**6, vmin = 0.0, vmax = 1000)
	cbar = plt.colorbar(format=majorFormatter)
	cbar.set_label('Cr(VI) [ppb]', rotation=270,labelpad=20)
	
	plt.xlim(min(xgrid),max(xgrid))
	plt.ylim(min(ygrid),max(ygrid))
	plt.title(filenames[i])

	if i == 0 or i == 3:
		plt.xlabel('x [m]')
	if i == 1 or i == 3:
		plt.ylabel('y [m]')

plt.tight_layout()
plt.savefig('example_1_'+str(int(timestep*dt))+'d.png',dpi=600)
plt.close()
