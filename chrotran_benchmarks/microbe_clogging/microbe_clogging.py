# Id: microbe_clogging.py, Mon 28 Aug 2017 05:03:53 PM MDT #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for biomass clogging and unclogging.
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from numpy import ma
import h5py

filename = 'microbe_clogging.h5'
f = h5py.File(filename, 'r')

timesteps = range(1,200)
dt = 1.0
myvars = ['biomass [mol_m^3]','Permeability_X [m^2]','Porosity']
majorFormatter = plt.matplotlib.ticker.FormatStrFormatter("%0.1e")

keys = f.keys()
timekeys = [key for key in keys if 'Time:' in key]
times = [float(timekey.split('Time:  ')[1].split(' d')[0]) for timekey in timekeys]
st = sorted(times)
keydict = dict(zip(times,timekeys))

# get coordinates
xgrid = f['Coordinates']['X [m]'].value
ygrid = f['Coordinates']['Y [m]'].value
# x, y = np.meshgrid(xgrid, ygrid)
x, y = np.mgrid[min(xgrid):max(xgrid):xgrid[1]-xgrid[0], min(ygrid):max(ygrid):ygrid[1]-ygrid[0]]

results = dict()
for myvar in myvars:
    results[myvar] = []
    for timestep in timesteps:
        # get values at the well (center node)
        results[myvar].append(f[keydict[timestep*dt]][myvar].value[1][1][0])

f,ax = plt.subplots(1,3,figsize=(10,4))
for i in range(0,3):
    ax[i].plot([timestep*dt for timestep in timesteps],results[myvars[i]])
    ax[i].set_ylabel(myvars[i])
    ax[i].set_xlabel('Time [day]')

ax[0].yaxis.set_major_formatter(majorFormatter)
ax[1].yaxis.set_major_formatter(majorFormatter)

plt.suptitle("microbe_clogging benchmark")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
f.savefig('microbe_clogging.png')
