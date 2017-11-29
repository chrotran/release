import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from numpy import ma
import h5py

filename = 'CT2.h5'
f = h5py.File(filename, 'r')

timestep = 50
dt = 4.0
myvars = ['biomass [mol_m^3]','Liquid X-Velocity [m_per_d]','Liquid Y-Velocity [m_per_d]']
qsf = 2 # quiver plot skip factor
majorFormatter = plt.matplotlib.ticker.FormatStrFormatter("%0.0e")

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
    results[myvar] = f[keydict[timestep*dt]][myvar].value
    results[myvar] = results[myvar].reshape(len(xgrid)-1,len(ygrid)-1)

f,ax = plt.subplots()

B = plt.pcolor(x, y, results['biomass [mol_m^3]'], cmap='Greens', vmin = 0, vmax = 2.0e5)
cbar = plt.colorbar(format=majorFormatter)
cbar.set_label('Biomass [g/m^3]', rotation=270,labelpad=20)
Q = plt.quiver(x[::qsf,::qsf], y[::qsf,::qsf], results['Liquid X-Velocity [m_per_d]'][::qsf,::qsf], results['Liquid Y-Velocity [m_per_d]'][::qsf,::qsf],units='x',scale=0.1,angles='xy')

ax.set_xlim(min(xgrid),max(xgrid))
ax.set_ylim(min(ygrid),max(ygrid))
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
 
f.savefig('example_2_'+str(int(timestep*dt))+'d.png',dpi=600)
plt.close()
