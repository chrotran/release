# Id: microbe_clogging.py, Mon 28 Aug 2017 05:03:53 PM MDT #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for biomass clogging and unclogging.
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from numpy import ma
import sys
sys.path.append('../python')
import pyfun as pf
import itertools as it

timesteps = range(0,50)
dt = 4.0
myvars = ['biomass [mol_m^3]','Permeability_X [m^2]','Porosity']
simbasename = "microbe_clogging"

# New simulation
filename = simbasename + '.h5'
results_new = pf.readh5_1d(filename=filename,timesteps=timesteps,dt=dt,myvars=myvars)

# Old simulation
filename = simbasename + '_gold.h5'
results_gold = pf.readh5_1d(filename=filename,timesteps=timesteps,dt=dt,myvars=myvars)

f,ax = plt.subplots(1,3,figsize=(10,4))
mycmap=plt.cm.jet(np.linspace(0,1,5))
majorFormatter = plt.matplotlib.ticker.FormatStrFormatter("%0.1e")
skipfactor = 3
for i in range(0,3):
	ax[i].plot([timestep*dt for timestep in timesteps],results_new[myvars[i]],c=mycmap[0])
	ax[i].scatter([timestep*dt for timestep in timesteps][0:-1:skipfactor],results_gold[myvars[i]][0:-1:skipfactor],c=mycmap[0])
	ax[i].set_ylabel(myvars[i])
	ax[i].set_xlabel('Time [day]')

ax[0].yaxis.set_major_formatter(majorFormatter)
ax[1].yaxis.set_major_formatter(majorFormatter)

plt.suptitle("microbe_clogging benchmark")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
f.savefig('microbe_clogging.png')

# ------------------------------------------------------------------------------
# Calculate regression test error
# ------------------------------------------------------------------------------
# For regression test, use gold standard file instead of python ODE

results_new['time']=[timestep*dt for timestep in timesteps]
results_gold['time']=[timestep*dt for timestep in timesteps]
results_new['biomass [mol_m^3] '] = results_new['biomass [mol_m^3]']
observation_list = ['']
variable_list = ['biomass [mol_m^3]']
pflo_plotvars = [[variable_list[0]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars  = ['biomass [mol_m^3]']
regression_result = pf.calc_regression(ts = 20.0,tol = 1.0e-9,results_ode=results_gold, results_pflotran=results_new, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)
