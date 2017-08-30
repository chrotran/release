# Id: abiotic.py, Mon 28 Aug 2017 05:03:53 PM MDT pandeys #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for abiotic reaction without retardation.
#------------------------------------------------------------------------------
import sys
import os
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import odespy

# ------------------------------------------------------------------------------
# Python numerical simulation
# ------------------------------------------------------------------------------
execfile('../python/pyfun.py')

# Parameters
pars = {'s'				: 1.0, # saturation
		'por'			: 0.25, # porosity
		'v_cell'		: 1.0, # m^3

		'alpha'			: 1.0, # [-]
		'B_min'			: 1.e-10, # [mol/m^3]
		'rho_b'			: 1.e20, # [g/L = M]

		'gamma_B'		: 0.0, # [L/mol/s]
		'gamma_CD'		: 1.0, # [L/mol/s]
		'gamma_X'		: 0.0, # [L/mol/s]

		'lambda_B1'		: 0.0, # [/s]
		'lambda_B2'		: 1.e-5, # [/s]
		'lambda_C'		: 0.0, # [/s]
		'lambda_D'		: 0.0, # [/s]
		'lambda_D_i'	: 150.e-2, # [/s]
		'lambda_D_m'	: 1.e-2, # [/s]

		'K_B'			: 5.e1, #  [mol/L_bulk]
		'K_C'			: 1.e-7, # [M]
		'K_D'			: 1.e-6, # [M]
		'K_I'			: 1.e-4, # [M]

		'S_C'			: 0.33, # [-]
		'S_D_1'			: 1.0, # [-]
		'S_D_2'			: 0.020833# [-]
}

# Solver options
sopt = {'T' :	0.1*3600, # end of simulation [s from hr]
		'N' :	1000, # no of timesteps
}

# Initial conditions
init = {'C'		: 1.e-2, # [M]
		'D_m'	: 1.e-2, # [M]
		'I'		: 1.e-20, # [M]
		'X'		: 1.e-20, # [M]
		'B'		: 1.e-20, # mol/m^3_bulk
		'D_i'	: 1.e-20, # mol/m^3_bulk
		'chubbite_vf' : 0.85, # [m^3/m^3_bulk]
}

u, t = run_ode()

L_water = pars['v_cell'] * pars['por'] * pars['s'] * 1.e3 # [L]
results_ode = {}
results_ode['time'] = t/3600 # [hr from s]
results_ode['C'] = u[:,0]/L_water
results_ode['D_m'] = u[:,1]/L_water
results_ode['I'] = u[:,2]/L_water
results_ode['X'] = u[:,3]/L_water
results_ode['B'] = u[:,4]/pars['v_cell']
results_ode['D_i'] = u[:,5]/pars['v_cell']
results_ode['chubbite_vf'] = u[:,6]

# ------------------------------------------------------------------------------
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "abiotic_mimt"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total molasses [M]', 'Total Cr(VI)', 'molasses_im [mol/m^3]']
observation_list = ['obs1']
results_pflotran =  getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
majorFormatter = plt.matplotlib.ticker.FormatStrFormatter("%0.1e")
mycmap=plt.cm.jet(np.linspace(0,1,5))

# First plot
fig = plt.figure(figsize=[10,5])
ax = fig.add_subplot(1, 2, 1)
xlims = [0,0.1]
skipfactor = 50 # skip data in ode results
fontsize = 9

pflo_plotvars = [variable_list[0:2], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['D_m','C']
legend_list = ['D_m - PFLOTRAN','Cr(VI) - PFLOTRAN', 'D_m - odespy','Cr(VI) - odespy']

plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

# Second plot
ax2 = fig.add_subplot(1, 2, 2)
pflo_plotvars = [[variable_list[2]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['D_i']
legend_list = ['D_i - PFLOTRAN','D_i - odespy']

plot_benchmarks(ax2, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [mol/m^3_bulk]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

plt.tight_layout()
plt.savefig(simbasename + '.png')
