# Id: microbe_enzymatic.py, Mon 28 Aug 2017 05:03:53 PM MDT #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for chromium reduction by biomass.
#------------------------------------------------------------------------------
import sys
import os
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import odespy
sys.path.append('../python')
import pyfun as pf

# ------------------------------------------------------------------------------
# Python numerical simulation
# ------------------------------------------------------------------------------
# Parameters
pars = {'s'				: 1.0, # saturation
		'por'			: 0.25, # porosity
		'v_cell'		: 1.0, # m^3_bulk

		'alpha'			: 1.0, # [-]
		'B_min'			: 1.e-10, # [mol/m^3_bulk]
		'rho_b'			: 1.e20, # [g/L = M]

		'gamma_B'		: 0.0, # [L/mol/s]
		'gamma_CD'		: 0.0, # [L/mol/s]
		'gamma_X'		: 0.0, # [L/mol/s]

		'lambda_B1'		: 0.0, # [/s]
		'lambda_B2'		: 0.0, # [/s]
		'lambda_C'		: 1.e-4, # [/s]
		'lambda_D'		: 0.0, # [/s]
		'lambda_D_i'	: 0.0, # [/s]
		'lambda_D_m'	: 0.0, # [/s]

		'K_B'			: 5.e1, #  [mol/m^3_bulk]
		'K_C'			: 1.e-5, # [M]
		'K_D'			: 1.e-6, # [M]
		'K_I'			: 1.e-4, # [M]

		'S_C'			: 0.33, # [-]
		'S_D_1'			: 1.0, # [-]
		'S_D_2'			: 0.020833# [-]
}

# Solver options
sopt = {'T' :	10*3600, # end of simulation [s from hr]
		'N' :	1000, # no of timesteps
}

# Initial conditions
init = {'C'		: 1.e-2, # [M]
		'D_m'	: 1.e-20, # [M]
		'I'		: 1.e-20, # [M]
		'X'		: 1.e-20, # [M]
		'B'		: 1.e0, # mol/m^3_bulk
		'D_i'	: 1.e-20, # mol/m^3_bulk
		'chubbite_vf' : 1-pars['por'], # [m^3/m^3_bulk]  #?????
}

chrotran_sandbox = pf.make_chrotran_sandbox(pars)
u, t = pf.run_ode_chrotran(init, pars, sopt, chrotran_sandbox)

# L_water = pars['v_cell'] * pars['por'] * pars['s'] * 1.e3 # [L]
results_ode = {}
results_ode['time'] = t/3600 # [hr from s]
results_ode['C'] = u[:,0]
results_ode['D_m'] = u[:,1]
results_ode['I'] = u[:,2]
results_ode['X'] = u[:,3]
results_ode['B'] = u[:,4]
results_ode['D_i'] = u[:,5]
results_ode['chubbite_vf'] = u[:,6]

# ------------------------------------------------------------------------------
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "microbe_enzymatic"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total Cr(VI)']
observation_list = ['obs1']
results_pflotran =  pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting/Regression
# ------------------------------------------------------------------------------
# First plot
fig = plt.figure(figsize=[5,5])
ax = fig.add_subplot(1, 1, 1)
# xlims = [0,5]
skipfactor = 25 # skip data in ode results
fontsize = 9

pflo_plotvars = [[variable_list[0]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['C']
legend_list = ['C - PFLOTRAN', 'C - odespy']

pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]", ylabel="Concentration [M]", skipfactor=skipfactor, fontsize=fontsize)

regression_result = pf.calc_regression(ts = 1.0e-1,tol = 1.0e-7,results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)

plt.suptitle("microbe_enzymatic benchmark")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(simbasename + '.png')
