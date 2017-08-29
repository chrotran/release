# Id: abiotic.py, Mon 28 Aug 2017 05:03:53 PM MDT pandeys #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for abiotic reaction without retardation.
#------------------------------------------------------------------------------
import sys
sys.path.append('../python/')
import pyfun as pf
import os
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

# ------------------------------------------------------------------------------
# Python numerical simulation
# ------------------------------------------------------------------------------
def chrotran_sandbox(u, t):
	C = u[0]/L_water # [M]
	D_m = u[1]/L_water # [M]
	I = u[2]/L_water # [M]
	X = u[3]/L_water # [M]
	B = u[4] # [mol/m^3_bulk]
	D_i = u[5] # [mol/m^3_bulk]

	D = D_m + D_i * immobile_to_water_vol # [mol/L water]
	mu_B = lambda_B1 * B * D/(D + K_D) * (K_B / (K_B + B))**alpha * K_I/(K_I + I) # [mol/m3 bulk/s]
	mu_CD = gamma_CD * C * D # mol/L/s

	# calculate changes to residuals
	R_C = -lambda_C*B*C/(K_C+C)*volume - S_C*mu_CD*volume*immobile_to_water_vol

	R_D_m = (-S_D_1 * (D_m/D) * mu_B * volume -
		lambda_D * (D_m/D) * B * volume - 
		S_D_2 * (D_m/D) * mu_CD * volume * immobile_to_water_vol - 
		lambda_D_i * D_m * volume * immobile_to_water_vol + 
		lambda_D_m * D_i * volume)

	return [R_C,R_D_m,0,0,0,0]

# Parameters
s = 1.0
por = 0.25 
volume = 1.0 # m^3
L_water = volume * por * s * 1.e3
immobile_to_water_vol = por * s * 1.e3 # L water / m3_bulk

alpha = 1.0
B_min = 1.e-10 # [mol/m^3]
gamma_B = 0.0 # [L/mol/s]
gamma_CD = 1.0 # [L/mol/s]
gamma_X = 0.0 # [L/mol/s]
lambda_B1 = 0.0 # [/s]
lambda_B2 = 1.e-5 # [/s]
lambda_C = 0.0 # [/s]
lambda_D = 1.e-5 # [/s]
lambda_D_m = 0.0 # [/s]
lambda_D_i = 0.0 # [/s]

K_B = 5.e1 #  [mol/L_bulk]
K_C = 1.e-7 # [M]
K_D = 1.e-6 # [M]
K_I = 1.e-4 # [M]
       
rho_b = 1.e20 # [g/L = M]

S_C = 0.33 # [-]
S_D_1 = 1.0 # [-]
S_D_2 = 0.020833# [-]

import odespy
# solver = odespy.RK4(f)
solver = odespy.CashKarp(chrotran_sandbox)
solver.set_initial_condition([
	1.e-2*por*s*1e3, # C
	1.e-2*por*s*1e3, # D_m
	1.e-20*por*s*1e3, # I
	1.e-20*por*s*1e3, # X
	1.e-20, # B
	1.e-20, # D_i
	])  

T = 0.1*3600 # end of simulation [s]
N = 1000 # no of timesteps
time_points = np.linspace(0,T,N+1)
u, t = solver.solve(time_points)

C = u[:,0]/L_water
D_m = u[:,1]/L_water
I = u[:,2]/L_water
X = u[:,3]/L_water
B = u[:,4]
D_i = u[:,5]

# ------------------------------------------------------------------------------
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "abiotic"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total molasses [M]', "Total Cr(VI)"]
observation_list = ['obs1']
plot_filename = 'test.png'
x_label = 'time [h]'
y_label = 'Concentration [M]'
xrange=([0,0.1])

majorFormatter = plt.matplotlib.ticker.FormatStrFormatter("%0.1e")
combined_var_obs_list = [variable_list, observation_list]
combined_var_obs_list = list(it.product(*combined_var_obs_list))

time, data = pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
legend_list = ['molasses - PFLOTRAN','Cr(VI) - PFLOTRAN', 'molasses - odespy','Cr(VI) - odespy']
xrange = [0,0.1]
fontsize = 9
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
lns = []

for item in combined_var_obs_list:
	ln, = ax.plot(time, data[item[0] + " " + item[1]], linestyle='-')
	lns.append(ln)
ax.yaxis.set_major_formatter(majorFormatter)

ln, =  ax.plot(t[::50]/3600,C[::50],c="g",ls=' ',marker = 'o')
lns.append(ln)
ln, =  ax.plot(t[::50]/3600,D_m[::50],c="b",ls=' ',marker = 'o')
lns.append(ln)

ax.set_xlim(xrange)
ax.set_xlabel("Time [hr]")
ax.legend(lns, legend_list, ncol=1, fancybox=True, shadow=False, prop={'size': str(fontsize)}, loc='best')

plt.savefig(simbasename + '.png')
