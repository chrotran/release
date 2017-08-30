import numpy as np
import itertools as it
import odespy

def chrotran_sandbox(u, t):
	L_water = pars['v_cell'] * pars['por'] * pars['s'] *1.e3 # [L]
	immobile_to_water_vol = pars['por'] * pars['s'] * 1.e3 # L water / m3_bulk

	C = u[0]/L_water # [M]
	D_m = u[1]/L_water # [M]
	I = u[2]/L_water # [M]
	X = u[3]/L_water # [M]
	B = u[4]/pars['v_cell'] # [mol/m^3_bulk]
	D_i = u[5]/pars['v_cell'] # [mol/m^3_bulk]
	chubbite_vf = u[6]

	D = D_m + D_i / immobile_to_water_vol # [mol/L water]
	mobile_mole_fraction = D_m/D
	immobile_mole_fraction = 1-mobile_mole_fraction

	mu_B = pars['lambda_B1'] * B * D/(D + pars['K_D']) * pars['K_B'] / (pars['K_B'] + B)**pars['alpha'] * pars['K_I']/(pars['K_I'] + I) # [mol/m3 bulk/s]
	mu_CD = pars['gamma_CD'] * C * D # mol/L/s

	# calculate changes to residuals
	R_C = -pars['lambda_C']*B*C/(pars['K_C']+C)*pars['v_cell'] - pars['S_C']*mu_CD*pars['v_cell']*immobile_to_water_vol

	R_D_m = (-pars['S_D_1'] * (mobile_mole_fraction) * mu_B * pars['v_cell'] -
		pars['lambda_D'] * (mobile_mole_fraction) * B * pars['v_cell'] -
		pars['S_D_2'] * (mobile_mole_fraction) * mu_CD * pars['v_cell'] * immobile_to_water_vol -
		pars['lambda_D_i'] * D_m * pars['v_cell'] * immobile_to_water_vol +
		pars['lambda_D_m'] * D_i * pars['v_cell'])

	R_D_i = (-pars['S_D_1'] * (immobile_mole_fraction) * mu_B * pars['v_cell'] -
		pars['lambda_D'] * (immobile_mole_fraction) * B * pars['v_cell'] -
		pars['S_D_2'] * (immobile_mole_fraction) * mu_CD * pars['v_cell'] * immobile_to_water_vol +
		pars['lambda_D_i'] * D_m * pars['v_cell'] * immobile_to_water_vol -
		pars['lambda_D_m'] * D_i * pars['v_cell'])

	return [R_C,R_D_m,0,0,0,R_D_i,0]

def run_ode():
	# solver = odespy.RK4(chrotran_sandbox)
	solver = odespy.CashKarp(chrotran_sandbox)
	solver.set_initial_condition([
		init['C']*pars['por']*pars['s']*1e3, # [moles]
		init['D_m']*pars['por']*pars['s']*1e3, # [moles]
		init['I']*pars['por']*pars['s']*1e3, # [moles]
		init['X']*pars['por']*pars['s']*1e3, # [moles]
		init['B']*pars['v_cell'], # [moles]
		init['D_i']*pars['v_cell'], # [moles]
		init['chubbite_vf'],
		])

	time_points = np.linspace(0,sopt['T'],sopt['N']+1)
	u, t = solver.solve(time_points)

	return u, t

def getobsdata(variable_list=[], observation_list=[], observation_filenames=[]):
	"""
	Get observation data from pflotran tec file
	"""
	combined_dict = {}
	for file in observation_filenames:
		variable = []
		f = open(file, 'r')
		title = f.readline()
		title = title.split(',')
		for i in title:
			variable.append(i.strip('"'))
			data = np.genfromtxt(file, skip_header=1)
			data = data.T.tolist()
			var_values_dict = dict(zip(variable, data))
			combined_dict.update(var_values_dict)

		for key in combined_dict.keys():
			if 'Time' in key:
				time = combined_dict[key]

		combined_var_obs_list = [variable_list, observation_list]
		combined_var_obs_list = list(it.product(*combined_var_obs_list))

	combined_dict_trimmed = {}
	combined_dict_trimmed['time']= time
	for item in combined_var_obs_list:
		for key in combined_dict.keys():
			if item[0] in key and item[1] in key:
				var_new = [v for v in combined_dict[key]]
				combined_dict_trimmed[item[0] + " " + item[1]] = var_new

	return combined_dict_trimmed

def plot_benchmarks(ax,results_ode = {}, results_pflotran = {}, ode_plotvars =[], pflo_plotvars = [], legend_list=[], xlabel='', ylabel='', xlims=[], ylims=[], skipfactor=1, fontsize=10):
	"""
	Plot data to an axis object
	"""
	lns = []
	ctr = 0
	for item in pflo_plotvars:
		ln, = ax.plot(results_pflotran['time'], results_pflotran[item[0] + " " + item[1]], linestyle='-',c=mycmap[ctr])
		lns.append(ln)
		ctr =+ 1

	ctr = 0
	for item in ode_plotvars:
		ln, =  ax.plot(results_ode['time'][::skipfactor],results_ode[item][::skipfactor],ls=' ',marker = 'o',c=mycmap[ctr])
		lns.append(ln)
		ctr =+ 1

	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if xlims: ax.set_xlim(xlims)
	if ylims: ax.set_ylim(ylims)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.legend(lns, legend_list, ncol=1, fancybox=True, shadow=False, prop={'size': str(fontsize)}, loc='best')
