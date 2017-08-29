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
	B = u[4] # [mol/m^3_bulk]
	D_i = u[5] # [mol/m^3_bulk]

	D = D_m + D_i * immobile_to_water_vol # [mol/L water]
	mu_B = pars['lambda_B1'] * B * D/(D + pars['K_D']) * (pars['K_B'] / (pars['K_B'] + B))**pars['alpha'] * pars['K_I']/(pars['K_I'] + I) # [mol/m3 bulk/s]
	mu_CD = pars['gamma_CD'] * C * D # mol/L/s

	# calculate changes to residuals
	R_C = -pars['lambda_C']*B*C/(pars['K_C']+C)*pars['v_cell'] - pars['S_C']*mu_CD*pars['v_cell']*immobile_to_water_vol

	R_D_m = (-pars['S_D_1'] * (D_m/D) * mu_B * pars['v_cell'] -
		pars['lambda_D'] * (D_m/D) * B * pars['v_cell'] -
		pars['S_D_2'] * (D_m/D) * mu_CD * pars['v_cell'] * immobile_to_water_vol -
		pars['lambda_D_i'] * D_m * pars['v_cell'] * immobile_to_water_vol +
		pars['lambda_D_m'] * D_i * pars['v_cell'])

	R_D_i = (-pars['S_D_1'] * (D_i/D) * mu_B * pars['v_cell'] -
		pars['lambda_D'] * (D_i/D) * B * pars['v_cell'] -
		pars['S_D_2'] * (D_i/D) * mu_CD * pars['v_cell'] * immobile_to_water_vol +
		pars['lambda_D_i'] * D_m * pars['v_cell'] * immobile_to_water_vol -
		pars['lambda_D_m'] * D_i * pars['v_cell'])
		
	return [R_C,R_D_m,0,0,0,R_D_i,0]

def run_ode():
	# solver = odespy.RK4(f)
	solver = odespy.CashKarp(chrotran_sandbox)
	solver.set_initial_condition([
		init['C']*pars['por']*pars['s']*1e3, # C
		init['D_m']*pars['por']*pars['s']*1e3, # D_m
		init['I']*pars['por']*pars['s']*1e3, # I
		init['X']*pars['por']*pars['s']*1e3, # X
		init['B'], # B
		init['D_i'], # D_i
		init['chubbite_vf'], # D_i
		])

	T = 0.1*3600 # end of simulation [s]
	N = 1000 # no of timesteps
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
	for item in combined_var_obs_list:
		for key in combined_dict.keys():
			if item[0] in key and item[1] in key:
				var_new = [v for v in combined_dict[key]]
				combined_dict_trimmed[item[0] + " " + item[1]] = var_new

	return time, combined_dict_trimmed

def plot_benchmarks():
	print("To do")
