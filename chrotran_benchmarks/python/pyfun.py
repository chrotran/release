import numpy as np
import itertools as it

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
		
	R_D_i = (-S_D_1 * (D_i/D) * mu_B * volume -
		lambda_D * (D_i/D) * B * volume - 
		S_D_2 * (D_i/D) * mu_CD * volume * immobile_to_water_vol + 
		lambda_D_i * D_m * volume * immobile_to_water_vol - 
		lambda_D_m * D_i * volume)

	return [R_C,R_D_m,0,0,0,R_D_i]

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
