import numpy as np
import itertools as it
import odespy
import matplotlib.pyplot as plt

liter_b_to_m3_b = 1e3

def make_chrotran_sandbox(pars):
    def chrotran_sandbox(u, t):
        C = u[0] # [M]
        D_m = u[1] # [M]
        I = u[2] # [M]
        X = u[3] # [M]
        B = u[4] # [mol/m3_bulk]
        D_i = u[5] # [mol/m3_bulk]
        chubbite_vf = u[6]

        theta0 = pars['por'] * pars['s'] * liter_b_to_m3_b # L_water / m3_bulk
        #theta0 = max(0.01,(1 -chubbite_vf)) * pars['s'] * liter_b_to_m3_b # L water / m3_bulk

        D = D_m + D_i / theta0 # [mol/L water]

        mmf = D_m/D # mobile mole fraction
        immf = 1-mmf

        mu_B = pars['lambda_B1'] * B * D/(D + pars['K_D']) * (pars['K_B'] / (pars['K_B'] + B))**pars['alpha'] * pars['K_I']/(pars['K_I'] + I) # [mol/m^3_bulk/s]

        mu_CD = pars['gamma_CD'] * C * D # [mol/L_water/s]

        # MOBILE DERIVATIVES [mol/L_water/s]
        dC_dt = (- pars['lambda_C']*(B/theta0)*C/(pars['K_C']+C)
                 - pars['S_C']*mu_CD
                 )

        dDm_dt = (- pars['S_D_1'] * mmf * (mu_B/theta0)
                  - pars['lambda_D'] * mmf * (B/theta0)
                  - pars['S_D_2'] * mmf * mu_CD
                  - pars['lambda_D_i'] * D_m
                  + pars['lambda_D_m'] * (D_i/theta0)
                  )
        
        dI_dt = 0.0
        
        dX_dt = - pars['gamma_X'] * X * (B/theta0)
        
        # IMMOBILE DERIVATIVES [mol/m^3_bulk/s]
        dB_dt = (  mu_B
                 - pars['lambda_B2'] * (B - pars['B_min'])
                 - pars['gamma_B'] * (B - pars['B_min']) * X
                 )

        dDi_dt = (- pars['S_D_1'] * immf * mu_B
                  - pars['lambda_D'] * immf * B
                  - pars['S_D_2'] * immf * mu_CD * theta0
                  + pars['lambda_D_i'] * D_m * theta0
                  - pars['lambda_D_m'] * D_i
                  )

        return [dC_dt, dDm_dt, dI_dt, dX_dt, dB_dt, dDi_dt, 0]
    
    return chrotran_sandbox

def run_ode(init, pars, sopt, function):
	# solver = odespy.RK4(function)
	solver = odespy.CashKarp(function)
	solver.set_initial_condition([
		init['C'], # [M]
		init['D_m'], # [M]
		init['I'], # [M]
		init['X'], # [M]
		init['B'], # [mol/m3_bulk]
		init['D_i'], # [mol/m3_bulk]
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

def plot_benchmarks(ax,results_ode = {}, results_pflotran = {}, ode_plotvars =[], pflo_plotvars = [], legend_list=[], xlabel='', ylabel='', xlims=[], ylims=[], skipfactor=1, fontsize=10, mycmap=plt.cm.jet(np.linspace(0,1,5)), majorFormatter=plt.matplotlib.ticker.FormatStrFormatter("%0.1e")):
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

	return lns
