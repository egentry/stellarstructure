import numpy as np
import astropy.constants as const
from scipy import optimize
import gasproperties
import opacity
import energygeneration
import modelparameters



def load1(m, P_c, T_c, mu, X, Y, Z, verbose=False):
	"""
		Determines initial conditions of *center*
			(Technically at center + epsilon of mass)

		Assumes:
			mass, m as independent variable

		Inputs:
			m 			1x1 float 	- initial mass step [g]
			P_c			1x1 float	- central Pressure [dyne cm^-2]
			T_c			1x1 float	- central Temperature [K]
			mu			1x1 float	- mean molecular weight
			X			1x1 float	- Hydrogen mass fraction
			Y			1x1 float	- Helium mass fraction
			Z			1x1 float	- Metals mass fraction


		Outputs:
			r 			1x1 float	- radius [cm]
			l 			1x1 float	- luminosity [ergs s^-1]
			P 			1x1 float	- Pressure [dyne cm^-2]
			T 			1x1 float	- Temperature [K]


		Warnings:

	"""

	beta 	= gasproperties.calculate_beta(P_c, T_c)
	rho 	= gasproperties.calculate_density(P_c * beta, T_c, mu)

	#approx r(m)
	r = np.power( np.divide(3. * m, 4 * np.pi * rho), 1./3) 


	#approx l(m)
	epsilon 	= energygeneration.calculate_epsilon(T_c, rho, X, Y, Z)
	l = epsilon * m


	#calculate P(m)
	P = P_c - np.divide(3 * const.G.cgs.value, 8 * np.pi) * \
	np.power(np.divide(4 * np.pi * rho, 3),4./3) * np.power(m, 2./3)


	#approx T(m)
	kappa 	= opacity.calculate_opacity(T_c,rho)
	del_rad = gasproperties.calculate_del_rad(m, l, P_c, T_c, kappa)
	del_ad 	= gasproperties.calculate_del_ad()

	if del_rad < del_ad:	#Eq 11.9, Kippenhan Weigert Weiss 2nd ed. 
		T = np.power(np.power(T_c, 4.) - \
			np.divide(1., 8 * np.pi * const.sigma_sb.cgs.value) * \
			np.power(np.divide(3, 4 * np.pi), 2./3) * kappa * \
			epsilon * np.power(rho, 4./3) * np.power(m, 2./3), 1./4)

	else:
		T = T_c * np.exp(-1 * np.power(np.pi/6, 1./3) * \
			const.G.cgs.value * del_ad * \
			np.divide(np.power(rho, 4./3) * np.power(m, 2./3), P_c))	# convective


	if verbose is True:
		print ""
		print "load1 results)"
		# print "log10 rho: ", log10_rho, "[log10 g cm^-3]"
		print "rho: ", rho, "[g cm^-3]"
		# print "logR: ", opacity.calculate_log_R(np.log10(T), np.log10(rho))
		# print "log kappa: ", log_kappa, " log [cm^2 g^-1]"
		print "kappa: ", kappa, "[cm^2 g^-1]"
		print "pressure: ", P, "[dyne cm^-2]"
		print "Temp: ", T, "[K]"
		print "r: ", r, "[cm]"
		print ""

	return r, l, P, T


def load2(M, R, L, mu,
		log10_rho_bound_min=-10., log10_rho_bound_max=-6, 
		verbose=False):
	"""
		Determines initial conditions of *surface*

		Assumes:
			mass, m as independent variable
			opacity dominated by thomson scattering

		Inputs:
			mu			1x1 float	- mean molecular weight
			M 			1x1 float	- total mass [g]
			L 			1x1 float	- total luminosity [cm]
			R 			1x1 float	- total radius [cm]

			# opacity_interp_fn 	function 	- opacity table interpolater
			# 	- defined in opacity.py

		Outputs:
			r 			1x1 float	- radius [cm]
			l 			1x1 float	- luminosity [ergs s^-1]
			P 			1x1 float	- Pressure [dyne cm^-2]
			T 			1x1 float	- Temperature [K]


		Warnings:

	"""
	g = const.G.cgs.value * M / (R**2.)

	# kappa = .20 * (2 / mu_e) # assuming thomson scattering; Eq 17.2, 4.30; Kippenhan Weigert Weiss 2nd ed.

	r = R

	l = L

	T = (L / (4 * np.pi * R**2. * const.sigma_sb.cgs.value) )**.25


	args_in 	= (g, T, mu)
	if verbose is True:
		visualize_solver_bounds(g, T, mu, 
			log10_rho_bound_min, log10_rho_bound_max, show=True)
	log10_rho 	= optimize.brentq(log_rho_to_solve, 
					log10_rho_bound_min, log10_rho_bound_max,
					args=args_in)
	# log10_rho 	= optimize.newton(log_rho_to_solve, log10_rho_initial,
	# 				args=args_in)
	rho 		= np.power(10., log10_rho)
	kappa 		= opacity.calculate_opacity(T, rho)
	log_kappa 	= np.log10(kappa)


	P 			= (2./3) * g / kappa	

	if verbose is True:
		print "load2 results)"
		print "log10 rho: ", log10_rho, "[log10 g cm^-3]"
		print "rho: ", rho, "[g cm^-3]"
		print "logR: ", opacity.calculate_log_R(np.log10(T), log10_rho)
		print "log kappa: ", log_kappa, " log [cm^2 g^-1]"
		print "kappa: ", kappa, "[cm^2 g^-1]"
		print "pressure: ", P, "[dyne cm^-2]"
		print "Temp: ", T, "[K]"



	return r, l, P, T

def visualize_solver_bounds(g, T, mu, 
		log10_rho_bound_min, log10_rho_bound_max, show=False):
	"""
		For choosing bound ranges, and debugging
	"""
	import matplotlib.pyplot as plt
	from matplotlib import rc
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	rc('text', usetex=True)

	label_fontsize  = 20
	legend_fontsize = 18
	ticks_fontsize  = 16
	plt.rcParams.update({'font.size': ticks_fontsize})


	test_rhos = np.arange(-10, -2, .01)
	test_vals = log_rho_to_solve(test_rhos, g, T, mu)
	plt.plot(10**test_rhos, test_vals)
	plt.plot(10**test_rhos, test_vals*0)
	plt.xscale('log')
	plt.xlabel(r"$ \rho $ [g cm$^{-3}$]", fontsize=label_fontsize)
	plt.ylabel(r"$1- \frac{P_\mathrm{opacity}}{P_{\mathrm{EOS}}}$", fontsize=label_fontsize)
	plt.ylim(-1,1)
	# plt.title(r"$M_\ast = 1 M_\odot$", fontsize=label_fontsize)
	# plt.savefig("plots/opacity_1Msolar.eps", bbox_inches='tight')
	if show is True:
		plt.show()

	return

def log_rho_to_solve(log10_rho, g, T, mu):
	"""
		Calculates the differences between pressure using:
			Eq 1) P = P_gas + P_rad
					and
			Eq 2) P = (2/3) * g/kappa (at surface)

		Used for solving for the rho,
			in order to determine kappa,
			in order to determine pressure

		Assumes:
			You know the equation of state (Eq 1 above)

		Inputs:
			log10_rho	1x1 float	- log10 of density [g cm^-3]
			args 		1x1 array	
				L 			1x1 float	- total luminosity [cm]
				R 			1x1 float	- total radius [cm]


		Outputs:
			difference 	1x1 float	- fractional difference between Eq 1 and Eq 2 (above)
				solver is trying to minimize difference

		Warnings:

	"""
	difference = 1- np.divide(opacity_pressure(g, T, log10_rho),
			gas_and_rad_pressure(T, log10_rho, mu))
	return difference

def opacity_pressure(g, T, log10_rho):
	# print 'in opacity pressure'
	# print opacity.calculate_log_R(np.log10(T), log10_rho)
	return (2./3) * np.divide(g, opacity.calculate_opacity(T, np.power(10., log10_rho)))

def gas_and_rad_pressure(T, log10_rho, mu):
	return gasproperties.calculate_rad_pressure(T) + \
		np.power(10., log10_rho) * const.k_B.cgs.value * T / (mu * const.m_p.cgs.value)

def test_pressures(log10_rho, g, T, mu):
	print "for log10_rho = ", log10_rho
	print "opacity presure: ", opacity_pressure(g, T, log10_rho)
	print "gas + rad pressure: ", gas_and_rad_pressure(T, log10_rho, mu)
	return

def test_load1_and_load2():
	P_c = modelparameters.P_c # core pressure, [dyne cm^-2]
	T_c = modelparameters.T_c # core temperature, [K]

	X 	= modelparameters.X
	Y 	= modelparameters.Y
	Z 	= modelparameters.Z
	mu 	= modelparameters.mu

	m 	= 1e-6

	r, l, P, T 	= load1(m, P_c, T_c, mu, X, Y, Z, verbose=True)

	print ""
	print "load1 results:"
	print "r = ", r, "[cm]"
	print "l = ", l, "[erg s^-1]"
	print "P = ", P, "[dyne cm^-2]"
	print "T = ", T, "[K]"
	print ""

	M_star 			= modelparameters.M_star
	R_star 			= modelparameters.R_star
	L_star 			= modelparameters.L_star
	g = const.G.cgs.value * M_star / (R_star**2.)

	beta 	= gasproperties.calculate_beta(P, T)
	rho 	= gasproperties.calculate_density(P * beta, T, mu)
	kappa 	= opacity.calculate_opacity(T, rho)
	# print rho
	# print kappa


	r, l, P, T = load2(M_star, R_star, L_star, mu, verbose=True)
	# r, l, P, T = load2(M_star, R_star, L_star, mu, verbose=False)


	print ""
	print "load2 results:"
	print "r = ", r, "[cm]"
	print "l = ", l, "[erg s^-1]"
	print "P = ", P, "[dyne cm^-2]"
	print "T = ", T, "[K]"
	print ""

	beta 	= gasproperties.calculate_beta(P, T)
	rho 	= gasproperties.calculate_density(P * beta, T, mu)
	kappa 	= opacity.calculate_opacity(T, rho)

	return

# test_load1_and_load2()
