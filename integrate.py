import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import derivs
import gasproperties
import opacity
import loadinitial
import modelparameters
from multiprocessing import Pool


def deriv_wrapper(y, x, X, Y, Z, mu):
	"""
		Creates a system of 1st order ode's to be solved

		Assumes:
			mass, m as independent variable

		Inputs:
			y 			1x4 float 	- consists of:
										- r(m), radius [cm]
										- l(m), luminosity [erg s^-1]
										- P(m), total pressure [dyne cm^-2]
										- T(m), Temperature [K]

			x 			1x1 float 	- consists of:
										- m, mass [g]

			params 		1x4 float 	- consists of:
										- X, hydrogen mass fraction
										- Y, helium mass fraction
										- Z, metals mass fraction
										- mu, mean molecular weight

		Outputs:
			dy_dx		1x4 float 	- consists of:
										- dr(m)/dm, radius derivative [cm g^-1]
										- dl(m)/dm, luminosity derivative [erg s^-1 g^-1]
										- dP(m)/dm, total pressure derivative [dyne cm^-2 g^-1]
										- dT(m)/dm, Temperature derivative [K g^-1]

		Warnings:

	"""
	m 			= x
	r, l, P, T 	= y

	beta 		= gasproperties.calculate_beta(P, T)
	rho 		= gasproperties.calculate_density(P * beta, T, mu)
	kappa 		= opacity.calculate_opacity(T, rho)

	dr_dm		= derivs.calculate_dr_dm(r, rho)
	dl_dm		= derivs.calculate_dl_dm(T, rho, X, Y, Z)
	dP_dm		= derivs.calculate_dP_dm(m, r)
	dT_dm		= derivs.calculate_dT_dm(m, r, l, P, T, kappa)

	dy_dx 		= [dr_dm, dl_dm, dP_dm, dT_dm]

	return dy_dx


def integrate_outwards(M_star, m_fitting_point, P_c, T_c, mu, X, Y, Z,
		n_steps=1e4,
		logspacing=True,
		file_suffix="",
		write=False):
	m0 				= 1e-8 * M_star

	beta 			= gasproperties.calculate_beta(P_c, T_c)
	rho 			= gasproperties.calculate_density(P_c * beta, T_c, mu)
	r0, l0, P0, T0 	= loadinitial.load1(m0, P_c, T_c, mu, X, Y, Z)
	y0 				= [r0, l0, P0, T0]

	mu 		= modelparameters.mu
	params 	= (X, Y, Z, mu)

	if logspacing is True:
		m = np.logspace(np.log10(m0), np.log10(m_fitting_point), n_steps)
	else:
		m = np.linspace(m0, m_fitting_point, n_steps)
	y, infodict = integ.odeint(deriv_wrapper, y0, m, 
		mxstep=500,
		args=params, full_output=True)
	r,l,P,T 	= y.transpose()
	sol 		= np.column_stack((m, y))
	if write is True:
		np.savetxt('data/sol_outwards' + file_suffix + '.dat', sol, 
			header=" \t\t m [m]\t\t\t\t\t r [cm]\t\t\t\t\t\t l [erg s^-1]\t\t\t\t\t P [dyne cm^-2]\t\t\t\t\t\t T [K]")


		plt.figure(1)

		plt.subplot(221)
		plt.plot(m / M_star, r)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$r(m)$")

		plt.subplot(222)
		plt.semilogy(m / M_star, l)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$\ell(m)$")

		plt.subplot(223)
		plt.semilogy(m / M_star, P)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$P(m)$")

		plt.subplot(224)
		plt.semilogy(m / M_star, T)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$T(m)$")


		plt.savefig("plots/stellar_model_outwards" + file_suffix + ".eps")
		plt.savefig("plots/stellar_model_outwards" + file_suffix + ".pdf")
		# plt.show()
		plt.close()

	return m, y, infodict

def integrate_inwards(M_star, m_fitting_point, R_star, L_star, mu, X, Y, Z,
	n_steps=1e4,
	logspacing=False,
	file_suffix="",
	write=False):

	r0, l0, P0, T0 	= loadinitial.load2(M_star, R_star, L_star, mu)
	y0 				= [r0, l0, P0, T0]
	mu 		= modelparameters.mu
	params 	= (X, Y, Z, mu)

	if logspacing is True:
		m = np.logspace(np.log10(m_fitting_point), np.log10(M_star), n_steps)
	else:
		m = np.linspace(m_fitting_point, M_star, n_steps)
	m 			= np.flipud(m) #reverse direction of integration
	y, infodict = integ.odeint(deriv_wrapper, y0, m, 
		mxstep=5000,
		args=params, full_output=True)
	r,l,P,T 	= y.transpose()
	sol 		= np.column_stack((m, y))

	if write is True:
		np.savetxt('data/sol_inwards' + file_suffix + '.dat', sol, 
			header=" \t\t m [m]\t\t\t\t\t r [cm]\t\t\t\t\t\t l [erg s^-1]\t\t\t\t\t P [dyne cm^-2]\t\t\t\t\t\t T [K]")


		plt.figure(1)

		plt.subplot(221)
		plt.plot(m / M_star, r)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$r(m)$")

		plt.subplot(222)
		plt.semilogy(m / M_star, l)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$\ell(m)$")

		plt.subplot(223)
		plt.semilogy(m / M_star, P)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$P(m)$")

		plt.subplot(224)
		plt.semilogy(m / M_star, T)
		plt.xlabel(r"$\frac{m}{M}$")
		plt.ylabel(r"$T(m)$")


		plt.savefig("plots/stellar_model_inwards" + file_suffix + ".pdf")
		# plt.show()
		plt.close()

	return m, y, infodict


def test():
	X 		= modelparameters.X
	Y 		= modelparameters.Y
	Z 		= modelparameters.Z
	mu 		= modelparameters.mu
	params 	= (X, Y, Z, mu)

	P_c 			= modelparameters.P_c # core pressure, [dyne cm^-2]
	T_c 			= modelparameters.T_c # core temperature, [K]
	M_star 			= modelparameters.M_star
	R_star 			= modelparameters.R_star
	L_star 			= modelparameters.L_star
	m_fitting_point	= modelparameters.m_fitting_point




	m_outward, y_outward, infodict_outward 		= integrate_outwards(M_star, 
		m_fitting_point, P_c, T_c, mu, X, Y, Z, n_steps = 5e1)

	m_inward, y_inward, infodict_inward 		= integrate_inwards(M_star, 
		m_fitting_point, R_star, L_star, mu, X, Y, Z, n_steps = 5e1)

	r_inward,	l_inward,	P_inward,	T_inward 	= y_inward.transpose()
	r_outward,	l_outward,	P_outward,	T_outward 	= y_outward.transpose()

	m_tot = np.concatenate((m_outward, np.flipud(m_inward)))
	r_tot = np.concatenate((r_outward, np.flipud(r_inward)))
	l_tot = np.concatenate((l_outward, np.flipud(l_inward)))
	P_tot = np.concatenate((P_outward, np.flipud(P_inward)))
	T_tot = np.concatenate((T_outward, np.flipud(T_inward)))



	plt.figure(1)

	plt.subplot(221)
	plt.plot(m_tot / M_star, r_tot)
	plt.xlabel(r"$\frac{m}{M}$")
	plt.ylabel(r"$r(m)$")

	plt.subplot(222)
	plt.semilogy(m_tot / M_star, l_tot)
	plt.xlabel(r"$\frac{m}{M}$")
	plt.ylabel(r"$\ell(m)$")

	plt.subplot(223)
	plt.semilogy(m_tot / M_star, P_tot)
	plt.xlabel(r"$\frac{m}{M}$")
	plt.ylabel(r"$P(m)$")

	plt.subplot(224)
	plt.semilogy(m_tot / M_star, T_tot)
	plt.xlabel(r"$\frac{m}{M}$")
	plt.ylabel(r"$T(m)$")


	plt.savefig("plots/stellar_model_total.pdf")
	# plt.show()
	plt.close()

	return (m_tot, r_tot, l_tot, P_tot, T_tot)
