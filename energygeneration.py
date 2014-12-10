import numpy as np

def calculate_epsilon(T, rho, X, Y, Z):
	"""
		Calculates energy generation rate from all sources
		Assumes:
			See assumptions of individual components
				- epsilon_n: assume only H burning
				- epsilon_nu, epsilon_g: assume negligible


		Inputs:
			T			- Temperature [K]
			rho 		- Density [g  cm^-3]
			X			- Hydrogen mass fraction
			Y			- Helium mass fraction
			Z			- Metals mass fraction

		Outputs:
			epsilon 	- energy generation rate [erg g^-1 s^-1]

		Warnings:
	"""
	epsilon_n 	= calculate_epsilon_n(T, rho, X, Y, Z)
	epsilon_nu 	= calculate_epsilon_nu()
	epsilon_g 	= calculate_epsilon_g()

	epsilon 	= epsilon_n - epsilon_nu + epsilon_g

	return epsilon

def calculate_epsilon_n(T, rho, X, Y, Z):
	"""
		Calculates energy generation rate from the pp chain + CNO cycle

		Assumes:
			only H burning

		Inputs:
			T			- Temperature [K]
			rho 		- Density [g  cm^-3]
			X			- Hydrogen mass fraction
			Y			- Helium mass fraction
			Z			- Metals mass fraction

		Outputs:
			epsilon 	- energy generation rate [erg g^-1 s^-1]


		Warnings:
	"""
	epsilon_pp 	= calculate_epsilon_pp(rho, T, X)
	epsilon_cno = calculate_epsilon_CNO(rho,T,X,Z)

	epsilon 	= epsilon_pp + epsilon_cno

	return epsilon

def calculate_epsilon_nu():
	"""
	Energy losses due to neutrino free-streaming

	Will be ignored for this project
	"""
	return 0.

def calculate_epsilon_g():
	"""
	Energy generation/losses due to gravitational contraction

	Will be ignored for this project
	"""
	return 0.

def calculate_epsilon_pp(rho, T, X1, f11=1.):
	"""
		Calculates energy generation rate from the pp chain

		Assumes:
			uses Eq 18.63, Kippenhan Weigert Weiss 2nd ed.


		Inputs:
			rho 		- Density [g  cm^-3]
			T			- Temperature [K]
			X1			- Hydrogen mass fraction
			f11			- electron shielding factor [default = 1.; weak sheilding]

		Outputs:
			epsilon_pp 	- energy generation rate [erg g^-1 s^-1]


		Warnings:
			psi is approximate between T = 1e7 and 3e7 K
			exact coefficients might vary between sources
	"""

	T9 	= T/1e9

	g11 = 1 + 3.82 * T9 + 1.51 * np.power(T9,2.) + 0.144 * np.power(T9, 3.) - 0.0114 * np.power(T9, 4.)

	epsilon_pp = 2.57e4 * calculate_psi(T) * f11 * g11 * rho * np.power(X1, 2.) * np.power(T9, -2./3) * \
		np.exp(-3.381 * np.power(T9, -1./3))

	return epsilon_pp

def calculate_psi(T):
	"""
		Calculates correction for pp chain due to He abundances

		Assumes:
			- uses Fig 18.7, Kippenhan Weigert Weiss 2nd ed.
			- ramp approximation for between T7=1 and T7=3.


		Inputs:
			T		- Temperature [K]

		Outputs:
			psi 	- energy generation correction factor


		Warnings:
			works poorly (within a factor of 2) for high He abundances
			very approximate
	"""

	T7 = T / 1e7

	if 	 T7 < 1.:
		psi = 1.
	elif T7 > 3.:
		psi = 2.
	else:
		psi = 1 + (T7-1)*0.5 

	return psi

def calculate_epsilon_CNO(rho, T, X1, Z):
	"""
		Calculates energy generation rate from the CNO chain

		Assumes:
			uses Eq 18.65, Kippenhan Weigert Weiss 2nd ed.


		Inputs:
			rho 		- Density [g  cm^-3]
			T			- Temperature [K]
			X1			- Hydrogen mass fraction
			Z			- metals mass fraction

		Outputs:
			epsilon_cno - energy generation rate [erg g^-1 s^-1]


		Warnings:
			exact coefficients might vary between sources
	"""

	T9 		= T/1e9

	XCNO 	= 0.7 * Z # approximate

	g141 	= 1 - 2.00 * T9 + 3.41 * np.power(T9,2.) - 2.43 * np.power(T9, 3.)

	epsilon_cno = 8.24e25 * g141 * XCNO * X1 * rho * np.power(T9, -2./3) * \
		np.exp((-15.231 * np.power(T9, -1./3)) - np.power(T9/0.8, 2))

	return epsilon_cno

def test_epsilons():
	X = .70
	Y = .28
	Z = .02

	T = 18e6
	rho = 80

	ratio = calculate_epsilon_pp(rho, T, X) / calculate_epsilon_CNO(rho, T, X, Z)

	print "epsilon_pp: ", calculate_epsilon_pp(rho, T, X)
	print "epsilon_CNO: ", calculate_epsilon_CNO(rho, T, X, Z)
	print ratio

	return


# test_epsilons()