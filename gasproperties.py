import numpy as np
import astropy.constants as const


def calculate_mu(X,Y,Z=None):
	"""
		Calculates mean molecular weight
		Assumes:
			Fully ionized gas
			Composition of Hydrogen, Helium, Metals

		Inputs:
			X 	- Hydrogen mass fraction
			Y 	- Helium mass fraction
			[Z]	- Metals mass fraction
			  	- if Z omitted, assumes X+Y+Z=1
		  	Inputs must sum to 1

		Outputs:
			mu - mean molecular weight

		Warnings:
	"""

	if Z is None:
		Z = 1 - X - Y
	elif (X + Y + Z) != 1:
		print "Incorrect values for X,Y,Z. Must sum to 1"
		raise Exception('calculate_mu incorrect inputs')

	mu = np.power((X / (1./2)) + (Y / (4./3)) + (Z / 2.), -1.)

	return mu

def calculate_mu_e(X,Y,Z=None):
	"""
		Calculates mean molecular weight of electrons
		Assumes:
			Fully ionized gas
			Composition of Hydrogen, Helium, Metals

		Inputs:
			X 	- Hydrogen mass fraction
			Y 	- Helium mass fraction
			[Z]	- Metals mass fraction
			  	- if Z omitted, assumes X+Y+Z=1
		  	Inputs must sum to 1

		Outputs:
			mu_e - mean molecular weight of electrons

		Warnings:
	"""

	if Z is None:
		Z = 1 - X - Y
	elif (X + Y + Z) != 1:
		print "Incorrect values for X,Y,Z. Must sum to 1"
		raise Exception('calculate_mu incorrect inputs')

	mu_e = np.divide(2., (1+X)) #Eq. 4.30, Kippenhan Weigert Weiss 2nd ed.

	return mu_e

def calculate_rad_pressure(T):
	"""
		Calculates radiative pressure of photon gas at temperature T
		Assumes:
			

		Inputs:
			T - temperature [K]

		Outputs:
			P_rad [dyne/cm^2]

		Warnings:
	"""

	if T <= 0:
		print "Temperature must be greater than 0"
		raise Exception('calculate_rad_pressure incorrect inputs')

	a = 4 * const.sigma_sb.cgs.value / const.c.cgs.value

	P_rad = a * np.power(T, 4.) / 3

	return P_rad


def calculate_beta(P_total,T):
	"""
		Calculates beta, ratio between gas pressure and total pressure
		Assumes:
			

		Inputs:
			P_total	- pressure [dyne cm^-2]
			T 		- temperature [K]

		Outputs:
			beta (can range between 0 and 1)

		Warnings:
	"""

	if P_total <= 0 or T <= 0:
		print "Inputs must be greater than 0"
		raise Exception('calculate_beta incorrect inputs')

	P_rad = calculate_rad_pressure(T)

	beta = np.divide(P_total - P_rad,  P_total)

	return beta

def calculate_density(P_gas, T, mu):
	"""
		Calculates density using ideal gas law
		Assumes:
			ideal gas law


		Inputs:
			P_gas	- pressure [dyne cm^-2]
			T 		- temperature [K]
			mu 		- mean molecular weight

		Outputs:
			rho 	- density [g cm^-3]

		Warnings:
	"""

	if (P_gas <= 0 or T <= 0) or mu <= 0:
		print "Inputs must be greater than 0"
		raise Exception('calculate_density incorrect inputs \n P_gas = ' + \
			str(P_gas) + "\n T = " + str(T))

	R 	= const.k_B.cgs.value / const.m_p.cgs.value # ideal gas constant

	rho = np.divide(P_gas * mu, R * T)

	return rho
	
def calculate_del_rad(m, l, P, T, kappa):
	"""
		Calculates temperature gradient of radiative heat transfer
			del = d ln T / d ln P

		Assumes:
			mass, m as independent variable

		Inputs:
			(all local conditions at an enclosed mass, m)
			m			1xN array	- mass [g]
			l			1xN array	- luminosity  [ergs s^-1]
			P			1xN array	- total Pressure [dyne cm^-2]
			T			1xN array	- Temperature [K]
			kappa		1xN array	- opacity [cm^2 g^-1]


		Outputs:
			del_rad 	1xN array 	- radiative pressure gradient
				= d ln T / d ln P

		Warnings:
			This may not be your actual pressure gradient
				Need to do: min(del_rad, del_ad)
	"""
	del_rad = np.divide(3., 64 * np.pi * const.sigma_sb.cgs.value  * const.G.cgs.value) * \
		np.divide(kappa * l * P, m * np.power(T, 4.))

	return del_rad

def calculate_del_ad():
	"""
		Calculates temperature gradient of radiative heat transfer
			del = d ln T / d ln P

		Assumes:
			mass, m as independent variable
			fixed composition (ionization levels not changing)

		Inputs:
			(all local conditions at an enclosed mass, m)


		Outputs:
			del_ad 	1x1 float 	- adiabatic pressure gradient
				= d ln T / d ln P

		Warnings:
			This may not be your actual pressure gradient
				Need to do: min(del_rad, del_ad)
	"""
	return .4

def calculate_del(m, l, P, T, kappa):
	"""
		Calculates temperature gradient
			del = d ln T / d ln P

		Assumes:
			mass, m as independent variable

		Inputs:
			(all local conditions at an enclosed mass, m)
			m			1xN array	- mass [g]
			l			1xN array	- luminosity  [ergs s^-1]
			P			1xN array	- total Pressure [dyne cm^-2]
			T			1xN array	- Temperature [K]
			kappa		1xN array	- opacity [cm^2 g^-1]


		Outputs:
			del_actual 	1xN array 	- overall pressure gradient
				= d ln T / d ln P

		Warnings:

	"""

	del_rad 	= calculate_del_rad(m, l, P, T, kappa)
	del_ad 		= calculate_del_ad()

	del_actual 	= np.minimum(del_rad, del_ad)

	return del_actual

def solve_problem_5():
	P_total_a 	= 10**(16.85)
	T_a 		= 10**(7.55)
	X_a			= 0
	Y_a			= .98
	mu_a 		= calculate_mu(X_a, Y_a)
	P_gas_a 	= P_total_a - calculate_rad_pressure(T_a)
	beta_a 		= P_gas_a / P_total_a

	P_total_b 	= 10**(16.87)
	T_b 		= 10**(6.91)
	X_b			= .70
	Y_b 		= .28
	mu_b 		= calculate_mu(X_b, Y_b)
	P_gas_b 	= P_total_b - calculate_rad_pressure(T_b)
	beta_b 		= P_gas_b / P_total_b



	print "5a)"
	rho_a = calculate_density(P_gas_a, T_a, mu_a)
	print "rho_a: ", rho_a, "[g cm^-3]"
	print "beta_a: ", beta_a

	print "5b)"
	rho_b = calculate_density(P_gas_b, T_b, mu_b)
	print "rho_b: ", rho_b, "[g cm^-3]"
	print "beta_b: ", beta_b


	return

# solve_problem_5()