import numpy as np
from scipy import interpolate
from itertools import izip


def read_in_opacity_data(filename):
	"""
		Reads in and parses the OPAL opacity table

		Assumes:
			Rectangular grid
				Missing data MUST be marked "N/A"

		Inputs:
			filename - string to 

		Outputs:
			log_T 		1xN array	- Temperature [K]
			log_R 		1xM array	- rho / T6^3 [g  cm^-3  K^-3]
			log_kappa 	NxM array 	- opacity [cm^2/g]

		Warnings:
			R is not density, instead:
				R = rho / T6^3
			Will probably contain np.nan values
	"""

	# print filename
	opacity_raw = np.genfromtxt(filename, skiprows=4, filling_values=np.nan, 
		missing_values='N/A')

	# print opacity_raw
	# print opacity_raw[-1][-1] + 1

	log_T = opacity_raw[1:,0]
	log_R = opacity_raw[0,1:]

	log_kappa = opacity_raw[1:,1:]


	return log_T, log_R, log_kappa


def interpolate_opacity_table(log_T, log_R, log_kappa):
	"""
		Takes a sampled grid of opacities, and interpolates between values

		Assumes:
			

		Inputs:
			log_T 		1xN array	- Temperature [K]
			log_R 		1xM array	- rho / T6^3 [g  cm^-3  K^-3]
			log_kappa 	NxM array 	- opacity [cm^2/g]

		Outputs:
			opacity_interp_fn	- pointer to interpolation function
								- log_opacity = f(log_T, log_R)

		Warnings:
			R is not density, instead:
				R = rho / T6^3
	"""
	R 	= np.power(10, log_R)
	T6	= np.power(10, log_T-6)

	rho_grid 		= np.outer(np.power(T6,3), R)
	log_rho_grid	= np.log10(rho_grid)
	log_T_grid		= np.transpose(np.tile(log_T, (len(R),1)))

	log_T_flat 		= log_T_grid.flatten()
	log_rho_flat 	= log_rho_grid.flatten()
	log_kappa_flat 	= log_kappa.flatten()

	in_bounds 		= np.isfinite(log_kappa_flat) #filter out nan's

	log_T_flat 		= log_T_flat[in_bounds]
	log_rho_flat 	= log_rho_flat[in_bounds]
	log_kappa_flat 	= log_kappa_flat[in_bounds]

	opacity_interp_fn = interpolate.interp2d(log_T_flat, log_rho_flat, log_kappa_flat,
						kind='cubic', bounds_error=True)

	# # # Check against metallicity table, to ensure indices are working right:
	# print log_T_flat[10]
	# print calculate_log_R(log_T_flat[10], log_rho_flat[10])
	# print opacity_interp_fn(log_T_flat[10], log_rho_flat[10])


	return opacity_interp_fn


def calculate_log_R(log_T, log_rho):
	"""
		Calculates log_R, from log_T and log_rho
			Useful for checking functions during development

		Assumes:
			

		Inputs:
			log_T	- Temperature [K]
			log_rho - Density [g  cm^-3]

		Outputs:
			log_R 	- rho / T6^3 [g  cm^-3  K^-3]


		Warnings:
			R is not density, instead:
				R = rho / T6^3
	"""
	return (log_rho - (3 * (log_T-6)))

def solve_question_6():
	filename = "opacity_data/opacity_solar_like_metallicity.dat"
	# X = .70
	# Y = .28
	# Z = .02

	log_T, log_R, log_kappa = read_in_opacity_data(filename)

	opacity_interp_fn_log	= interpolate_opacity_table(log_T, log_R, log_kappa)


	print 'Part a)'
	log_T_a = 6.3
	log_rho_a = 0.3
	log_R_a = calculate_log_R(log_T_a, log_rho_a) 

	print 'log_T = ', log_T_a, '[K]'
	print 'log_rho = ', log_rho_a, '[g cm^-3]'
	print 'log_R = ', log_R_a
	print 'log_kappa = ', opacity_interp_fn_log(log_T_a, log_rho_a), 'cm^2 g^-1'
	print 'kappa = ', np.power(10,opacity_interp_fn_log(log_T_a, log_rho_a)), 'cm^2 g^-1'
	print ''


	print 'Part b)'
	log_T_b = 5
	log_rho_b = -4
	log_R_b = calculate_log_R(log_T_b, log_rho_b) 

	print 'log_T = ', log_T_b, '[K]'
	print 'log_rho = ', log_rho_b, '[g cm^-3]'
	print 'log_R = ', log_R_b
	print 'log_kappa= ', opacity_interp_fn_log(log_T_b, log_rho_b), 'cm^2 g^-1'
	print 'kappa = ', np.power(10,opacity_interp_fn_log(log_T_b, log_rho_b)), 'cm^2 g^-1'


	return


def get_opacity_interp_fn():

	filename = "opacity_data/opacity_solar_like_metallicity.dat"
	# X = .70
	# Y = .28
	# Z = .02

	log_T, log_R, log_kappa = read_in_opacity_data(filename)

	opacity_interp_fn_log	= interpolate_opacity_table(log_T, log_R, log_kappa)

	def opacity_interp_fn(Ts, rhos, opacity_interp_fn_log=opacity_interp_fn_log):
		"""
			For use in producing an output closure for get_opacity_interp_fn

			Assumes:
				Inputs will either have same length
					or one will have length 1, with that input being repeated

			Inputs:
				Ts			1xN float	- Temperature [K]
				rhos		1xM float	- density [g cm^-3]
					where (N==M) or N==1 or M==1

				opacity_interp_fn_log - function reference
					- e.g. log10_kappa = opacity_interp_fn_log(log10_T, log10_rho)
					- produced by interpolate_opacity_table

			Outputs:
				opacity_interp_fn 	- function handle
					- e.g.  kappa = opacity_interp_fn(T, rho)


			Warnings:
				If input sizes don't match the specifications, 
					UserWarning will be thrown

				Multi-dimensional arrays might result in unexpected behavior

				Intended to be visible from outside modules
		"""
		if np.size(Ts) is 1 and np.size(rhos) is 1:
			return np.float(np.power(10., 
				opacity_interp_fn_log(np.log10(Ts), np.log10(rhos))))
		elif np.size(Ts) is 1:
			return np.array([np.power(10., 
				opacity_interp_fn_log(np.log10(Ts), np.log10(rho))) 
				for rho in rhos]).flatten()	
		elif np.size(rhos) is 1:
			return np.array([np.power(10., 
				opacity_interp_fn_log(np.log10(T), np.log10(rhos))) 
				for T in Ts]).flatten()	
		else:
			raise UserWarning('inputs of unequal size; must match or one must be size 1')

	return opacity_interp_fn

calculate_opacity = get_opacity_interp_fn()


# solve_question_6()
