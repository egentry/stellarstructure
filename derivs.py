import numpy as np
import astropy.constants as const
import gasproperties
import loadinitial
import opacity
import energygeneration
import modelparameters


def calculate_dr_dm(r, rho):
	dr_dm = np.divide(1., 4 * np.pi * np.power(r,2.) * rho )
	return dr_dm

def calculate_dl_dm(T, rho, X, Y, Z):
	dl_dm = energygeneration.calculate_epsilon(T, rho, X, Y, Z)
	return dl_dm

def calculate_dP_dm(m, r):
	dP_dm = np.divide(-1. * const.G.cgs.value * m, 
		4 * np.pi * np.power(r,4.))
	return dP_dm

def calculate_dT_dm(m, r, l, P, T, kappa):
	del_actual 	= gasproperties.calculate_del(m, l, P, T, kappa)

	dT_dm = np.divide(-1 * const.G.cgs.value * m * T * del_actual, 
		4 * np.pi * np.power(r, 4.) * P )

	return dT_dm



def calculate_all_derivs(m, r, l, P, T, rho, kappa, X, Y, Z):
	dr_dm = calculate_dr_dm(r, rho)
	dl_dm = calculate_dl_dm(T, rho, X, Y, Z)
	dP_dm = calculate_dP_dm(m, r)
	dT_dm = calculate_dT_dm(m, r, l, P, T, kappa)

	return [dr_dm, dl_dm, dP_dm, dT_dm]

def print_inner_derivs(mass_scale):
	X 		= modelparameters.X
	Y 		= modelparameters.Y
	Z 		= modelparameters.Z
	mu 		= modelparameters.mu

	P_c 			= modelparameters.P_c # core pressure, [dyne cm^-2]
	T_c 			= modelparameters.T_c # core temperature, [K]

	m = 1e-6

	r, l, P, T 	= loadinitial.load1(m, P_c, T_c, mu, X, Y, Z)
	beta 		= gasproperties.calculate_beta(P_c, T_c)
	rho 		= gasproperties.calculate_density(P_c * beta, T_c, mu)
	kappa 		= opacity.calculate_opacity(T, rho)

	dr_dm, dl_dm, dP_dm, dT_dm = calculate_all_derivs(m, r, l, P, T, rho, kappa, X, Y, Z)

	# print ""
	# print "core derivatives: "
	# print "dr/dm = ", dr_dm, "[cm g^-1]"
	# print "dl/dm = ", dl_dm, "[ergs s^-1 g^-1]"
	# print "dP/dm = ", dP_dm, "[dyne cm^-2 g^-1]"
	# print "dT/dm = ", dT_dm, "[K g^-1]"
	# print ""

	# print ""
	# print "core values: "
	# print "m = ", m, "[g]"
	# print "r = ", r, "[cm]"
	# print "l = ", l, "[ergs s^-1]"
	# print "P = ", P, "[dyne cm^-2]"
	# print "T = ", T, "[K]"
	# print "rho = ", rho, "[g cm^-3]"
	# print "del_rad = ", gasproperties.calculate_del_rad(m,l,P,T,kappa), "[cm^2 g^-1]"
	# print "del_ad = ", gasproperties.calculate_del_ad(), "[cm^2 g^-1]"
	# print "del = ", gasproperties.calculate_del(m,l,P,T,kappa), "[cm^2 g^-1]"
	# print ""

	return

def print_outer_derivs():
	M_star 			= modelparameters.M_star
	R_star 			= modelparameters.R_star
	L_star 			= modelparameters.L_star
	g = const.G.cgs.value * M_star / (R_star**2.)

	X 		= modelparameters.X
	Y 		= modelparameters.Y
	Z 		= modelparameters.Z
	mu 		= modelparameters.mu

	r, l, P, T = loadinitial.load2(M_star, R_star, L_star, mu)
	beta 		= gasproperties.calculate_beta(P, T)
	rho = gasproperties.calculate_density(P * beta, T, mu)
	kappa = opacity.calculate_opacity(T,rho)

	dr_dm, dl_dm, dP_dm, dT_dm = calculate_all_derivs(M_star, r, l, P, T, rho, kappa, X, Y, Z)

	print ""
	print "surface derivatives: "
	print "dr/dm = ", dr_dm, "[cm g^-1]"
	print "dl/dm = ", dl_dm, "[ergs s^-1 g^-1]"
	print "dP/dm = ", dP_dm, "[dyne cm^-2 g^-1]"
	print "dT/dm = ", dT_dm, "[K g^-1]"
	print ""

	return


# for i in xrange(100):
# 	print_inner_derivs(1e6)
# print_outer_derivs()
