import numpy as np
import gasproperties
import astropy.constants as const
import gasproperties

M_solar = 1.988e33 	# [g]
R_solar = 6.955e10 	# [cm]
L_solar = 3.83e33 	# [erg s^-1]

mass_scale = 1.5

# P_c 		= 2.2e17 # core pressure, [dyne cm^-2]
# T_c 		= 1.5e7 # core temperature, [K]

# M_star 		= mass_scale * M_solar
# R_star 		= (mass_scale**.57) * R_solar
# L_star 		= (mass_scale**4.)	* L_solar

#using Hansen + Kawaler + Trimble, Table 2.5 + 2.6
M_star 		= 2 * M_solar
L_star 		= L_solar * np.power(10., 1.262)
T_eff		= np.power(10., 3.992)
R_star 		= np.power(np.divide(L_star, 
	4 * np.pi  * const.sigma_sb.cgs.value * T_eff**4.), .5)

T_c 		= 21.09e6
P_c			= np.power(10., 17.21)

adhoc_P_c_scale = .75
adhoc_T_c_scale = 1.25
adhoc_R_star_scale = 1.25
adhoc_L_star_scale = 1.

# R_star 	= R_star * adhoc_R_star_scale
# L_star 	= L_star * adhoc_L_star_scale
# P_c 	= P_c  * adhoc_P_c_scale
# T_c 	= T_c * adhoc_T_c_scale

# print T_c
# print P_c

m_fitting_point = M_star * .5


X 		= .70
Y 		= .28
Z 		= .02
mu 		= gasproperties.calculate_mu(X,Y,Z)
mu_e	= gasproperties.calculate_mu_e(X,Y,Z)

