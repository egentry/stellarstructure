import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import modelparameters
import integrate
from multiprocessing import Pool
from scipy import optimize


def get_fractional_errors(R_star, L_star, P_c, T_c):
	"""
		Pass in "guess" conditions.
		Will then calculate inward and outward errors,

		Returns:
			[Data array]
			dY - over/undershoots (+/-, going outward)
				[dx handled outside this]
	"""

	# R_star, L_star, P_c, T_c = x

	P_c_0		= modelparameters.P_c # core pressure, [dyne cm^-2]
	T_c_0 		= modelparameters.T_c # core temperature, [K]
	R_star_0 	= modelparameters.R_star
	L_star_0 	= modelparameters.L_star

	print ""
	print "R: " + str(R_star / R_star_0)
	print "L: " + str(L_star / L_star_0)
	print "P: " + str(P_c / P_c_0)
	print "T: " + str(T_c / T_c_0)


	X 		= modelparameters.X
	Y 		= modelparameters.Y
	Z 		= modelparameters.Z
	mu 		= modelparameters.mu
	params 	= (X, Y, Z, mu)

	M_star 	= modelparameters.M_star
	m_fitting_point	= modelparameters.m_fitting_point

	pool = Pool(2)
	outward_results = pool.apply_async(integrate.integrate_outwards, 
		[M_star, m_fitting_point, P_c, T_c, mu, X, Y, Z] )

	inward_results  = pool.apply_async(integrate.integrate_inwards, 
		[M_star, m_fitting_point, R_star, L_star, mu, X, Y, Z] )

	m_outward, y_outward, infodict_outward 	= outward_results.get()

	m_inward, y_inward, infodict_inward 	= inward_results.get()

	dr = y_inward[-1,0] - y_outward[-1,0]
	dl = y_inward[-1,1] - y_outward[-1,1]
	dP = y_inward[-1,2] - y_outward[-1,2]
	dT = y_inward[-1,3] - y_outward[-1,3]

	dY = np.array([dr, dl, dP, dT])

	print ''
	print 'fractional errors:'
	print "dR: " + str(dr / y_inward[-1,0])
	print "dL: " + str(dl / y_inward[-1,1])
	print "dP: " + str(dP / y_inward[-1,2])
	print "dT: " + str(dT / y_inward[-1,3])

	return dY

def get_new_guess(R_star, L_star, P_c, T_c, step_scale = .04):	
		dR = R_star * step_scale
		dL = L_star * step_scale
		dP = P_c 	* step_scale
		dT = T_c 	* step_scale


		dY = get_fractional_errors(R_star, L_star, P_c, T_c)
		if np.max(np.absolute(dY)) < .01:
			return np.zeros(4), dY

		dY_dR = (dY - get_fractional_errors(R_star + dR, L_star, P_c, T_c)) / step_scale
		dY_dL = (dY - get_fractional_errors(R_star, L_star + dL, P_c, T_c)) / step_scale
		dY_dP = (dY - get_fractional_errors(R_star, L_star, P_c + dP, T_c)) / step_scale
		dy_dT = (dY - get_fractional_errors(R_star, L_star, P_c, T_c + dT)) / step_scale

		partial_deriv_matrix = np.column_stack((dY_dR, dY_dL, dY_dP, dy_dT))
		partial_deriv_matrix_inv = np.linalg.inv(partial_deriv_matrix)

		dX = np.dot(partial_deriv_matrix_inv, dY)

		print ''
		print "Proposed change:"
		print dX
		print "dR: " + str(dX[0])
		print "dL: " + str(dX[1])
		print "dP: " + str(dX[2])
		print "dT: " + str(dX[3])

		return dX, dY



M_star		= modelparameters.M_star
P_c_0		= modelparameters.P_c # core pressure, [dyne cm^-2]
T_c_0 		= modelparameters.T_c # core temperature, [K]
R_star_0 	= modelparameters.R_star
L_star_0 	= modelparameters.L_star

dX_scale = .25 # helps prevent overshooting into unstable regions

dX, dY = get_new_guess(R_star_0, L_star_0, P_c_0, T_c_0)
R_star 	= R_star_0 	* (1 + dX[0]*dX_scale)
L_star 	= L_star_0 	* (1 + dX[1]*dX_scale)
P_c 	= P_c_0 	* (1 + dX[2]*dX_scale)
T_c 	= T_c_0 	* (1 + dX[3]*dX_scale)
	
dX, dY = get_new_guess(R_star_0 * (1 -.02), L_star_0 * (1-.2), 
	P_c_0 * (1-.04), T_c_0 * (1-.01))
i=0
while (np.sum(np.equal(dX, 0)) is not 4) and (i < 50) :
	i = i + 1
	dX, dY = get_new_guess(R_star, L_star, P_c, T_c)

	R_star 	= R_star 	* (1 + dX[0]*dX_scale)
	L_star 	= L_star 	* (1 + dX[1]*dX_scale)
	P_c 	= P_c 		* (1 + dX[2]*dX_scale)
	T_c 	= T_c 		* (1 + dX[3]*dX_scale)

print ''
print "success!"
print 'iterations: ' + str(i)
print "dX: " + str(dX)
print "Delta R: " + str(R_star)
print "Delta L: " + str(L_star)
print "Delta Tc: " + str(T_c)
print "Delta Pc: " + str(P_c)

X 		= modelparameters.X
Y 		= modelparameters.Y
Z 		= modelparameters.Z
mu 		= modelparameters.mu


M_sol 	= modelparameters.M_solar
R_sol 	= modelparameters.R_solar
L_sol 	= modelparameters.L_solar

# R_star 	= R_star_0 	* 1.14651564384
# L_star 	= L_star_0 	* 0.829448956503
# P_c 	= P_c_0 	* 0.96758319806
# T_c 	= T_c_0 	* 0.970207578518

print 'initial'
print 'R_star: ' + str(R_star_0/R_sol)
print 'L_star: ' + str(L_star_0/L_sol)
print 'P_c: ' + str(np.log10(P_c_0))
print 'T_c: ' + str(np.log10(T_c_0))

print 'converged:'
print 'R_star: ' + str(R_star/R_sol)
print 'L_star: ' + str(L_star/L_sol)
print 'P_c: ' + str(np.log10(P_c))
print 'T_c: ' + str(np.log10(T_c))





m_fitting_point	= modelparameters.m_fitting_point

pool = Pool(2)
results_outwards = pool.apply_async(integrate.integrate_outwards, [M_star, 
	m_fitting_point, P_c_0, T_c_0, mu, X, Y, Z], {"n_steps":1e4, "write":True, "file_suffix":"_final"})

results_inwards = pool.apply_async(integrate.integrate_inwards, [M_star, 
	m_fitting_point, R_star_0, L_star_0, mu, X, Y, Z], {"n_steps":1e4, "write":True, "file_suffix":"_final"})



m_outward, y_outward, infodict_outward 		= results_outwards.get()

m_inward, y_inward, infodict_inward 		= results_inwards.get()

r_inward,	l_inward,	P_inward,	T_inward 	= y_inward.transpose()
r_outward,	l_outward,	P_outward,	T_outward 	= y_outward.transpose()

m_tot = np.concatenate((m_outward, np.flipud(m_inward)))
r_tot = np.concatenate((r_outward, np.flipud(r_inward)))
l_tot = np.concatenate((l_outward, np.flipud(l_inward)))
P_tot = np.concatenate((P_outward, np.flipud(P_inward)))
T_tot = np.concatenate((T_outward, np.flipud(T_inward)))


sol_tot 		= np.column_stack((m_tot, r_tot, l_tot, P_tot, T_tot))

np.savetxt('data/sol_final.dat', sol_tot, 
	header=" \t\t m [m]\t\t\t\t\t r [cm]\t\t\t\t\t\t l [erg s^-1]\t\t\t\t\t P [dyne cm^-2]\t\t\t\t\t\t T [K]")



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


plt.savefig("plots/stellar_model_final.pdf")
# plt.show()
plt.close()





























