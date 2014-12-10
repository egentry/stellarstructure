import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import modelparameters


filename_ezweb		= 'data/ezweb/ezweb_07067/structure_00000.txt'
filename_initial	= 'data/sol_initial.dat'
filename_converged	= 'data/sol_final.dat'




m_ezweb, r_ezweb, l_ezweb, P_ezweb, rho_ezweb, T_ezweb = np.loadtxt(filename_ezweb,
	unpack=True,
	usecols=(0,1,2,3,4,5))
M_star_ezweb			= m_ezweb[0]
R_star_ezweb			= r_ezweb[0]
L_star_ezweb			= l_ezweb[0]

# print m_ezweb
# print r_ezweb


m_initial, r_initial, l_initial, P_initial, T_initial = np.loadtxt(filename_initial, skiprows=1, unpack=True)
M_star_initial			= m_initial[-1]
R_star_initial			= r_initial[-1]
L_star_initial			= l_initial[-1]

m_converged, r_converged, l_converged, P_converged, T_converged = np.loadtxt(filename_converged,
	skiprows=1, unpack=True)
M_star_converged		= m_converged[-1]
R_star_converged		= r_converged[-1]
L_star_converged		= l_converged[-1]

M_star 			= M_star_converged
m_fitting_point	= modelparameters.m_fitting_point


plt.figure(1)

plt.subplot(221)
plt.loglog(m_converged / M_star_converged, r_converged / R_star_converged)
plt.loglog(m_ezweb / M_star_ezweb, r_ezweb / R_star_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.ylim(1e-3, 3)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$r / R_\ast$")

plt.subplot(222)
plt.loglog(m_converged / M_star_converged, l_converged/ L_star_converged)
plt.loglog(m_ezweb / M_star_ezweb, l_ezweb / L_star_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.ylim(1e-5, 3)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$\ell / L_\ast$")

plt.subplot(223)
plt.loglog(m_converged / M_star_converged, P_converged)
plt.loglog(m_ezweb / M_star_ezweb, P_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$P$ [dyne cm$^{-2}$]")

plt.subplot(224)
plt.loglog(m_converged / M_star_converged, T_converged)
plt.loglog(m_ezweb / M_star_ezweb, T_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$T$ [K]")

plt.tight_layout()
plt.savefig("plots/stellar_model_converged_good.eps", bbox='tight')
plt.savefig("plots/stellar_model_converged_good.pdf", bbox='tight')
# plt.show()
plt.close()



plt.figure(1)

plt.subplot(221)
plt.loglog(m_initial / M_star_initial, r_initial / R_star_initial)
# plt.loglog(m_ezweb / M_star_ezweb, r_ezweb / R_star_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.ylim(1e-3, 3)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$r / R_\ast$")

plt.subplot(222)
plt.loglog(m_initial / M_star_initial, l_initial/ L_star_initial)
# plt.loglog(m_ezweb / M_star_ezweb, l_ezweb / L_star_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.ylim(1e-5, 3)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$\ell / L_\ast$")

plt.subplot(223)
plt.loglog(m_initial / M_star_initial, P_initial)
# plt.loglog(m_ezweb / M_star_ezweb, P_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$P$ [dyne cm$^{-2}$]")

plt.subplot(224)
plt.loglog(m_initial / M_star_initial, T_initial)
# plt.loglog(m_ezweb / M_star_ezweb, T_ezweb, linestyle='dashed')
plt.axvline(x= m_fitting_point/M_star, linestyle='dotted')
plt.xlim(1e-6, 1)
plt.xlabel(r"$\frac{m}{M_\ast}$")
plt.ylabel(r"$T$ [K]")

plt.tight_layout()
plt.savefig("plots/stellar_model_initial_good.eps", bbox='tight')
plt.savefig("plots/stellar_model_initial_good.pdf", bbox='tight')
# plt.show()
plt.close()