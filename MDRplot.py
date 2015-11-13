# 
# 	Short script for calculations associated with
#	Morphology dependent resonances.
#	The theory is based on:
#
#	[1] Johnson; JOSA A/Vol. 10, No. 2, (Feb. 1993)
#
#	Author: Patrick Stegmann (PST)
#	Date: 2015-11-02
#
#	Revisions:
#	----------
#
#	- Radial potential plot (2015-11-02,PST)
#	- Radial wave function evaluation & plot (2015-11-11,PST)
#	- Numerical solution to check semi-analytic approach (2015-11-12,PST)
#	- Print number and location of MDR resonances (2015-11-13,PST)
#
#	Programmer's Comments:
#	----------------------
#	
#	1. The Bessel functions from Scipy always return a
#		tuple with shape [1][n]. Undocumented "feature".
# 

# import modules
import numpy as np 				# NumPy
import scipy as sp         		# Scipy

# define constants
global m,n,k,a
m = 1.47 						# refractive index
n = 40 							# radial mode order
k = 35.286219 					# wave number == resonace value for semianalytic approach
a = 1.0 						# particle radius
r = np.linspace(0.1,3.0,3000) 	# radial coordinate

# calculate potential
def potential(radi):
	if radi <= a:
		V = k*k*(1-m**2) + n*(n-1)/(radi*radi)
	elif radi > a:
		V = n*(n-1)/(radi*radi)
	else:
		V = 1.0
	return V

V = np.zeros(r.shape,dtype=float)

vecV = np.vectorize(potential)
V = vecV(r)
#for ii in range(1,r.size):
	#V[ii] = potential(r[ii])

import matplotlib as mpl 		# Matplotlib (2D/3D)
import matplotlib.pyplot as plt # Matplotlib's pyplot
from pylab import * 			# Matplotlib's pylab
fig0 = figure(0)
plt.plot(r,V/np.amax(V),'k-',linewidth=2)
axes = plt.gca()
axes.set_xlim([0.5,2.0])
axes.set_ylim([0.,0.035])
plt.text(a, 1800, 'scattering particle radius')
plt.text(0.41, k*k, '$k^2$',fontsize=15)
plt.xlabel('radial coordinate $r$',fontsize=20)
plt.ylabel('potential $V_n(r)$', fontsize=20)
plt.grid('on'), plt.box('on')
fig0.show()

#=====================================================================================================
# Semianalytical approach towards the radial ODE.

from scipy import special		# Scipy's special functions

# radial wave function TE (eq. 12)
radwavefunc = np.zeros(r.shape,dtype=float)

# First derivative of the Riccati-Bessel function
# of the first kind based on the spherical Bessel
# function.
def Psiprime(n,x):
	def jp(n,x):
		jderiv = sp.special.sph_jn(n,x)[0][n-1] \
				-(n+1)/x*sp.special.sph_jn(n,x)[0][n]
		return jderiv
	pp = sp.special.sph_jn(n,x)[0][n]+x*jp(n,x)
	return pp

# First derivative of the Riccati-Bessel function
# of the second kind based on the corresponding
# spherical Bessel function.
def Xiprime(n,x):
	def yp(n,x):
		yderiv = sp.special.sph_yn(n,x)[0][n-1] \
				-(n+1)/x*sp.special.sph_yn(n,x)[1][n]
		return yderiv

	pp = sp.special.sph_yn(n,x)[0][n]+x*yp(n,x)
	return -1.0*pp

# Bn from continuity of the first derivative
Bn = (m*Psiprime(n,m*k*a)\
	-sp.special.riccati_jn(n,(m*k*a))[0][n]/sp.special.riccati_jn(n,(m*k*a))[0][n]*Psiprime(n,k*a))\
	/(Xiprime(n,k*a)\
	-sp.special.riccati_yn(n,k*a)[0][n]/sp.special.riccati_jn(n,k*a)[0][n]*Psiprime(n,k*a))

print 'Bn = ', Bn

# beta_n from continuity of the radial wave function at the interface
beta_n = (sp.special.riccati_jn(n,m*k*a)[0][n]/Bn - sp.special.riccati_yn(n,k*a)[0][n]) \
		/(sp.special.riccati_jn(n,(k*a))[0][n])

print 'beta_n = ', beta_n

for ii in range(0,r.size):
	if r[ii] > a:
		radwavefunc[ii] = Bn*(sp.special.riccati_yn(n,(k*r[ii]))[0][n] \
						+ beta_n*sp.special.riccati_jn(n,(k*r[ii]))[0][n])
	elif r[ii] <= a:
		radwavefunc[ii] = sp.special.riccati_jn(n,(m*k*r[ii]))[0][n]
fig4 = plt.figure(4)
scalingfactor = 1000.
plt.plot(r,radwavefunc,'r-',linewidth=2)
fig4.show()

#=======================================================================================================
# Numerical ODE solution approach

fig2 = plt.figure(2)
#plt.plot(r,V/np.amax(V),'k-',linewidth=2)
axes = plt.gca()
#axes.set_xlim([0.5,2.0])
#axes.set_ylim([0.,2000.])
plt.text(a, 1800, 'scattering particle radius')
plt.text(0.41, k*k, '$k^2$',fontsize=15)
plt.xlabel('radial coordinate $r$',fontsize=20)
plt.ylabel('potential $V_n(r)$', fontsize=20)
plt.grid('on'), plt.box('on')

from scipy.integrate import odeint

# Define Right-Hand-Side (RHS) function of the ODE \
# x is the dependent variable. \
# y is the 2-element vector of the solution y[1] and its derivative y[0].
def RHS(y,x):
	# k and n same as above.
	k = 34.611195 # resonance value for numerical integration in accordance with ref. [1]
	n = 40.
	# define refractive index profile:
	def mr(r):
		if r < 1.0:
			return 1.47
		elif r > 1.:
			return 1.00
	return [(n*(n+1)/(x*x)-k*k*mr(x)*mr(x))*y[1],\
			y[0]]

y0 = [0.01,0.]					# initial condition \
								# ODE solver cannot be startet with [0.,0.] or at xr = 0.

xr = np.arange(0.01,3.0,0.01)	# dependent variable, stepsize was chosen arbitrarily.
Sn = odeint(RHS, y0, xr)		# solve ODE
#plt.plot(xr,Sn[:,0],linewidth=2)
plt.plot(xr,Sn[:,1]/np.amax(Sn[:,1]),linewidth=2)
plt.plot(r,radwavefunc/np.amax(radwavefunc)+0.5,linewidth=2)
fig2.show()

#=======================================================================================================
# Estimation of resonance values of the Mie parameter

xB = (n+0.5)/m  # potential well bottom
xT = n+0.5		# potential well top

miepara = np.linspace(xB,xT,20000)
betavalues = np.zeros(miepara.shape,dtype=float)

def beta_residuum(m,n,x):

	# First derivative of the Riccati-Bessel function
	# of the first kind based on the spherical Bessel
	# function.
	def Psiprime(n,x):
		def jp(n,x):
			jderiv = sp.special.sph_jn(n,x)[0][n-1] \
					-(n+1)/x*sp.special.sph_jn(n,x)[0][n]
			return jderiv
		pp = sp.special.sph_jn(n,x)[0][n]+x*jp(n,x)
		return pp

	# First derivative of the Riccati-Bessel function
	# of the second kind based on the corresponding
	# spherical Bessel function.
	def Xiprime(n,x):
		def yp(n,x):
			yderiv = sp.special.sph_yn(n,x)[0][n-1] \
					-(n+1)/x*sp.special.sph_yn(n,x)[1][n]
			return yderiv

		pp = sp.special.sph_yn(n,x)[0][n]+x*yp(n,x)
		return -1.0*pp

	# Bn from continuity of the first derivative
	Bn = (m*Psiprime(n,m*x)\
		-sp.special.riccati_jn(n,(m*x))[0][n]/sp.special.riccati_jn(n,(m*x))[0][n]*Psiprime(n,x))\
		/(Xiprime(n,x)\
		-sp.special.riccati_yn(n,x)[0][n]/sp.special.riccati_jn(n,x)[0][n]*Psiprime(n,x))

	# beta_n from continuity of the radial wave function at the interface
	beta_n = (sp.special.riccati_jn(n,m*x)[0][n]/Bn - sp.special.riccati_yn(n,x)[0][n])
	#		/(sp.special.riccati_jn(n,x)[0][n]) 
	#beta_n = ((sp.special.riccati_yn(n,x)[0][n]*Bn) / sp.special.riccati_jn(n,m*x)[0][n])
	return abs(beta_n)

resonancecount = 0

for ii in range(0,miepara.size):
	betavalues[ii] = beta_residuum(m,n,miepara[ii])
	epsilon = 0.005
	if betavalues[ii] < 0.0 + epsilon and betavalues[ii] > 0.0 - epsilon:
		print miepara[ii]
		resonancecount += 1

print 'Number of MDRs / shape resonances: ', resonancecount

fig3 = plt.figure(3)
plt.plot(miepara,betavalues)
axes = plt.gca()
axes.set_yscale('log')
#axes.set_xlim([0.5,2.0])
#axes.set_ylim([0.,1.])
plt.grid('on'), plt.box('on')
fig3.show()

raw_input() # wait for input to avoid immediate destruction of figures.