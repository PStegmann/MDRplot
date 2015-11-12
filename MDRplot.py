# 
# 	Short script for calculations associated with
#	Morphology dependent resonances.
#	The theory is based on:
#
#	Johnson; JOSA A/Vol. 10, No. 2, (Feb. 1993)
#
#	Author: Patrick Stegmann (PST)
#	Date: 2015-11-02
#
#	Revisions:
#	----------
#
#	- Radial potential plot (2015-11-02,PST)
#	- Radial wave function evaluation & plot (2015-11-11,PST)
#
#	Programmer's Comments:
#	----------------------
#	
#	1. The Bessel functions from Scipy always return a
#		tuple with shape [1][n]. Undocumented "feature".
# 

# import modules
import numpy as np 			# NumPy
import scipy as sp         		# Scipy
from scipy import special
import matplotlib as mpl 		# Matplotlib (2D/3D)
import matplotlib.pyplot as plt 	# Matplotlib's pyplot
from pylab import * 			# Matplotlib's pylab

# define constants
global m,n,k,a
m = 1.47 						# refractive index
n = 40 							# order
k = 33. 						# wave number
a = 1.0 						# particle radius
r = np.linspace(0.5,2.0,3000) 	# radial coordinate

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

plt.plot(r,V,'k-',linewidth=2)
axes = plt.gca()
axes.set_xlim([0.5,2.0])
axes.set_ylim([0.,2000.])
plt.text(a, 1800, 'scattering particle radius')
plt.text(0.41, k*k, '$k^2$',fontsize=15)
plt.xlabel('radial coordinate $r$',fontsize=20)
plt.ylabel('potential $V_n(r)$', fontsize=20)
plt.grid('on'), plt.box('on')

# radial wave function TE (eq. 12)
radwavefunc = np.zeros(r.shape,dtype=float)

# First derivative of the Riccati-Bessel function
# of the first kind based on the spherical Bessel
# function.
def Psiprime(n,x):
	def jp(n,x):
		jderiv = sp.special.sph_jn(n,x)[1][n-1] \
				-(n+1)/x*sp.special.sph_jn(n,x)[1][n]
		return jderiv
	pp = sp.special.sph_jn(n,x)[1][n]+x*jp(n,x)
	return pp

# First derivative of the Riccati-Bessel function
# of the second kind based on the corresponding
# spherical Bessel function.
def Xiprime(n,x):
	def yp(n,x):
		yderiv = sp.special.sph_yn(n,x)[1][n-1] \
				-(n+1)/x*sp.special.sph_yn(n,x)[1][n]
		return yderiv

	pp = sp.special.sph_yn(n,x)[1][n]+x*yp(n,x)
	return -1.0*pp

# Bn from continuity of the first derivative
Bn = (m*Psiprime(n,m*k*a)\
	-sp.special.riccati_jn(n,(m*k*a))[1][n]/sp.special.riccati_jn(n,(m*k*a))[1][n]*Psiprime(n,k*a))\
	/(Xiprime(n,k*a)\
	-sp.special.riccati_yn(n,k*a)[1][n]/sp.special.riccati_jn(n,k*a)[1][n]*Psiprime(n,k*a))

print 'Bn = ', Bn

# beta_n from continuity of the radial wave function at the interface
beta_n = (sp.special.riccati_yn(n,m*k*a)[1][n]/Bn - sp.special.riccati_jn(n,k*a)[1][n]) \
		/(sp.special.riccati_yn(n,(k*a))[1][n])

print 'beta_n = ', beta_n

for ii in range(1,r.size):
	if r[ii] > a:
		radwavefunc[ii] = Bn*(sp.special.riccati_jn(n,(k*r[ii]))[1][n] \
						+ beta_n*sp.special.riccati_yn(n,(k*r[ii]))[1][n])
	elif r[ii] <= a:
		radwavefunc[ii] = sp.special.riccati_jn(n,(m*k*r[ii]))[1][n]

scalingfactor = 1000.
plt.plot(r,scalingfactor*radwavefunc+scalingfactor,'r-',linewidth=2)

plt.show()
