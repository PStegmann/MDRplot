# 
# 	Short script for calculations associated with
#	Morphology dependent resonances.
#
#	Author: Patrick Stegmann (PST)
#	Date: 2015-11-02
#
#	Revisions:
#	----------
#
#	- Radial potential plot (2015-11-02,PST)
#	- Radial wave function evaluation & plot (2015-11-02,PST)
#

# import modules
import numpy as np 				# NumPy
import scipy as sp 				# SciPy
import matplotlib as mpl 		# Matplotlib (2D/3D)
import matplotlib.pyplot as plt # Matplotlib's pyplot
from pylab import * 			# Matplotlib's pylab

# define constants
global m,n,k,a
m = 1.47 						# refractive index
n = 40 							# order
k = 33.0 						# wave number
a = 1.0 						# particle radius
r = np.linspace(0.5,1.5,1000) 	# radial coordinate

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
plt.text(a, 1800, 'scattering particle radius')
plt.text(0.41, k*k, '$k^2$',fontsize=15)
plt.xlabel('radial coordinate $r$',fontsize=20)
plt.ylabel('potential $V_n(r)$', fontsize=20)
plt.grid('on'), plt.box('on')
plt.show()