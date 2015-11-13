# MDRplot
	A short python script for MDR calculations.

 	Python script for calculations associated with
	Morphology dependent resonances.
	The theory is based on:

	[1] Johnson; JOSA A/Vol. 10, No. 2, (Feb. 1993)

	Author: Patrick Stegmann (PST)
	Date: 2015-11-02

	Revisions:
	----------

	- Radial potential plot (2015-11-02,PST)
	- Radial wave function evaluation & plot (2015-11-11,PST)
	- Numerical solution to check semi-analytic approach (2015-11-12,PST)
	- Print number and location of MDR resonances (2015-11-13,PST)

	Programmer's Comments:
	----------------------
	
	1. The Bessel functions from Scipy always return a
		tuple with shape [1][n]. Undocumented "feature".
