"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio for three different temperatures.

The results are intended to produce results in the same format as the results 
from Heywood's Fundamentals of Internal Combustion Engines, 1988, Fig 3.11

The starting point of this file was from:

http://cantera.org/docs/sphinx/html/cython/examples/reactors_combustor.html
"""

#  Alan Michel Thomason
#  alanmthomason@icloud.com
#  2016-02-06

import sys
import cantera.interrupts
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import csv

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the gasture.


T = 700 #K
P = 10 * 101E3 #Pa

#Note that the cti file below is a crude combination of two cti files.  This is 
# because the octane reactions needed are found in the Curran et al research, the 
# NO reactions are found in the GRI30 files, and other fuels can be found in the LLNL
# C1C4 reactions.

gas=ct.Solution('Curran_GRI30_NO_LLNLC1C4NO.cti')


print('Calculating the equilibrium gas mixture temperature and pressure resulting from the combustion of C8H18 (octane)')
print('Note that this calculation requires several minutes and the result will pop up in a ')
print(' seperate window which might not show up on top')

# gaseous fuel species
fuel_species = 'IC8H18'

# air composition
air_N2_O2_ratio = 3.76
A = 12.5

# equivalence ratio range
phi_min = 0
phi_max = 1.4
npoints = 20

##############################################################################

# create some arrays to hold the data
phi = np.zeros(npoints)
tad_v = np.zeros(npoints)
tad_p = np.zeros(npoints)
pad = np.zeros(npoints)
xeq = np.zeros((npoints,gas.n_species))
xeq_max = np.zeros(gas.n_species+1)
phi_xeq_max = np.zeros(gas.n_species)
index_big2small = np.zeros(gas.n_species+1)

# find fuel, nitrogen, and oxygen indices
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

if gas.n_atoms(fuel_species,'O') > 0 or gas.n_atoms(fuel_species,'N') > 0:
    raise "ERROR: only hydrocarbon fuels are supported."

stoich_O2 = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')


for i in range(npoints):
    phi[i] = phi_min + (phi_max - phi_min)*i/(npoints - 1)
    X = np.zeros(gas.n_species)
    X[ifuel] = phi[i]
    X[io2] = stoich_O2
    X[in2] = stoich_O2*air_N2_O2_ratio

    # set the gas state
    gas.TPX = T, P, X

    # equilibrate the mixture adiabatically at constant P
    gas.equilibrate('HP')
    tad_p[i]=gas.T
	
	# set the gas state
    gas.TPX = T, P, X
	
	# equilibrate the mixture adiabatically at constant P
    gas.equilibrate('UV')

    tad_v[i] = gas.T
    pad[i] = gas.P / 101325


fig = plt.figure(facecolor="white")
ax1 = fig.add_subplot(1,1,1)

ax1.plot(phi,tad_v, label = 'T_v', color = 'r')
ax1.plot(phi,tad_p, label = 'T_p', color = 'g')
ax1.set_xlabel('Equivalence Ratio')
ax1.set_ylabel('Adiabatic Flame Temperature [K]')
ax1.set_yticks(np.arange(0,3600,400))
ax1.minorticks_on()
ax1.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
ax1.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)
plt.legend(loc = 'upper left')

ax2 = ax1.twinx() 
ax2.plot(phi, pad, label = 'P_v', linestyle = '--', color = 'b')  
ax2.set_ylabel('Pressure, atm') 
ax2.axis([phi[0],phi[npoints-1],0,80])
ax2.minorticks_on()

plt.legend(loc = 'lower right')

plt.show()


