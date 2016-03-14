"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio for three different temperatures.

The results are intended to produce results in the same format as the results 
from Heywood's Fundamentals of Internal Combustion Engines, 1988, Fig 3.10

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

#Duplicate the conditions for Heywood
T = np.zeros(4)
T = [0,1750,2250,2750] #[K]
P = 30 * 101325 #Pa

#Note that the cti file below is a crude combination of two cti files.  This is 
# because the octane reactions needed are found in the Curran et al research, the 
# NO reactions are found in the GRI30 files, and other fuels can be found in the LLNL
# C1C4 reactions.

gas=ct.Solution('Curran_GRI30_NO_LLNLC1C4NO.cti')

# gaseous fuel species
fuel_species = 'IC8H18'

# air composition
air_N2_O2_ratio = 3.76
A = 12.5

# equivalence ratio range
phi_min = 0.2
phi_max = 1.4
npoints = 100

##############################################################################

# create some arrays to hold the data
phi = np.zeros(npoints)
tad = np.zeros(npoints)
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

print('Calculating the equilibrium gas mixture resulting from the combustion of C8H18 (octane)')
print(' at one pressure and three different temperatures.')
print('Note that this calculation requires several minutes and the result will pop up in a ')
print(' seperate window which might not show up on top')

fig = plt.figure(facecolor="white")

for j in range(1,4):
#    print(j)
    temp_num = T[j]
    temp_str = str(temp_num)
    graph = fig.add_subplot(1,3,j)

    for i in range(npoints):
        phi[i] = phi_min + (phi_max - phi_min)*i/(npoints - 1)
        X = np.zeros(gas.n_species)
        X[ifuel] = phi[i]
        X[io2] = stoich_O2
        X[in2] = stoich_O2*air_N2_O2_ratio

        # set the gas state
        gas.TPX = T[j], P, X

        # equilibrate the mixture at constant temperature and pressure
        gas.equilibrate('TP')

        tad[i] = gas.T
        xeq[i,:] = gas.X
        for k in range(gas.n_species):
            if xeq[i,k]>xeq_max[k]:
                xeq_max[k]= xeq[i,k]
                phi_xeq_max[k]=phi[i]

    xeq_max[gas.n_species]=0
    for k in range(gas.n_species):
        index_big2small[k]=gas.n_species
        for m in range(gas.n_species):
            if (xeq_max[m] >= xeq_max[index_big2small[k]]):
                present = False
                for n in range(k):
                    if index_big2small[n]==m:
                        present = True
                if present == False:
                    index_big2small[k] = m

    plt.plot(phi,tad)
    plt.semilogy(phi,xeq)
    plt.axis([phi[0],phi[npoints-1],1.0e-4,1])
    plt.minorticks_on()
    plt.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
    plt.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Mole Fraction')
    plt.title('Pressure = ' + str(P/101325) + 'BarA' + '\n' + 'Temperature = '+ str(T[j]) +'K')

    for k in range (16):
        label_offset_x = 20
        label_offset_y = 5
        if phi_xeq_max[index_big2small[k]]>(phi_min+phi_max)/2:
            label_offset_x = label_offset_x * -2
        if gas.species_name(index_big2small[k]) == 'CO2':
            label_offset_y = -8* label_offset_y
        graph.annotate(gas.species_name(index_big2small[k]), xy=(phi_xeq_max[index_big2small[k]], xeq_max[index_big2small[k]]),  xycoords='data',
        xytext=(label_offset_x, label_offset_y), textcoords='offset points',
            size=16,
                )

print('Calculations complete, please check for graphical output window')
plt.show()


