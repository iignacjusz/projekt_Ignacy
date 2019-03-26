#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 14:37:48 2019

@author: ignacy
"""


import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# Use reaction mechanism GRI-Mech 3.0
gas1 = ct.Solution('gri30.xml')
gas2 = ct.Solution('gri30.xml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified
# equivalence ratio

#First option
print('Set equivalence ratio for methane/air mixture')

equiv_ratio1 = input()  # lean combustion for example: "0.5"
equiv_ratio1 = float(equiv_ratio1)
print('Set Temperature') #300
T1=input()
T1=float(T1)
print('Now pressure')#101325
p1=input()
p1=float(p1)
gas1.TP = T1, p1
gas1.set_equivalence_ratio(equiv_ratio1, 'CH4:1.0', 'O2:1.0, N2:3.76')


#Second option
print('Now set equivalence ratio for methane/oxygen')

equiv_ratio2 = input()  # lean combustion for example: "0.5"
equiv_ratio2 = float(equiv_ratio2)
print('Set Temperature') #300
T2=input()
T2=float(T2)
print('Now pressure')#101325
p2=input()
p2=float(p2)
gas2.TP = T2, p2
gas2.set_equivalence_ratio(equiv_ratio2, 'CH4:1.0', 'O2:1.0')




inlet1 = ct.Reservoir(gas1) #gas1
inlet2 = ct.Reservoir(gas2) #gas2
# Create the combustor, and fill it initially with a mixture consisting of the
# equilibrium products of the inlet mixture. This state corresponds to the state
# the reactor would reach with infinite residence time, and thus provides a good
# initial condition from which to reach a steady-state solution on the reacting
# branch.

#gas1
gas1.equilibrate('HP')
combustor1 = ct.IdealGasReactor(gas1)
combustor1.volume = 1.0

#gas2
gas2.equilibrate('HP')
combustor2 = ct.IdealGasReactor(gas2)
combustor2.volume = 1.0

# Create a reservoir for the exhaust
exhaust1 = ct.Reservoir(gas1) #gas1
exhaust2 = ct.Reservoir(gas2) #gas2

# Use a variable mass flow rate to keep the residence time in the reactor
# constant (residence_time = mass / mass_flow_rate). The mass flow rate function
# can access variables defined in the calling scope, including state variables
# of the Reactor object (combustor) itself.

#gas1
def mdot(t):
    return combustor1.mass / residence_time1

inlet1_mfc = ct.MassFlowController(inlet1, combustor1, mdot=mdot)

#gas2
def mdot(t):
    return combustor2.mass / residence_time2

inlet2_mfc = ct.MassFlowController(inlet2, combustor2, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'master'
# MassFlowController, with an additional pressure-dependent term. By explicitly
# including the upstream mass flow rate, the pressure is kept constant without
# needing to use a large value for 'K', which can introduce undesired stiffness.
outlet1_mfc = ct.PressureController(combustor1, exhaust1, master=inlet1_mfc, K=0.01) #gas1
outlet2_mfc = ct.PressureController(combustor2, exhaust2, master=inlet2_mfc, K=0.01) #gas2

# the simulation
sim1 = ct.ReactorNet([combustor1]) #gas1
sim2 = ct.ReactorNet([combustor2]) #gas2

# Run a loop over decreasing residence times, until the reactor is extinguished,
# saving the state after each iteration.
states1 = ct.SolutionArray(gas1, extra=['tres1'])
states2 = ct.SolutionArray(gas2, extra=['tres2'])

residence_time2=0.1
residence_time1 = 0.1  # starting residence time
while combustor1.T > 500 or combustor2.T > 500:
    sim1.set_initial_time(0.0)  # reset the integrator
    sim2.set_initial_time(0.0)
    sim1.advance_to_steady_state()
    sim2.advance_to_steady_state()
    if combustor1.T > T1:
        print('tres1 = {:.2e}; T1 = {:.1f}'.format(residence_time1, combustor1.T))
        states1.append(combustor1.thermo.state, tres1=residence_time1)
        residence_time1 *= 0.9  # decrease the residence time for the next iteration
    if combustor2.T > T2:
        print('tres2 = {:.2e}; T2 = {:.1f}'.format(residence_time2, combustor2.T))
        states2.append(combustor2.thermo.state, tres2=residence_time2)
        residence_time2 *= 0.9 

Q1 = - np.sum(states1.net_production_rates * states1.partial_molar_enthalpies, axis=1)
Q2 = - np.sum(states2.net_production_rates * states2.partial_molar_enthalpies, axis=1)

f, ax1 = plt.subplots(1,1)
ax1.plot(states1.tres1, Q1, '.-', color='C0')
ax2 = ax1.twinx()
ax2.plot(states1.tres1[:-1], states1.T[:-1], '.-', color='C1')
ax1.set_xlabel('residence time1 [s]')
ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
ax2.set_ylabel('temperature [K]', color='C1')
plt.title('Methan/Air')
f.tight_layout()
plt.show()

f, ax1 = plt.subplots(1,1)
ax1.plot(states2.tres2, Q2, '.-', color='C0')
ax2 = ax1.twinx()
ax2.plot(states2.tres2[:-1], states2.T[:-1], '.-', color='C1')
ax1.set_xlabel('residence time2 [s]')
ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
ax2.set_ylabel('temperature [K]', color='C1')
plt.title('Methan/Oxygen')
f.tight_layout()
plt.show()