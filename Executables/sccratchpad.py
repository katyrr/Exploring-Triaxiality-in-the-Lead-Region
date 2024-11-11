#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 17:45:40 2024

@author: katyrr
"""




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:07:25 2024

@author: katyrr


run from Outputs directory in terminal, with command: 
    python3 ../../../Executables/write_inputs.py

or run from console with: 
    runfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/Tl207/MO/Outputs')

to debug (use breakpoints) run from console with: 
    debugfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/tl207/MO/Outputs')

"""




import numpy as np
import matplotlib.pyplot as plt


#eps_to_test = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
#eps_to_test = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
#eps_to_test = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
eps_min = 0.0
eps_max = 0.83
eps_to_test = np.linspace(eps_min, eps_max, num=6)

eps_to_plot = [] 
gamma_to_plot = []
fermi_energies_to_plot = [] 

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_thetamin(0)   # Start angle in degrees
ax.set_thetamax(60)  # End angle in degrees

theta_ticks = np.arange(0, 70, 10)  # Create tick locations in degrees
ax.set_xticks(np.radians(theta_ticks))  # Convert to radians for set_xticks

for l in range(len(eps_to_test)):
    for p in range(l+1):
        eps = eps_to_test[l]
        eps_to_plot.append(eps)
        if l==0:
            gamma = 0
        else:
            gamma = p*(60/l)*np.pi/180
        
        gamma_to_plot.append(gamma)
        fermi_energies_to_plot.append(eps)
        plt.polar(gamma, eps, 'kx')

        

cax = ax.tricontourf(gamma_to_plot, eps_to_plot, fermi_energies_to_plot, levels=10)

plt.colorbar(cax)
plt.show()
    



'''
string_list = ["1.1", "2.2", "3.3"]
float_list = [float(num) for num in string_list]

gamma_to_test = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
eps_to_test = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]
fermi_energies =  [5.9771, 6.1252, 6.1799, 6.1011, 6.1113, 6.0125, 5.9858, 6.1180, 6.1791, 6.1069, 6.1088, 6.0142, 6.0315, 6.1011, 6.1770, 6.1168, 6.1045, 6.0195, 6.0148, 6.0896, 6.1734, 6.1272, 6.1022, 6.0294, 6.0257, 6.1312, 6.1685, 6.1372, 6.1043, 6.0524, 6.0658, 6.1375, 6.1623, 6.1465, 6.1114, 6.0770, 6.1026, 6.1231, 6.1549, 6.1549, 6.1231, 6.1026]


gamma_to_test = np.array(gamma_to_test)
eps_to_test = np.array(eps_to_test)
fermi_energies = np.array(fermi_energies) 

gamma_to_plot, eps_to_plot = np.meshgrid(gamma_to_test, eps_to_test)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_thetamin(0)   # Start angle in degrees
ax.set_thetamax(60)  # End angle in degrees

theta_ticks = np.arange(0, 70, 10)  # Create tick locations in degrees
ax.set_xticks(np.radians(theta_ticks))  # Convert to radians for set_xticks

# (eps,gamma) and (-eps,60-gamma) correspnod to the same shape - covert to positive eps so that it can be plotted in polar coordinates
for r in range(len(gamma_to_test)):
    for c in range(len(eps_to_test)):
        if(eps_to_plot[c][r]<0):
            eps_to_plot[c][r] *= -1
            gamma_to_plot[c][r] = 60-gamma_to_plot[c][r]
        gamma_to_plot[c][r] *= np.pi/180
        plt.polar(gamma_to_plot[c][r], eps_to_plot[c][r], 'wx')
            
fermi_energies_to_plot = np.array(fermi_energies).reshape((len(eps_to_test),len(gamma_to_test)))

cax = ax.contourf(gamma_to_plot, eps_to_plot, fermi_energies_to_plot, 10)

plt.colorbar(cax)
plt.show()
'''
