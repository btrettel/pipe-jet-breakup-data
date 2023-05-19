#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# pipe-jet-breakup-data is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

# Not recommended but this seems to be the easiest way.
# https://docs.python.org/3/faq/programming.html#what-are-the-best-practices-for-using-import-in-a-module
from jetbreakup3 import *
#import jetbreakup3 as jb

# # bhunia_impingement_1993, fig. 3-11, bottom
# d_0   = 5.8e-3 # m
# We_l0 = 1444
# x_s   = 22.4

# bhunia_impingement_1993, fig. 3-12, top right
d_0   = 2.7e-3 * ureg.meter
We_l0 = 6906
x_s   = 20.

rho_l = 1000. * ureg.kilogram/(ureg.meter**3.)
sigma = 72.e-3 * ureg.newton/ureg.meter

# TODO: Make function to compute velocity from Weber number
Ubar_0 = np.sqrt(We_l0 * sigma / (rho_l * d_0))
print(Ubar_0.to_base_units())
assert Ubar_0.check('[length]/[time]')

# #freq_of_first_zero = (np.pi / (x_s * d_0))**(2./3.) * Ubar_0**(5./3.) * (rho_l/sigma)**(1./3.)
# freq_of_first_zero = ((np.pi / (x_s * d_0))**2. * Ubar_0**5. * (rho_l/sigma))**(1./3.)
# print(freq_of_first_zero.to_base_units())
# assert freq_of_first_zero.check('1/[time]')

f_arr = np.array([])
G_arr = np.array([])
with open('../data/bhunia_impingement_1993/bhunia_impingement_1993_fig_3-12_top_right.csv') as csv_file:
   csv_reader = csv.DictReader(csv_file, delimiter=',')
   
   for row in csv_reader:
      #print(row['Ubar_0*kappa_x'], row['~G'])
      f_arr = np.append(f_arr, float(row['Ubar_0*kappa_x']))
      G_arr = np.append(G_arr, float(row['~G']))
f_arr = f_arr / ureg.second

plt.loglog(f_arr.to_base_units(), G_arr)
#plt.axvline(freq_of_first_zero.to_base_units())
plt.grid()
plt.savefig('3-12.png')

v_s0p_est = 0.05 * Ubar_0

kappa_x_arr = f_arr / Ubar_0
print('Lowest wavenumber:', min(kappa_x_arr).to_base_units())
print('Damping ratio for lowest wavenumber:', (v_s0p_est / 4.) * np.sqrt(rho_l / (min(kappa_x_arr) * sigma)))
print('Highest wavenumber:', max(kappa_x_arr).to_base_units())
print('Damping ratio for highest wavenumber:', (v_s0p_est / 4.) * np.sqrt(rho_l / (max(kappa_x_arr) * sigma)))
