#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

from jetbreakup import *

#revno = lastchangedrevision[0:6]
revno = None

# MAYBE: Make a function to plot points in the smooth pipe regime diagram. Use this to plot all photos in a regime diagram?

macros_regime = open('../outputs/macros/regime.tex', 'w')

input_filename = 'pipe-jet-breakup-data'

print

# TODO: Add lightly colored boxes for fire hoses, diesel injectors, jet sprinklers, etc. to these plots.

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

i = 0
combined_regime_array = []
for regime_photo, regime_turb in zip(df_jet_breakup['regime photo'], df_jet_breakup['regime turb']):
   if isinstance(regime_photo, basestring):
      #print i, regime_photo
      if regime_photo == 'bursting':
         regime_photo = 'downstream transition'
      
      if (regime_photo == 'Rayleigh') and (regime_turb == 'laminar'):
         regime_photo = 'laminar Rayleigh'
      
      if (regime_photo == 'Rayleigh') and (regime_turb == 'transitional'):
         regime_photo = 'transitional Rayleigh'
      
      if (regime_photo == 'Rayleigh') and (regime_turb == 'turbulent'):
         regime_photo = 'turbulent Rayleigh'
      
      if regime_photo == 'second wind-induced':
         regime_photo = 'turbulent surface breakup'
      
      if regime_photo == 'first wind-induced':
         regime_photo = 'downstream transition'
      
      combined_regime_array.append(regime_photo)
   else:
      combined_regime_array.append(np.nan)
   
   i = i + 1

i = 0
for regime_L_b, regime_turb in zip(df_jet_breakup['regime L_b'], df_jet_breakup['regime turb']):
   if isinstance(regime_L_b, basestring):
      if (regime_L_b == 'Rayleigh') and (regime_turb == 'laminar'):
         regime_L_b = 'laminar Rayleigh'
      
      if (regime_L_b == 'Rayleigh') and (regime_turb == 'transitional'):
         regime_L_b = 'transitional Rayleigh'
         
      if (regime_L_b == 'Rayleigh') and (regime_turb == 'turbulent'):
         regime_L_b = 'turbulent Rayleigh'
      
      if regime_L_b == 'second wind-induced':
         regime_L_b = 'turbulent surface breakup'
      
      if regime_L_b == 'first wind-induced':
         regime_L_b = 'downstream transition'
      
      if isinstance(combined_regime_array[i], basestring):
         if not(combined_regime_array[i] == regime_L_b):
            raise ValueError('Inconsistent regimes:', i, regime_L_b, combined_regime_array[i])
      else:
         combined_regime_array[i] = regime_L_b
   
   i = i + 1

We_l0_p = []
for We_l0, Tu_0 in zip(df_jet_breakup['We_l0'], df_jet_breakup['I_0']):
   We_l0_p.append(We_l0 * Tu_0**2.)

marker_array = ['o', 's', 'D', 'x', '+', '>', '<', '^', 'v', '.']
#color_array = ['blue', 'green', 'red', 'orange', 'y']

regime_all_df = df_jet_breakup
regime_all_df['regime combined'] = combined_regime_array
regime_all_df = regime_all_df[regime_all_df['regime combined'].notnull()]

assert(not('bursting' in regime_all_df['regime combined']))
summary_table(regime_all_df)

macros_regime.write(r'\newcommand{\totalregimenum}{\num{'+str(len(regime_all_df))+'}}\n')
macros_regime.write(r'\newcommand{\photoregimenum}{\num{'+str(len(regime_all_df[regime_all_df['regime photo'].notnull()]))+'}}\n')
macros_regime.write(r'\newcommand{\xbavgregimenum}{\num{'+str(len(regime_all_df[regime_all_df['regime L_b'].notnull()]))+'}}\n')

# https://tex.stackexchange.com/a/273112/9945
macros_regime.write(r'\newcommand{\Relorange}{\numrange[round-mode=places,round-precision=1,retain-zero-exponent=true]{'+'%.1E' % Decimal(np.amin(regime_all_df['Re_l0']))+'}{'+'%.1E' % Decimal(np.amax(regime_all_df['Re_l0']))+'}}\n')
macros_regime.write(r'\newcommand{\Welorange}{\numrange[round-mode=places,round-precision=1,retain-zero-exponent=true]{'+'%.1E' % Decimal(np.amin(regime_all_df['We_l0']))+'}{'+'%.1E' % Decimal(np.amax(regime_all_df['We_l0']))+'}}\n')

df_key_turb = regime_all_df[regime_all_df['regime turb'] == 'turbulent']
macros_regime.write(r'\newcommand{\Tuorange}{\numrange{'+str(round(np.amin(df_key_turb['I_0']), 3))+'}{'+str(round(np.amax(df_key_turb['I_0']), 3))+'}}\n')

regime_all_df_with_eisenklam = regime_all_df
regime_all_df = regime_all_df[regime_all_df['key'] != 'eisenklam_flow_1958'] # Removed as it seems to have atypically high transition Reynolds numbers.
regime_all_df = regime_all_df[regime_all_df['key'] != 'sterling_instability_1975'] # Removed as it seems to have atypically high transition Reynolds numbers.

#summary_table(regime_all_df[regime_all_df['key'] == 'grant_newtonian_1965'])
#summary_table(regime_all_df[regime_all_df['key'] == 'kusui_liquid_1969'])

regime_all_df = regime_all_df[regime_all_df['rho_s'] > 500]
#regime_all_df = regime_all_df[regime_all_df['photo filename'].notnull()]
regime_all_df = regime_all_df[regime_all_df['roughness'] == 'smooth']

# overall regime plot for smooth pipes

regime_all_df_with_trans = regime_all_df
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'F2S']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'R2F']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'R2B']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'S2A']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'R2S']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'F2S']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'F2R']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'transitional Rayleigh']
regime_all_df = regime_all_df[regime_all_df['regime combined'] != 'bursting']

regime_array = []
for regime in regime_all_df['regime combined']:
   if not(regime in regime_array):
      regime_array.append(regime)

print len(regime_all_df), 'data points with breakup regime classified (smooth pipe nozzles).'
macros_regime.write(r'\newcommand{\smoothregimenum}{\num{'+str(len(regime_all_df))+'}}\n')

#R2S_df = regime_all_df_with_trans[regime_all_df_with_trans['regime combined'] == 'R2S']
#avgTu = np.average(R2S_df['I_0'])
avgTu = 0.05
Re_Rto2WI = np.logspace(np.log(Re_turb) / np.log(10.), 6., 1e2)
#Tu_vec = I_fully_developed_array(Re_Rto2WI)
Tu_vec = avgTu * np.ones(100)
We_Rto2WI = We_l0_crit_TR(Tu_vec)
plt.loglog(We_Rto2WI, Re_Rto2WI, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--') #, label=r'$\mathrm{We}_\text{T,crit}$')

Re_2WItoA = Re_Rto2WI
##We_2WItoA = 0.12 * Tu_vec**(-1.26) * (1000./1.2)**1.53
#We_2WItoA = 0.12 * Tu_vec**(-5./4.) * (1000./1.2)**(3./2.)
We_2WItoA = We_l0_crit(Tu_vec, 1000./1.2)
plt.loglog(We_2WItoA, Re_2WItoA, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
#We_Kolmogorov = 0.01022 * Tu_vec**(-3./4.) * Re_2WItoA**(5./4.)
#plt.loglog(We_Kolmogorov, Re_2WItoA, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

# straight line
#We_corner_RtoF = 1.e2
# We_RtoF = np.logspace(np.log(We_corner_RtoF) / np.log(10.), 6., 1e2)
# Re_RtoF = Re_trans * (We_RtoF / We_corner_RtoF)**(-0.5)
#print 'Re_x_trans =', 13.4 * Re_trans * We_corner_RtoF**0.5 + 3. * We_corner_RtoF

# Sterling
#def We_crit_Sterling(We_l0, *data):
   #Re_l0 = data
   #rho_g = 1.2
   #rho_l = 1000.
   #return We_l0 * (rho_g / rho_l) - (1.2 + 3.41 * (We_l0**(0.5) / Re_l0)**(0.9))
#We_corner_RtoF = fsolve(We_crit_Sterling, 1.e3, args=(Re_trans))[0]
#We_RtoF = np.logspace(np.log(We_corner_RtoF) / np.log(10.), 6., 1e2)
#Re_RtoF = (3.91 * We_RtoF**0.54) / (We_RtoF * (1.2/1000.) - 1.2)**(1.11)

#Re_RtoF = Re_trans * (We_RtoF / We_corner_RtoF)**(-0.2) # similar to Grant
#We_RtoF = np.logspace(3., 6., 1e2)
#Re_RtoF = 125 * We_RtoF**(-0.19)

# new Feb. 2020 theory
Re_x_trans  = 1.7e5
C_LR        = 8.5
We_low      = ((sqrt(Re_trans**2. - 12. * Re_x_trans / C_LR) - Re_trans) / 6.)**2.
We_high     = Re_x_trans / (3. * C_LR)
We_RtoF_all = np.logspace(np.log(We_low) / np.log(10.), np.log(We_high) / np.log(10.), 1e2)
Re_RtoF_all = Re_l0_crit_DT(We_RtoF_all, Re_x_trans=Re_x_trans, C_LR=C_LR)
#plt.loglog(We_RtoF, Re_RtoF, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--') #, label=r'$\mathrm{We}_\text{T,crit}$')

Re_stable = Re_x_trans / 1.e3
We_stop = np.interp(Re_stable, Re_RtoF_all[::-1], We_RtoF_all[::-1])
We_RtoF_obs = np.logspace(np.log(We_low) / np.log(10.), np.log(We_stop) / np.log(10.), 1e2)
We_RtoF_far = np.logspace(np.log(We_stop) / np.log(10.), np.log(We_high) / np.log(10.), 1e2)
Re_RtoF_obs = Re_l0_crit_DT(We_RtoF_obs, Re_x_trans=Re_x_trans, C_LR=C_LR)
Re_RtoF_far = Re_l0_crit_DT(We_RtoF_far, Re_x_trans=Re_x_trans, C_LR=C_LR)

plt.loglog(We_RtoF_obs, Re_RtoF_obs, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
plt.loglog(We_RtoF_far, Re_RtoF_far, marker=None, color='k', zorder=4, linewidth=0.5, linestyle=':')
plt.loglog([We_stop, 1.e6], [Re_stable, Re_stable], marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

plt.axhline(Re_trans, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
#plt.axvline(4., marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
#dripping (old):
#plt.loglog(np.array([4., 4.]), np.array([1.e0, Re_trans]), marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
#plt.loglog(np.array([4., 4.]), np.array([Re_turb, 1.e6]), marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
plt.axhline(Re_turb, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

def dripping_We_func(We_lo):
   d_0        = rho_l * (Re_lo * nu_l)**2. / (sigma * We_lo)
   Bo_lo      = np.sqrt(rho_l * g * d_0**2. / (2. * sigma))
   Bo_louter  = np.sqrt(rho_l * g * (d_0 + t)**2. / (2. * sigma))
   We_lo_crit = 4. * (Bo_louter / Bo_lo) * (1. + K * Bo_lo * Bo_louter - ((1. + K * Bo_lo * Bo_louter)**2. - 1)**(1./2.))**2.
   return We_lo_crit - We_lo

rho_l       = rho_water(20.)
nu_l        = nu_water(20.)
sigma       = sigma_water(20.)
K           = 0.37  # clanet_transition_1999
t           = 0.    # tube thickness
Re_lo_arr   = np.logspace(0., np.log(500.) / np.log(10.), 1e3)
We_lo_arr   = np.array([])
We_lo_guess = 4.
for Re_lo in Re_lo_arr:
   We_lo = fsolve(dripping_We_func, We_lo_guess)[0]
   We_lo_arr = np.append(We_lo_arr, We_lo)
   We_lo_guess = We_lo

# This is probably slightly inaccurate but smooths out the part that probably has numerical error.
Re_lo_arr = np.append(Re_lo_arr, max(Re_lo_arr))
We_lo_arr = np.append(We_lo_arr, 1.)

plt.loglog(We_lo_arr, Re_lo_arr, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

lam_Re_center  = np.exp(0.5 * (np.log(Re_trans) + np.log(1.e0)))
turb_Re_center = np.exp(0.5 * (np.log(Re_turb) + np.log(1.e6)))

plt.text(2.e0, 3.e1, 'dripping', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
plt.text(np.exp(0.5 * (np.log(We_high) + np.log(4.))), lam_Re_center, 'laminar Rayleigh', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(np.exp(0.5 * (np.log(1.e0) + np.log(We_Rto2WI[0]))), turb_Re_center, 'turbulent Rayleigh', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(5.e4, np.exp(0.5 * (np.log(Re_trans) + np.log(Re_stable))), 'downstream transition', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(np.exp(0.5 * (np.log(We_high) + np.log(1.e6))), np.exp(0.5 * (np.log(3.e0) + np.log(Re_stable))), r'$\frac{x_\text{tr}}{d_0} > 10^3$', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(np.exp(0.5 * (np.log(We_2WItoA[0]) + np.log(1.e6))), turb_Re_center, 'atomization', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(1.e3, np.exp(0.5*(np.log(Re_trans) + np.log(Re_turb))), 'transitional at nozzle exit (critical $\mathrm{Re}_{\ell0}$ varies)', backgroundcolor='w', fontsize='small', verticalalignment='center', horizontalalignment='center', bbox=dict(boxstyle='square,pad=0.0',fc='white', ec='none'))
plt.text(np.exp(0.5 * (np.log(We_Rto2WI[0]) + np.log(We_2WItoA[0]))), turb_Re_center, 'turbulent\nsurf.\ breakup', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
plt.text(1.3e6, turb_Re_center, 'turb.\ at nozzle exit', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='left')
plt.text(1.3e6, lam_Re_center, 'laminar at nozzle exit', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='left')
plt.text(1.e3, 1.5e0, r'\textbf{do not use to determine regime}', backgroundcolor='w', verticalalignment='bottom', horizontalalignment='center', color='#777777')
plt.xlabel(r'$\mathrm{We}_{\ell0}$')
plt.ylabel(r'$\mathrm{Re}_{\ell0}$')
plt.grid()
#plt.tight_layout() # TODO: comment out for talk
axes = plt.gca()
if not(revno is None):
   box1 = TextArea('Rev. '+revno, textprops=dict(color="k"))
   anchored_box = AnchoredOffsetbox(loc=4,
                                    child=box1,
                                    frameon=False,
                                    bbox_to_anchor=(1., -0.20),
                                    bbox_transform=axes.transAxes)
   axes.add_artist(anchored_box)
fig = plt.gcf()
plt.xlim([1.e0, 1.e6])
plt.ylim([1.e0, 1.e6])
fig.set_size_inches(6., 4., forward=True) # report
plt.savefig('../outputs/figures/regime_map_low_atm_density.png')
plt.savefig('../outputs/figures/regime_map_low_atm_density.pgf', bbox_inches="tight")
fig.set_size_inches(5.5, 4., forward=True) # paper
plt.savefig('../outputs/figures/regime_map_low_atm_density_paper.pgf', bbox_inches="tight")
plt.close()

i = 0
for regime in regime_array:
   if regime == 'R2F':
      regime_print = 'Rayleigh to first wind-induced transition'
   elif regime == 'F2S':
      regime_print = 'first wind-induced to second wind-induced transition'
   elif regime == 'S2A':
      regime_print = 'second wind-induced to atomization transition'
   elif regime == 'R2B':
      regime_print = 'Rayleigh to bursting transition'
   else:
      regime_print = regime
   this_regime_all_df = regime_all_df[regime_all_df['regime combined'] == regime]
   plt.loglog(this_regime_all_df['We_l0'], this_regime_all_df['Re_l0'], linestyle='None', marker=marker_array[i], label=regime_print)
   
   i = i + 1

#d_0 = 0.25 * 2.54e-2 # m
#d_0 = 2.519e-3 # Clanet dripping tests
d_0 = 6.e-3 # m
rho_l = rho_water(T_std)
nu_l  = nu_water(T_std)
sigma = sigma_water(T_std)
Re_0_arr = np.logspace(0., 6., 1e2)
We_0_arr = (rho_l * nu_l**2.) / (sigma * d_0) * Re_0_arr**2.
plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None, label=r'$d_0 =$ 6~mm, water', zorder=1)

# "conventional" case
# # n-dodecane
# d_0   = 3.e-3 # m ==> Won't work any longer due to \sigma being wrong before
# rho_l = 749.5 # kg/m^3, https://en.wikipedia.org/wiki/Dodecane
# nu_l  = 1.34e-3 / rho_l # m^2/s, https://en.wikipedia.org/wiki/Dodecane
# sigma = 25.35e-3 # N/m, http://www.surface-tension.de/

# # water
# d_0 = 0.15e-3 # m
# # d_0 = 2.519e-3 # Clanet dripping tests
# rho_l = rho_water(T_std)
# nu_l  = nu_water(T_std)
# sigma = sigma_water(T_std)

# pure isopropyl alcohol
d_0   = 1.e-3 # m
rho_l = 0.786 * rho_water(T_std) # kg/m^3
nu_l  = 2.036e-3 / rho_l # m^2/s
sigma = 21.74e-3 # N/m

Re_0_arr = np.logspace(0., 6., 1e2)
We_0_arr = (rho_l * nu_l**2.) / (sigma * d_0) * Re_0_arr**2.
#plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None, label=r"$d_0 =$ 3~mm, $n$-dodecane (``conventional'' progression)", zorder=1)
plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None, label=r"$d_0 =$ 1~mm, isopropyl alcohol (``conventional'' progression)", zorder=1)

# Spray A
#d_0 = 90.e-6 # m, https://ecn.sandia.gov/diesel-spray-combustion/target-condition/spray-ab/
d_0 = 50.e-6 # m
rho_l = 749.5 # kg/m^3, https://en.wikipedia.org/wiki/Dodecane
nu_l  = 1.34e-3 / rho_l # m^2/s, https://en.wikipedia.org/wiki/Dodecane
sigma = 25.35e-3 # N/m, http://www.surface-tension.de/
Re_0_arr = np.logspace(0., 6., 1e2)
We_0_arr = (rho_l * nu_l**2.) / (sigma * d_0) * Re_0_arr**2.
plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None, label=r'$d_0 =$ 50~$\mu$m, $n$-dodecane', zorder=1)
#plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None, label=r'Spray A, diesel', zorder=1)

## Like Reitz (1978)
#d_0 = 0.34e-3 # m
#Re_0_arr = np.logspace(0., 6., 1e2)
#We_0_arr = (rho_l * nu_l**2.) / (sigma * d_0) * Re_0_arr**2.
#plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None, label=r'$d_0 =$ 0.34 mm, water', zorder=1)

#We_T_crit = 8.
Tu_vec = I_fully_developed_array(Re_Rto2WI)
We_Rto2WI = We_l0_crit_TR(Tu_vec)
plt.loglog(We_Rto2WI, Re_Rto2WI, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--') #, label=r'$\mathrm{We}_\text{T,crit}$')

#Re_2WItoA = Re_Rto2WI
##We_2WItoA = 0.12 * Tu_vec**(-1.26) * (1000./1.2)**1.53
#We_2WItoA = 0.12 * Tu_vec**(-5./4.) * (1000./1.2)**(3./2.)
We_2WItoA = We_l0_crit(Tu_vec, 1000./1.2)
plt.loglog(We_2WItoA, Re_2WItoA, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

# reusing arrays from in the map without the data
plt.loglog(We_RtoF_obs, Re_RtoF_obs, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
plt.loglog(We_RtoF_far, Re_RtoF_far, marker=None, color='k', zorder=4, linewidth=0.5, linestyle=':')
plt.loglog([We_stop, 1.e6], [Re_stable, Re_stable], marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

plt.axhline(Re_trans, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
#plt.loglog(np.array([4., 4.]), np.array([1.e0, Re_trans]), marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
#plt.loglog(np.array([4., 4.]), np.array([Re_turb, 1.e6]), marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
plt.axhline(Re_turb, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

# rho_l = rho_water(T_std)
# nu_l  = nu_water(T_std)
# sigma = sigma_water(T_std)
rho_l = 749.5 # kg/m^3, https://en.wikipedia.org/wiki/Dodecane
nu_l  = 1.34e-3 / rho_l # m^2/s, https://en.wikipedia.org/wiki/Dodecane
sigma = 25.35e-3 # N/m, http://www.surface-tension.de/
Re_lo_arr   = np.logspace(0., np.log(160.) / np.log(10.), 1e3)
We_lo_arr   = np.array([])
We_lo_guess = 4.
for Re_lo in Re_lo_arr:
   We_lo = fsolve(dripping_We_func, We_lo_guess)[0]
   We_lo_arr = np.append(We_lo_arr, We_lo)
   We_lo_guess = We_lo

# This is probably slightly inaccurate but smooths out the part that probably has numerical error.
Re_lo_arr = np.append(Re_lo_arr, max(Re_lo_arr))
We_lo_arr = np.append(We_lo_arr, 1.)

plt.loglog(We_lo_arr, Re_lo_arr, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')

plt.text(1.3e6, turb_Re_center, 'turb.\ at nozzle exit', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='left')
plt.text(1.3e6, lam_Re_center, 'laminar at nozzle exit', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='left')

plt.text(1.e3, 1.5e0, r'\textbf{do not use to determine regime}', backgroundcolor='w', verticalalignment='bottom', horizontalalignment='center', color='#777777')

#plt.legend(loc='best') #, fancybox=True, framealpha=0.5, numpoints=1)
plt.xlabel(r'$\mathrm{We}_{\ell0}$')
plt.ylabel(r'$\mathrm{Re}_{\ell0}$')
plt.grid()
#plt.tight_layout() # TODO: comment out for talk
axes = plt.gca()
if not(revno is None):
   box1 = TextArea('Rev. '+revno, textprops=dict(color="k"))
   anchored_box = AnchoredOffsetbox(loc=4,
                                    child=box1,
                                    frameon=False,
                                    bbox_to_anchor=(1., -0.20),
                                    bbox_transform=axes.transAxes)
   axes.add_artist(anchored_box)
fig = plt.gcf()
plt.xlim([1.e0, 1.e6])
plt.ylim([1.e0, 1.e6])

fig.set_size_inches(5.5, 4., forward=True) # talk
plt.savefig('../outputs/figures/regime_map_low_atm_density_with_data_and_lines_talk.pgf', bbox_inches="tight")

fig.set_size_inches(6., 4., forward=True) # report

# Move legend out of the plot to have extra space.
# https://stackoverflow.com/a/43439132
plt.legend(bbox_to_anchor=(0., -0.15), loc='upper left', frameon=False, fontsize=11)

plt.savefig('../outputs/figures/regime_map_low_atm_density_with_data_and_lines.png')
plt.savefig('../outputs/figures/regime_map_low_atm_density_with_data_and_lines.pgf', bbox_inches="tight")
fig.set_size_inches(5.5, 4., forward=True) # paper
plt.savefig('../outputs/figures/regime_map_low_atm_density_with_data_and_lines_paper.pgf', bbox_inches="tight")
plt.close()

d_0   = 6.e-3 # m
rho_l = rho_water(T_std)
nu_l  = nu_water(T_std)
sigma = sigma_water(T_std)

stability_curve(d_0, rho_l, nu_l, sigma, 'large')

#d_0 = (1./32.) * 2.54e-2 # m
#d_0 = 0.05e-3 # m
#stability_curve(d_0, rho_l, nu_l, sigma, 'conventional')

# pure isopropyl alcohol
d_0   = 1.e-3 # m
rho_l = 0.786 * rho_water(T_std) # kg/m^3
nu_l  = 2.036e-3 / rho_l # m^2/s
sigma = 21.74e-3 # N/m
stability_curve(d_0, rho_l, nu_l, sigma, 'conventional')

d_0 = 50.e-6 # m
rho_l = 749.5 # kg/m^3, https://en.wikipedia.org/wiki/Dodecane
nu_l  = 1.34e-3 / rho_l # m^2/s, https://en.wikipedia.org/wiki/Dodecane
sigma = 25.35e-3 # N/m, http://www.surface-tension.de/
stability_curve(d_0, rho_l, nu_l, sigma, 'engine')

# Older regime diagram in new coordinates
# WON'T: Put Magnotti's review version of the earlier regime diagram on the same We_l0 vs. Re_l0 coordinates you do as a comparison.
# reitz_atomization_1978 pdf p. 181 (table), pdf p. 189 (breakup length chart), pdf p. 191 (Ohnesorge diagram), pdf p. 193 (density ratio effect chart)

rho_l = rho_water(T_std)
nu_l  = nu_water(T_std)
sigma = sigma_water(T_std)
rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
rho_s = rho_l / rho_g

plt.axvline(8., marker=None, color='k', linewidth=0.8, linestyle='--') # dripping to Rayleigh

#We_Rto1WI = np.logspace(np.log(1.2 * rho_s) / np.log(10.) + 0.001, 6., 1e2)
#Re_Rto1WI = 3.91 * We_Rto1WI / ((We_Rto1WI / rho_s - 1.2)**1.11)
#plt.loglog(We_Rto1WI, Re_Rto1WI, marker=None, color='k', linewidth=0.8, linestyle='--') # R to 1WI
# TODO: Note in review that the We_g < 0.4 case is always satisfied for water at atmospheric conditions.

plt.axvline(0.4 * rho_s, marker=None, color='k', linewidth=0.8, linestyle='--') # R to 1WI

plt.axvline(13. * rho_s, marker=None, color='k', linewidth=0.8, linestyle='--') # 1WI to 2WI
plt.axvline(40.3 * rho_s, marker=None, color='k', linewidth=0.8, linestyle='--') # 2WI to atomization

plt.text(2.9e0, 1.e3, 'dripping', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
plt.text(5.9e1, 1.e3, 'Rayleigh', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(1.9e3, 1.e3, 'first wind-induced', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
plt.text(2.e4, 1.e3, 'second wind-induced', backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
plt.text(1.9e5, 1.e3, 'atomization', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.xlabel(r'$\mathrm{We}_{\ell0}$')
plt.ylabel(r'$\mathrm{Re}_{\ell0}$')
plt.grid()
plt.xlim([1.e0, 1.e6])
plt.ylim([1.e0, 1.e6])
axes = plt.gca()
axes.set_xscale('log')
axes.set_yscale('log')
if not(revno is None):
   box1 = TextArea('Rev. '+revno, textprops=dict(color="k"))
   anchored_box = AnchoredOffsetbox(loc=4,
                                    child=box1,
                                    frameon=False,
                                    bbox_to_anchor=(1., -0.20),
                                    bbox_transform=axes.transAxes)
   axes.add_artist(anchored_box)
fig = plt.gcf()
fig.set_size_inches(6., 4., forward=True) # report
plt.savefig('../outputs/figures/old_regime_map_new_coords.png')
plt.savefig('../outputs/figures/old_regime_map_new_coords.pgf', bbox_inches="tight")
fig.set_size_inches(5.5, 3., forward=True) # paper
plt.savefig('../outputs/figures/old_regime_map_new_coords_paper.pgf', bbox_inches="tight")
plt.close()

Re_0_arr   = np.logspace(0., 6., 1e2)
Oh_ItoII   = 135. * Re_0_arr**(-1.27)
Oh_IItoIII = 741. * Re_0_arr**(-1.22)
plt.loglog(Re_0_arr, Oh_ItoII, color='k', linewidth=0.8, linestyle='--')
plt.loglog(Re_0_arr, Oh_IItoIII, color='k', linewidth=0.8, linestyle='--')
plt.text(1.e2, 1.e-1, 'I', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(7.e2, 1.e-1, 'II', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.text(5.e3, 1.e-1, 'III', backgroundcolor='w', verticalalignment='center', horizontalalignment='center')
plt.xlabel(r'$\mathrm{Re}_{\ell0}$')
plt.ylabel(r'$\mathrm{Oh}_{\ell0}$')
plt.grid()
plt.xlim([1.e0, 1.e6])
plt.ylim([1.e-3, 1.e1])
fig = plt.gcf()
fig.set_size_inches(6., 4., forward=True) # report
plt.savefig('../outputs/figures/ohnesorge_diagram.png')
plt.savefig('../outputs/figures/ohnesorge_diagram.pgf', bbox_inches="tight")
fig.set_size_inches(5.5, 3., forward=True) # paper
plt.savefig('../outputs/figures/ohnesorge_diagram_paper.pgf', bbox_inches="tight")
plt.close()

## DATA DEBUGGING:

## dripping points
#regime_df_point = regime_all_df[regime_all_df['We_l0'] < 4.]

#print(regime_df_point['key'])
##print(regime_df_point['photo filename'])
##print(regime_df_point['L_b/d_0 page fig'])

## identify 1WI points with photos
#print "First wind-induced points:"
#regime_all_df_1WI = regime_all_df[regime_all_df['regime combined'] == 'first wind-induced']
#regime_all_df_1WI = regime_all_df_1WI[regime_all_df_1WI['photo filename'].notnull()]
#print regime_all_df_1WI['photo filename']
#print regime_all_df_1WI['regime turb']

#print "Bursting points:"
#regime_all_df_bursting = regime_all_df_with_trans[regime_all_df_with_trans['regime combined'] == 'bursting']
#print regime_all_df_bursting['photo filename']
#print 'Re_l0:'
#print regime_all_df_bursting['Re_l0']
#print 'We_l0:'
#print regime_all_df_bursting['We_l0']
#regime_all_df_bursting = regime_all_df_with_eisenklam[regime_all_df_with_eisenklam['regime combined'] == 'bursting']
#regime_all_df_bursting = regime_all_df_bursting[regime_all_df_bursting['key'] == 'eisenklam_flow_1958']
#print regime_all_df_bursting['photo filename']
#print 'Re_l0:'
#print regime_all_df_bursting['Re_l0']
#print 'We_l0:'
#print regime_all_df_bursting['We_l0']

##print
##print "Near-transitional TSB points:"
##regime_all_df_TSB = regime_all_df[regime_all_df['regime combined'] == 'turbulent surface breakup']
###print len(regime_all_df_TSB)
##regime_near_transitional_TSB = regime_all_df_TSB[regime_all_df_TSB['Re_l0'] > 3.e3]
##regime_near_transitional_TSB = regime_near_transitional_TSB[regime_near_transitional_TSB['Re_l0'] < 1.e4]
##regime_near_transitional_TSB = regime_near_transitional_TSB[regime_near_transitional_TSB['We_l0'] > 3.e4]
##print regime_near_transitional_TSB['key']
##print regime_near_transitional_TSB['liquid']

#print
#print "Possibly misclassified TSB points:"
#regime_all_df_atomization = regime_all_df[regime_all_df['regime photo'] == 'second wind-induced']
#regime_all_df_atomization = regime_all_df_atomization[regime_all_df_atomization['We_l0'] > 2.e4]
#print regime_all_df_atomization['key']
##print regime_all_df_atomization['liquid']
##print regime_all_df_atomization['L_0/d_0']
#print regime_all_df_atomization['photo filename']

macros_regime.close()
