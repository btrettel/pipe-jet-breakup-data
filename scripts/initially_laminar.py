#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# pipe-jet-breakup-data is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

from jetbreakup import *

#revno = lastchangedrevision[0:6]
revno = None

macros_initially_laminar = open('../outputs/macros/initially_laminar.tex', 'w')

print

# TODO: Calculate Ohnesorge number range for laminar Rayleigh data to compare against those who claim C_LR is a function of the Ohnesorge number.
# TODO: Write C_LR and Re_x_tr to a pickle file and use those values instead of hardcoding.
# TODO: Compare estimated downstream transition breakup length to data.

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

# configuration
C_TR = 3.27

LR_df       = df_jet_breakup
photos_df   = df_jet_breakup
#xbavgs_df   = df_jet_breakup
boundary_df = df_jet_breakup

LR_df = LR_df[LR_df['regime L_b'] == 'Rayleigh']
LR_df = LR_df[LR_df['regime turb'] == 'laminar']
#LR_df['Oh_l0'] = LR_df['We_l0']**(1./2.) / LR_df['Re_l0']
#LR_df = LR_df[LR_df['Oh_l0'] > 0.01]

Rayleigh_Weber_half_pow = LR_df['We_l0']**(1./2.) + 3. * LR_df['We_l0'] / LR_df['Re_l0']
C_LR = Rayleigh_Weber_half_pow.dot(LR_df['L_b/d_0']) / Rayleigh_Weber_half_pow.dot(Rayleigh_Weber_half_pow)

#C_LR = 13.4

print 'C_LR = ', C_LR

LR_df['L_b/d_0 predicted'] = C_LR * (LR_df['We_l0']**(1./2.) + 3. * LR_df['We_l0'] / LR_df['Re_l0'])
plot_with_keys(LR_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_laminar_Rayleigh')
print 'R^2 =', coeff_of_determination(LR_df['L_b/d_0 predicted'], LR_df['L_b/d_0'])
# TODO: Note that predicted vs. actual plot error does not include statistical error.

i = 0
combined_regime_array = []
for regime_photo, regime_turb in zip(df_jet_breakup['regime photo'], df_jet_breakup['regime turb']):
   if isinstance(regime_photo, basestring):
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

photos_df['regime combined'] = combined_regime_array
photos_df = photos_df[photos_df['regime combined'] == 'downstream transition']
photos_df = photos_df[photos_df['key'] != 'eisenklam_flow_1958'] # Removed as it seems to have atypically high transition Reynolds numbers.
photos_df = photos_df[photos_df['x_trans/d_0'].notnull()]

#summary_table(photos_df)

#print photos_df['photo filename']
#print photos_df['x_trans/d_0']
Re_x_tr_arr = photos_df['Re_l0'] * photos_df['x_trans/d_0']
#print Re_x_tr_arr
Re_x_tr = np.average(Re_x_tr_arr)
print 'Re_x_tr =', Re_x_tr
print 'min Re_x_tr =', np.min(Re_x_tr_arr)
print 'max Re_x_tr =', np.max(Re_x_tr_arr)

boundary_df         = boundary_df[boundary_df['regime L_b'] == 'R2F']
boundary_df         = boundary_df[boundary_df['rho_s'] < 1500.]
boundary_df         = boundary_df[boundary_df['key'] != 'sterling_instability_1975']
Re_x_tr_implied_arr = C_LR * (boundary_df['Re_l0'] * boundary_df['We_l0']**(1./2.) + 3. * boundary_df['We_l0'])
Re_x_tr_implied = np.average(Re_x_tr_implied_arr)
print 'Re_x,tr,implied =', Re_x_tr_implied
print 'min Re_x,tr,implied =', np.min(Re_x_tr_implied_arr)
print 'max Re_x,tr,implied =', np.max(Re_x_tr_implied_arr)
boundary_df['Re_x,tr implied'] = Re_x_tr_implied_arr
#print
#print Re_x_tr_implied / Re_x_tr

# max_boundary_df = boundary_df[boundary_df['Re_x,tr implied'] == np.max(Re_x_tr_implied_arr)]
# print max_boundary_df['key']
# print max_boundary_df['L_b/d_0']
# print max_boundary_df['liquid']

# # start 2020-02-16 breakup length equation approximation
# We_l0_asympt = Re_x_tr / (3. * C_DTpeak)
# print 'We_l0_asympt =', We_l0_asympt

# # pure isopropyl alcohol
# d_0   = 1.e-3 # m
# rho_l = 0.786 * rho_water(T_std) # kg/m^3
# mu_l  = 2.036e-3 # Pa*s?
# nu_l  = mu_l / rho_l # m^2/s
# sigma = 21.74e-3 # N/m

# def Ubar_0_crit_func(Ubar_0):
   # return C_LR * ((rho_l * Ubar_0**2. * d_0 / sigma)**(1./2.) + 3. * rho_l * nu_l * Ubar_0 / sigma) - Re_x_tr * nu_l / (Ubar_0 * d_0)

# Ubar_0_crit = fsolve(Ubar_0_crit_func, 1.5)[0]
# print 'Ubar_0,crit =', Ubar_0_crit, 'm/s'

# Re_l0_crit = Ubar_0_crit * d_0 / nu_l
# We_l0_crit = rho_l * Ubar_0_crit**2. * d_0 / sigma
# Oh_l0      = mu_l / np.sqrt(rho_l * sigma * d_0)

# print 'Re_l0,crit =', Re_l0_crit
# print 'We_l0,crit =', We_l0_crit

# Re_l0_arr = np.linspace(Re_l0_crit, Re_trans, 5e2)
# We_l0_arr = (rho_l * nu_l**2.) / (sigma * d_0) * Re_l0_arr**2.

# xbavgs_LR_arr = C_LR * (We_l0_arr**(1./2.) + 3. * We_l0_arr / Re_l0_arr)
# x_trs_arr  = Re_x_tr / Re_l0_arr
# xbavgs_arr = x_trs_arr + np.log(1./(np.exp(C_LR * (x_trs_arr / xbavgs_LR_arr - 1.)) + np.exp(-C_TR) - np.exp(-C_LR))) * (We_l0_arr**(1./2.) + 3. * We_l0_arr / Re_l0_arr)

# Re_l0_peak = np.sqrt(Re_x_tr / (C_DTpeak * Oh_l0 * (1. + 3. * Oh_l0)))

# plt.plot(Re_l0_arr, xbavgs_arr, linestyle='-')
# plt.axvline(Re_l0_peak, linestyle='--')
# plt.xlabel(r'$\mathrm{Re}_{\ell0}$')
# plt.ylabel(r'$\langle x_\text{b} \rangle / d_0$')
# plt.grid()
# #plt.xlim([0., 30000])
# #plt.ylim([0., 1.])
# #plt.legend()
# fig = plt.gcf()
# fig.set_size_inches(6., 4., forward=True) # report
# plt.savefig('../outputs/figures/downstream_transition_xbavg.png')
# # plt.savefig('../outputs/figures/downstream_transition_xbavg.pgf', bbox_inches="tight")
# # fig.set_size_inches(5., 3., forward=True) # paper
# # plt.savefig('../outputs/figures/downstream_transition_xbavg_paper.pgf', bbox_inches="tight")
# plt.close()
# # end 2020-02-16 breakup length equation approximation

macros_initially_laminar.close()

# TODO: Unfortunately determining where transition occurred was difficult and the transition locations are not expected to be . Some photos had poor resolution. The process of determining where the flow transitioned also was fairly subjective. The transition location also was only found from photographs representing a single point in time. It is likely that the transition location fluctuates in time and this fluctuation does not appear in the photos.

# TODO: Make a regression for \xbavg in the downstream transition regime

# TODO: Fix this hack that gets citations working in the legends.
os.system("cd ../outputs/figures/ ; for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")
