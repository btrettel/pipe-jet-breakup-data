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
from scipy.optimize import fsolve
#import scipy.stats

#revno = lastchangedrevision[0:6]
revno = None

def We_crit_func(We_lo, *data):
   Tubar_0, Re_lo, C = data
   return We_lo**(3./2.) + 3. * We_lo**2. / Re_lo - C * Tubar_0**(-3.)

print

macros_TR = open('../outputs/macros/turbulent_rayleigh.tex', 'w')

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

TR_df = df_jet_breakup
TR_df = TR_df[TR_df['regime L_b'] == 'Rayleigh']
TR_df = TR_df[TR_df['regime turb'] == 'turbulent']
TR_df = TR_df[TR_df['Re_l0'] > Re_turb]
TR_df = TR_df[TR_df['key'] != 'mansour_effect_1994'] # Removing because it is inconsistent with the others. Likely due to not using the same definition of the breakup length.

# C_TR_arr = TR_df['L_b/d_0'] / Rayleigh_Weber_half_pow

# A = np.column_stack([np.ones(len(TR_df)), np.log(TR_df['We_l0'])])
# B = np.log(C_TR_arr)

# result, _, _, _ = np.linalg.lstsq(A, B)
# a, alpha_We_l0_C_TR = result

# C_TR_0 = exp(a)

# print 'C_TR_0 =', C_TR_0
# print 'alpha_We_l0_C_TR =', alpha_We_l0_C_TR
# #print 'alpha_Re_l0_C_TR =', alpha_Re_l0_C_TR

# #C_TR_predicted = C_TR_0 * TR_df['We_l0']**alpha_We_l0_C_TR * TR_df['Re_l0']**alpha_Re_l0_C_TR
# C_TR_predicted = C_TR_0 * TR_df['We_l0']**alpha_We_l0_C_TR

# #print 'C_TR R^2 =', coeff_of_determination(C_TR_predicted, C_TR_arr)

#Rayleigh_Weber_half_pow = TR_df['We_l0']**(1./2.) + 3. * TR_df['We_l0'] / TR_df['Re_l0']
Rayleigh_Weber_half_pow = TR_df['We_l0']**(1./2.)
C_TR = Rayleigh_Weber_half_pow.dot(TR_df['L_b/d_0']) / Rayleigh_Weber_half_pow.dot(Rayleigh_Weber_half_pow)

# TR_RHS = 1./np.sinh(TR_df['L_b/d_0'] / TR_df['We_l0']**(1./2.))
# f      = friction_factor_smooth_array(TR_df['Re_l0'])
# TR_LHS = (f*TR_df['We_l0']/8.)**(1./2.)
# C_v    = TR_LHS.dot(TR_RHS)/TR_LHS.dot(TR_LHS)

TR_RHS = 1./np.sinh(TR_df['L_b/d_0'] / TR_df['We_l0']**(1./2.))
TR_LHS = TR_df['I_0']*TR_df['We_l0']**(1./2.)
C_v    = TR_LHS.dot(TR_RHS)/TR_LHS.dot(TR_LHS)

print 'C_v =', C_v

summary_table(TR_df)

avgTubar_0 = np.average(TR_df['I_0'])
print 'Tubar =', avgTubar_0
print 'old ln(a/delta_0) model =', np.log((2./3.) * avgTubar_0**(-2.))
# LATER: Would be nice to include the uncertainty from Tubar being estimated.

# mansour_effect_1994 p. 134: "All breakup length data reported an average over 20 frames or more." The standard deviation of the breakup length is estimated from phinney_breakup_1973.

coefficient_Mansour = 4.8
print
print 'ln(a/delta_0) Mansour      =', coefficient_Mansour
print
print 'ln(a/delta_0) regression   =', C_TR

delta_0_s = 1./(2.*np.cosh(C_TR))
delta_0_s = 1./(2.*np.exp(C_TR))
print 'delta_0/d_0 regression =', format(delta_0_s, ".2e")

# WON'T: Add statistical error from R^2

#TR_df['L_b/d_0 predicted'] = np.log((2./3.) * TR_df['I_0']**(-2.)) * (TR_df['We_l0']**(1./2.) + 3. * TR_df['We_l0'] / TR_df['Re_l0'])
#TR_df['L_b/d_0 predicted'] = C_TR * TR_df['We_l0']**(1./2.)
const_CTR_R2 = coeff_of_determination(C_TR * TR_df['We_l0']**(1./2.), TR_df['L_b/d_0'])

#TR_df['L_b/d_0 predicted'] = C_TR * (TR_df['We_l0']**(1./2.) + 3. * TR_df['We_l0'] / TR_df['Re_l0'])
TR_df['L_b/d_0 predicted'] = np.arcsinh(1./(C_v * TR_df['I_0'] * TR_df['We_l0']**(1./2.))) * TR_df['We_l0']**(1./2.)
#TR_df['L_b/d_0 predicted'] = np.arcsinh(1./(C_v * (f*TR_df['We_l0']/8.)**(1./2.))) * TR_df['We_l0']**(1./2.)
print 'R^2 =', coeff_of_determination(TR_df['L_b/d_0 predicted'], TR_df['L_b/d_0'])
plot_with_keys(TR_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_turbulent_Rayleigh')

with open('../outputs/data/TR_xbavg.pickle', 'w') as f:
   pickle.dump([C_TR, C_v], f)

macros_TR.write(r'\newcommand{\CTRtheory}{'+roundstr(np.log((2./3.) * avgTubar_0**(-2.)))+'}\n')
macros_TR.write(r'\newcommand{\CTRnum}{'+roundstr(C_TR)+'}\n')
macros_TR.write(r'\newcommand{\CvTR}{'+roundstr(C_v)+'}\n')
macros_TR.write(r'\newcommand{\CTRrsquared}{'+roundstr(coeff_of_determination(TR_df['L_b/d_0 predicted'], TR_df['L_b/d_0']))+'}\n')
macros_TR.write(r'\newcommand{\constCTRrsquared}{'+roundstr(const_CTR_R2)+'}\n')
macros_TR.write(r'\newcommand{\CTRN}{\num{'+str(len(TR_df))+'}}\n\n')

print 'R^2 (TSB regression) =', coeff_of_determination(breakup_length(TR_df['I_0'], TR_df['We_l0']), TR_df['L_b/d_0'])

# TODO: Determine \Welocrit from data.

i = 0
combined_regime_array = []
for regime_photo, regime_turb in zip(df_jet_breakup['regime photo'], df_jet_breakup['regime turb']):
   if isinstance(regime_photo, basestring):
      #print i, regime_photo
      if regime_photo == 'bursting':
         regime_photo = 'F2S'
      
      if (regime_photo == 'Rayleigh') and (regime_turb == 'transitional'):
         regime_photo = 'Rayleigh (transitional)'
      
      if (regime_photo == 'Rayleigh') and (regime_turb == 'turbulent'):
         regime_photo = 'Rayleigh (turbulent)'
      
      combined_regime_array.append(regime_photo)
   else:
      combined_regime_array.append(np.nan)
   
   i = i + 1

i = 0
for regime_L_b, regime_turb in zip(df_jet_breakup['regime L_b'], df_jet_breakup['regime turb']):
   if isinstance(regime_L_b, basestring):
      if (regime_L_b == 'Rayleigh') and (regime_turb == 'transitional'):
         regime_L_b = 'Rayleigh (transitional)'
         
      if (regime_L_b == 'Rayleigh') and (regime_turb == 'turbulent'):
         regime_L_b = 'Rayleigh (turbulent)'
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

TR_df = df_jet_breakup
TR_df['regime combined'] = combined_regime_array
TR_df['We_l0_p'] = We_l0_p
TR_df = TR_df[TR_df['regime combined'].notnull()]
#TR_df = TR_df[TR_df['key'] != 'phinney_breakup_1973'] # Removed until I can understand why it seems wrong.
TR_df = TR_df[TR_df['key'] != 'eisenklam_flow_1958'] # Removed as it seems to have atypically high transition Reynolds numbers.
TR_df = TR_df[TR_df['rho_s'] > 500]
#TR_df = TR_df[TR_df['photo filename'].notnull()]
TR_df = TR_df[TR_df['roughness'] == 'smooth']
regime_all_df = TR_df
TR_df = TR_df[TR_df['regime turb'] == 'turbulent']
#TR_df = TR_df[TR_df['I_0'] < 0.06] # reduces the amount of data too much; use smooth pipes only instead even though this removes Chen
regime_array = []
for regime in TR_df['regime combined']:
   if not(regime in regime_array):
      regime_array.append(regime)

print len(TR_df), 'data points with breakup regime classified.'

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
   this_TR_df = TR_df[TR_df['regime combined'] == regime]
   plt.loglog(this_TR_df['We_l0_p'], this_TR_df['Re_l0'], linestyle='None', marker=marker_array[i], label=regime_print)
   
   i = i + 1

We_T_crit = 8.

d_arr = np.array([0.066, 0.191, 0.347]) * 2.54e-2 # m
rho_l = rho_water(33)
nu_l  = nu_water(33)
sigma = sigma_water(33)

Ubar_arr = np.linspace(1., 30., 1e2)

We_0_arr = rho_l * Ubar_arr**2. * d_arr[0] / sigma
We_1_arr = rho_l * Ubar_arr**2. * d_arr[1] / sigma
We_2_arr = rho_l * Ubar_arr**2. * d_arr[2] / sigma

Re_0_arr = Ubar_arr * d_arr[0] / nu_l
Re_1_arr = Ubar_arr * d_arr[1] / nu_l
Re_2_arr = Ubar_arr * d_arr[2] / nu_l

Tu_0_arr = I_fully_developed_array(Re_0_arr)
Tu_1_arr = I_fully_developed_array(Re_1_arr)
Tu_2_arr = I_fully_developed_array(Re_2_arr)

We_T_0_arr = We_0_arr * Tu_0_arr**2.
We_T_1_arr = We_1_arr * Tu_1_arr**2.
We_T_2_arr = We_2_arr * Tu_2_arr**2.

Ubar_crit_0 = sqrt((sigma / (rho_l * d_arr[0])) * (We_T_crit / 0.05**2.))
Ubar_crit_1 = sqrt((sigma / (rho_l * d_arr[1])) * (We_T_crit / 0.05**2.))
Ubar_crit_2 = sqrt((sigma / (rho_l * d_arr[2])) * (We_T_crit / 0.05**2.))
Ubar_crit_1_highturb = sqrt((sigma / (rho_l * d_arr[1])) * (We_T_crit / 0.10**2.))

print 'Ubar_crit_0 =', Ubar_crit_0, 'm/s'
print 'Ubar_crit_1 =', Ubar_crit_1, 'm/s'
print 'Ubar_crit_2 =', Ubar_crit_2, 'm/s'
print 'Ubar_crit_1_highturb =', Ubar_crit_1_highturb, 'm/s'

plt.loglog(We_T_0_arr, Re_0_arr, linestyle='-', marker=None)
plt.loglog(We_T_1_arr, Re_1_arr, linestyle='-', marker=None)
plt.loglog(We_T_2_arr, Re_2_arr, linestyle='-', marker=None)

d_real = 2. * 2.54e-2;
We_real_arr = rho_l * Ubar_arr**2. * d_real / sigma
Re_real_arr = Ubar_arr * d_real / nu_l
Tu_real_arr = I_fully_developed_array(Re_real_arr)
We_T_real_arr = We_real_arr * Tu_real_arr**2.
plt.loglog(We_T_real_arr, Re_real_arr, linestyle='-', marker=None)

plt.axvline(We_T_crit, marker=None, color='k', zorder=-1, linewidth=0.8, linestyle='--', label=r'$\mathrm{We}_\text{T,crit}$')
plt.axhline(Re_trans, marker=None, color='k', zorder=-1, linewidth=0.8, linestyle='--', label=r'approx.\ onset of turbulence')
plt.axhline(Re_turb, marker=None, color='k', zorder=-1, linewidth=0.8, linestyle='--', label=r'approx.\ fully turbulent')
plt.legend(loc='best', fancybox=True, framealpha=0.5, numpoints=1)
#plt.xlabel(r'$\mathrm{We}_{\ell0}$')
plt.xlabel(r'$\mathrm{We}_\text{T}$')
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
axes.legend().set_visible(False)
plt.savefig('../outputs/figures/regime_map_turb_low_atm_density.png')
#axes.legend().set_visible(True)
#fig.set_size_inches(12., 8., forward=True) # poster
#fig.set_size_inches(4., 2.5, forward=True) # talk
fig.set_size_inches(6., 4., forward=True) # report
# https://stackoverflow.com/a/43439132
#plt.legend(bbox_to_anchor=(0., -0.15), loc='upper left', frameon=False,fontsize=11)
plt.savefig('../outputs/figures/regime_map_turb_low_atm_density.pgf', bbox_inches="tight")
plt.close()

Tubar_0 = 0.06
C       = 20.**(3./2.)
Re_arr = np.logspace(np.log(Re_turb)/np.log(10.), 6., 1e2)
We_hat             = C**(2./3.) * Tubar_0**(-2.)
We_crit_approx_arr = We_hat / (1. + 3. * We_hat**(1./2.) / Re_arr)**(2./3.)

We_crit_arr = np.array([])
for Re_lo, We_crit_approx in zip(Re_arr, We_crit_approx_arr):
   We_crit = fsolve(We_crit_func, We_crit_approx, args=(Tubar_0, Re_lo, C))[0]
   We_crit_arr = np.append(We_crit_arr, We_crit)
   #print We_hat, We_crit, We_crit_approx

plt.axvline(We_hat, marker=None, color='k', zorder=-1, linewidth=0.8, linestyle='--')
plt.semilogy(We_crit_arr, Re_arr, linestyle='-', marker=None, color='b')
plt.semilogy(We_crit_approx_arr, Re_arr, linestyle='-', marker=None, color='r')
#plt.loglog(We_crit_approx_arr_2, Re_arr, linestyle='-', marker=None, color='k')
plt.axhline(Re_turb, marker=None, color='k', zorder=-1, linewidth=0.8, linestyle='--', label=r'approx.\ fully turbulent')
#plt.xlim([1.e2,1.e4])
#plt.loglog(We_0_arr, Re_0_arr, linestyle='-', marker=None)
#plt.loglog(We_1_arr, Re_1_arr, linestyle='-', marker=None)
#plt.loglog(We_2_arr, Re_2_arr, linestyle='-', marker=None)
plt.xlabel(r'$\mathrm{We}_{\ell0}$')
plt.ylabel(r'$\mathrm{Re}_{\ell0}$')
plt.grid()
plt.savefig('../outputs/figures/turb_Rayleigh_transition.png')
plt.close()

#x_0 = We_hat**(3./2.)
#alpha = 4./3.
#x = np.linspace(We_crit_arr[0]**(3./2.), We_hat**(3./2.), 1e2)
#plt.plot(x, x**(4./3.), linestyle='-', marker=None, color='b')
#plt.plot(x, x*We_hat**(1./2.), linestyle='-', marker=None, color='r')
#plt.plot(x, (alpha - 1.) * x_0**(alpha - 2.) * x**2. - alpha * x_0**(alpha - 1.), linestyle='-', marker=None, color='k')
#plt.xlabel(r'$x$')
#plt.ylabel(r'$f(x)$')
#plt.grid()
#plt.savefig('../outputs/figures/Rto2WI_function_approximation.png')
#plt.close()

# Analyze vliem_influence_1975 data.
# "Vliem" case is the one with the highest Tu, restriction number 1 with L_0/d_0 = 10.

print

Tu_pipe_centerline  = np.sqrt(0.00224)
Tu_Vliem_centerline = np.sqrt(0.0095)

print Tu_pipe_centerline, Tu_Vliem_centerline, Tu_Vliem_centerline/Tu_pipe_centerline

T_Vliem           = 15 # C, p. 20
Re_l0_Vliem       = 5000.
Ubar_0_Vliem      = 1.25 # m/s, p. 42
sigma_Vliem       = 73e-3 # N/m, p. 42
rho_l_Vliem       = 1000. # kg/m^3, p. 42
d_0_Vliem         = 4.e-3 # m
pipe_u_prime_plus = 2.65546295174 # approximate value at r/r_0 = 0.901313801666
pipe_max_u_prime  = pipe_u_prime_plus * np.sqrt(friction_factor_smooth(Re_l0_Vliem) / 8.)
Vliem_max_u_prime = np.sqrt(0.048)
We_l0_Vliem       = rho_l_Vliem * Ubar_0_Vliem**2. * d_0_Vliem / sigma_Vliem
C_TR_Vliem        = 0.19 / (d_0_Vliem * We_l0_Vliem**0.5)
# Seems that Vliem is among those with a higher value of C_TR.

Tubar_0_pipe      = I_fully_developed(friction_factor_smooth(Re_l0_Vliem))
Tubar_0_Vliem_est = Tubar_0_pipe * (Vliem_max_u_prime / pipe_max_u_prime)
print Tubar_0_pipe, Tubar_0_Vliem_est, Tubar_0_Vliem_est/Tubar_0_pipe

C_TR_pipe  = C_TR_from_file(Tubar_0_pipe, We_l0_Vliem)
C_TR_Vliem = C_TR_from_file(Tubar_0_Vliem_est, We_l0_Vliem)
print C_TR_pipe, C_TR_Vliem, C_TR_Vliem/C_TR_pipe

macros_TR.write(r'\newcommand{\Vliemcenterlinefactor}{'+roundstr(Tu_Vliem_centerline/Tu_pipe_centerline, 2)+'}\n')
macros_TR.write(r'\newcommand{\Vliempeakfactor}{'+roundstr(Vliem_max_u_prime / pipe_max_u_prime, 2)+'}\n')
macros_TR.write(r'\newcommand{\Vliembarfactor}{'+roundstr(Tubar_0_Vliem_est/Tubar_0_pipe, 2)+'}\n')
macros_TR.write(r'\newcommand{\VliemCTRfactor}{'+roundstr(C_TR_Vliem/C_TR_pipe, 2)+'}\n')

# TODO: Fix this hack that gets citations working in the legends.
os.system("cd ../outputs/figures/ ; for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")

macros_TR.close()
