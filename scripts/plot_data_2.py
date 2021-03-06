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

macros_reg = open('../outputs/macros/regression_2.tex', 'w')

print

if not(os.path.exists('../outputs/data/atomization_boundary.pickle')):
   print 'Run atomization.py to get the atomization regime boundary first.'
   exit(-1)

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

# # # # # # # # # # # # # # # # # # # # ##
# Comparison against data not from pipes #
## # # # # # # # # # # # # # # # # # # # #

# Issues with this data:
# 1. Divergence from previous correlation for high \langle x_b \rangle / d_0. This could be caused by a few things. Most notable are the exponent on Tu for the \langle x_b \rangle correlation. If the exponent is in the range of about -0.25 to -0.35 then the data matches the correlation much better (the correlation runs through the center of the data, although the data has large spread), but this is inconsistent with the most credible estimate for Tu that I can make for kusui_liquid_1969. Another possibility is that this data is in the atomization regime. However, the values of TA would suggest they are in the second wind induced regime. Yet another possibility is that the integral scales for this data are far off that for pipe jets.

ervine_effect_1980_csv = pd.read_csv('../data/ervine_effect_1980/ervine_effect_1980_fig_8_2.csv', sep=',', header=0)

d_0_ervine_effect_1980    = ervine_effect_1980_csv['d (m)']
I_0_ervine_effect_1980    = ervine_effect_1980_csv['TI']
Ubar_0_ervine_effect_1980 = ervine_effect_1980_csv['v (m/s)']
L_bs_ervine_effect_1980   = ervine_effect_1980_csv['L_b/d']

rho_l = rho_water(T_std)
sigma = sigma_water(T_std)
nu_l  = nu_water(T_std)
#P_v   = P_v_water(T_std)
rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g

We_l0_ervine_effect_1980 = rho_l * Ubar_0_ervine_effect_1980**2 * d_0_ervine_effect_1980 / sigma
Re_l0_ervine_effect_1980 = Ubar_0_ervine_effect_1980 * d_0_ervine_effect_1980 / nu_l
#Fr_0_ervine_effect_1980  = Ubar_0**2 / (g * d_0)
rho_s_ervine_effect_1980 = rho_l / rho_g
#nu_s_ervine_effect_1980  = nu_l / nu_g
#K_c_ervine_effect_1980   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
#Ma_g_ervine_effect_1980  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)

df_ervine_effect_1980 = pd.DataFrame({'We_l0': We_l0_ervine_effect_1980})

df_ervine_effect_1980['key'] = 'ervine_effect_1980'

df_ervine_effect_1980['L_b/d_0']          = L_bs_ervine_effect_1980
df_ervine_effect_1980['L_b/d_0 page fig'] = 'p. 436, fig. 8'

#df_ervine_effect_1980['We_l0']     = We_l0_ervine_effect_1980
df_ervine_effect_1980['Re_l0']     = Re_l0_ervine_effect_1980
df_ervine_effect_1980['I_0']       = I_0_ervine_effect_1980
#df_ervine_effect_1980['Fr_0']      = Fr_0_ervine_effect_1980
df_ervine_effect_1980['rho_s']     = rho_s_ervine_effect_1980
#df_ervine_effect_1980['nu_s']      = nu_s_ervine_effect_1980
#df_ervine_effect_1980['K_c']       = K_c_ervine_effect_1980
#df_ervine_effect_1980['Ma_g']      = Ma_g_ervine_effect_1980
#df_ervine_effect_1980['est rho_s'] = True
#df_ervine_effect_1980['est nu_s']  = True

df_other_comparison = df_ervine_effect_1980

mckeogh_air_1980_csv = pd.read_csv('../data/mckeogh_air_1980/mckeogh_air_1980_fig_3.csv', sep=',', header=0)

d_0_mckeogh_air_1980  = mckeogh_air_1980_csv['d (m)']
I_0_mckeogh_air_1980  = mckeogh_air_1980_csv['TI']
Q_mckeogh_air_1980    = mckeogh_air_1980_csv['Q (x1e3 m^3/s)'] * 1e-3
L_bs_mckeogh_air_1980 = mckeogh_air_1980_csv['L_b (m)'] / d_0_mckeogh_air_1980

A_mckeogh_air_1980      = (np.pi / 4) * d_0_mckeogh_air_1980**2
Ubar_0_mckeogh_air_1980 = Q_mckeogh_air_1980 / A_mckeogh_air_1980

We_l0_mckeogh_air_1980 = rho_l * Ubar_0_mckeogh_air_1980**2 * d_0_mckeogh_air_1980 / sigma
Re_l0_mckeogh_air_1980 = Ubar_0_mckeogh_air_1980 * d_0_mckeogh_air_1980 / nu_l
#Fr_0_mckeogh_air_1980  = Ubar_0**2 / (g * d_0)
rho_s_mckeogh_air_1980 = rho_l / rho_g
#nu_s_mckeogh_air_1980  = nu_l / nu_g
#K_c_mckeogh_air_1980   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
#Ma_g_mckeogh_air_1980  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)

df_mckeogh_air_1980 = pd.DataFrame({'We_l0': We_l0_mckeogh_air_1980})

df_mckeogh_air_1980['key'] = 'mckeogh_air_1980'

df_mckeogh_air_1980['L_b/d_0']          = L_bs_mckeogh_air_1980
df_mckeogh_air_1980['L_b/d_0 page fig'] = 'p. 1583, fig. 3'

#df_mckeogh_air_1980['We_l0']     = We_l0_mckeogh_air_1980
df_mckeogh_air_1980['Re_l0']     = Re_l0_mckeogh_air_1980
df_mckeogh_air_1980['I_0']       = I_0_mckeogh_air_1980
#df_mckeogh_air_1980['Fr_0']      = Fr_0_mckeogh_air_1980
df_mckeogh_air_1980['rho_s']     = rho_s_mckeogh_air_1980
#df_mckeogh_air_1980['nu_s']      = nu_s_mckeogh_air_1980
#df_mckeogh_air_1980['K_c']       = K_c_mckeogh_air_1980
#df_mckeogh_air_1980['Ma_g']      = Ma_g_mckeogh_air_1980
#df_mckeogh_air_1980['est rho_s'] = True
#df_mckeogh_air_1980['est nu_s']  = True

df_other_comparison = pd.concat([df_other_comparison, df_mckeogh_air_1980])

df_other_comparison['We_l0_crit'] = We_l0_crit(df_other_comparison['I_0'], df_other_comparison['rho_s'])
df_other_comparison = df_other_comparison[df_other_comparison['We_l0'] < df_other_comparison['We_l0_crit']] # Exclude atomization points.

#df_other_comparison = df_other_comparison[df_other_comparison['Re_l0'] > 1.e5] # Use only likely turbulent points.
#df_other_comparison = df_other_comparison[df_other_comparison['I_0'] >= 0.01] # Use only likely turbulent points.
#df_other_comparison['We_l0_crit_TR'] = We_l0_crit_TR(df_other_comparison['I_0'])
#df_other_comparison = df_other_comparison[df_other_comparison['We_l0'] > df_other_comparison['We_l0_crit_TR']]

with open('../outputs/data/TSB_xbavg.pickle') as f:
   C_TSB, alpha_We, alpha_Re_l0_2WI, alpha_Tu_2WI, alpha_rho_s_2WI = pickle.load(f)

df_other_comparison['L_b/d_0 predicted'] = C_TSB * df_other_comparison['We_l0']**alpha_We * df_other_comparison['Re_l0']**alpha_Re_l0_2WI * df_other_comparison['I_0']**alpha_Tu_2WI * df_other_comparison['rho_s']**alpha_rho_s_2WI
plot_with_keys(df_other_comparison, 'other comparison', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno)

macros_reg.write(r'\newcommand{\xbavgregaltrsquared}{'+roundstr(coeff_of_determination(df_other_comparison['L_b/d_0 predicted'], df_other_comparison['L_b/d_0']))+'}\n')
macros_reg.write(r'\newcommand{\xbavgregaltN}{\num{'+str(len(df_other_comparison))+'}}\n')
macros_reg.write(r'\newcommand{\xbavgTurangealt}{$\num{'+str(round(100.*np.amin(df_other_comparison['I_0']), 1))+r'}\% \leq \Tu_\text{c0} \leq \num{'+str(round(100.*np.amax(df_other_comparison['I_0']), 1))+'}\%$}\n')
macros_reg.write('\n')

# Reasons for poor fit:
# 1. Inaccuracy of TI measurement due to it being a pressure measurement.
# 2. Centerline TI is not equal to plane average TI.
# 3. Integral scale and other variables are not constant between data points, even for the same nozzle. Others: velocity profile/boundary layer thickness.
# 4. Reynolds number not high enough for turbulent flow. The trend is fairly different for Re > 1e5 only, but this data set has artificially enhanced turbulence, so likely Re_crit is fairly low.
# 5. Removing the lowest turbulence intensity data did not seem to help.

#######################
# spray angle (theta) #
#######################

# TODO: Check for air entrainment effects on \theta. t/d_i, distance, wall, coflow

theta_df = df_jet_breakup[df_jet_breakup['theta'].notnull()]
theta_df = theta_df[theta_df['regime turb'] == 'turbulent']
#theta_df = theta_df[theta_df['Re_l0'] > 5000]
theta_df = theta_df[theta_df['Ma_g'] < 0.3]
#theta_df = theta_df[theta_df['rho_s'] < 500]
#theta_df = theta_df[theta_df['key'] == 'arai_break-up_1985']
#theta_df = theta_df[theta_df['key'] == 'skrebkov_turbulentnyye_1963']
#theta_df = theta_df[theta_df['rho_s'] > 500]

macros_reg.write(r'\newcommand{\thetaikeysall}{\citep{')
key_array = []
for key in theta_df['key']:
   if not(key in key_array):
      key_array.append(key)
      if len(key_array) == 1:
         macros_reg.write(key)
      else:
         macros_reg.write(','+key)
macros_reg.write('}}\n')

theta_df = theta_df[theta_df['key'] != 'reitz_atomization_1978']
theta_df = theta_df[theta_df['key'] != 'arai_break-up_1985']
theta_df = theta_df[theta_df['key'] != 'rupe_dynamic_1962']
theta_df = theta_df[theta_df['key'] != 'grant_newtonian_1965']
theta_df = theta_df[theta_df['key'] != 'hoyt_pipe-exit_1980']
theta_df = theta_df[theta_df['key'] != 'hoyt_effect_1985']
theta_df = theta_df[theta_df['key'] != 'wu_liquid_1992']
#theta_df_1 = theta_df_1[theta_df_1['key'] != 'skrebkov_turbulentnyye_1963']

theta_df_1 = theta_df
theta_df_1['We_l0_crit'] = We_l0_crit(theta_df_1['I_0'], theta_df_1['rho_s'])
theta_df_1 = theta_df[theta_df['We_l0'] < theta_df['We_l0_crit']] # Exclude atomization points.

theta_df_2 = theta_df[theta_df['key'] == 'skrebkov_turbulentnyye_1963']

parameter_space_plots(theta_df, 'theta', revno=revno)

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Weber number only (smooth pipes only), TSB regime #
## # # # # # # # # # # # # # # # # # # # # # # # # ##

latex_summary_table(theta_df, 'theta')

A = np.column_stack([np.ones(len(theta_df_2)), np.log(theta_df_2['I_0'])])
tan_theta = np.tan(theta_df_2['theta']/2.)
B = np.log(tan_theta)

result, _, _, _ = np.linalg.lstsq(A, B)
a, alpha_Tubar_0_theta = result

A = np.column_stack([np.ones(len(theta_df_1)), np.log(theta_df_1['We_l0'])])
tan_theta = np.tan(theta_df_1['theta']/2.)
B = np.log(tan_theta * theta_df_1['I_0']**(-alpha_Tubar_0_theta))

result, _, _, _ = np.linalg.lstsq(A, B)
a, alpha_We_l0_theta = result

C_theta = exp(a)

print 'C_theta = ', C_theta
print 'alpha_We_l0_theta = ', alpha_We_l0_theta
print 'alpha_Tubar_0_theta = ', alpha_Tubar_0_theta

theta_df_1['theta predicted'] = 2. * np.arctan(C_theta * theta_df_1['We_l0']**alpha_We_l0_theta * theta_df_1['I_0']**alpha_Tubar_0_theta)
tan_theta_predicted = np.tan(theta_df_1['theta predicted'] / 2.)
print 'tan theta/2 R^2 =', coeff_of_determination(tan_theta_predicted, tan_theta)
plot_with_keys(theta_df_1, 'correlation', 'theta predicted', 'theta', plot_type='linear', add_line=True, revno=revno)
print 'theta R^2 =', coeff_of_determination(theta_df_1['theta predicted'], theta_df_1['theta'])
print

with open('../outputs/data/TSB_thetai.pickle', 'w') as f:
   pickle.dump([C_theta, alpha_We_l0_theta, alpha_Tubar_0_theta], f)

macros_reg.write(r'\newcommand{\thetaikeysused}{\citep{')
key_array = []
for key in theta_df_1['key']:
   if not(key in key_array):
      key_array.append(key)
      if len(key_array) == 1:
         macros_reg.write(key)
      else:
         macros_reg.write(','+key)
macros_reg.write('}}\n')

#macros_reg.write(r'\newcommand{\thetaireg}{\tan\left(\frac{\thetai}{2}\right) = '+roundstr(C_theta)+r' \Tubarexp{'+roundstr(C_I_0)+'} \Welo^{'+roundstr(alpha_We_l0_theta)+'}}\n')
macros_reg.write(r'\newcommand{\thetaireg}{\tan\left(\frac{\thetai}{2}\right) = '+roundstr(C_theta)+r' \Tubarexp{'+roundstr(alpha_Tubar_0_theta)+'} \Welo^{'+roundstr(alpha_We_l0_theta)+r'}}'+'\n')
macros_reg.write(r'\newcommand{\thetairegrsquared}{'+roundstr(coeff_of_determination(theta_df_1['theta predicted'], theta_df_1['theta']))+'}\n')
macros_reg.write(r'\newcommand{\thetairegN}{\num{'+str(len(theta_df_1))+'}}\n')
macros_reg.write('\n')

print '\n##########\n'

# # # # # # # # # # # # # # # # # # # # ##
# Comparison against data not from pipes #
## # # # # # # # # # # # # # # # # # # # #

# reitz_atomization_1978 w/ kent_nozzle_1983 (only use non-cavitating points because kent_nozzle_1983 was for air), ervine_behaviour_1987?, balewski_experimental_2008

with open('../outputs/data/atomization_boundary.pickle', 'r') as f:
   C_TSBtoA, alpha_Tu_2WItoA, alpha_rho_s_2WItoA = pickle.load(f)

theta_df = df_jet_breakup
theta_df = theta_df[theta_df['regime L_b'] != 'Rayleigh']
theta_df = theta_df[theta_df['regime photo'].notnull()]
theta_df = theta_df[theta_df['regime photo'] != 'Rayleigh']
theta_df = theta_df[theta_df['regime turb'] == 'turbulent']
theta_df = theta_df[theta_df['Re_l0'] > Re_turb]
theta_df = theta_df[theta_df['theta'].notnull()]

Tubar_0_bar = np.average(theta_df['I_0'])

#print 'C_TSBtoA = ', C_TSBtoA
#print 'C_theta =', C_theta
#print 'alpha_Tu_2WItoA =', alpha_Tu_2WItoA
#print 'alpha_Tubar_0_theta =', alpha_Tubar_0_theta

C_theta_boundary     = C_theta * C_TSBtoA**alpha_We_l0_theta
alpha_theta_boundary = alpha_Tu_2WItoA*alpha_We_l0_theta+alpha_Tubar_0_theta

print 'C_theta_boundary =', C_theta_boundary
print 'alpha_theta_boundary =', alpha_theta_boundary
print

macros_reg.write(r'\newcommand{\thetaiboundary}{\tan\left(\frac{\thetai}{2}\right)_\text{crit} = '+roundstr(C_theta_boundary)+r' \left(\frac{\rhol}{\rhog}\right)^{'+roundstr(alpha_We_l0_theta)+r'} \Tubarexp{'+roundstr(alpha_theta_boundary)+'}}\n')
macros_reg.write('\n')

for theta, Tubar_0, rho_s, regime_photo, photo_filename in zip(theta_df['theta'], theta_df['I_0'], theta_df['rho_s'], theta_df['regime photo'], theta_df['photo filename']):
   tan_theta_crit = C_theta_boundary * rho_s**alpha_We_l0_theta * Tubar_0**alpha_theta_boundary
   tan_theta = np.tan(theta / 2.)
   
   if (tan_theta / tan_theta_crit) <= 0.5:
      regime_theta = 'second wind-induced'
   elif ((tan_theta / tan_theta_crit) > 0.5) and ((tan_theta / tan_theta_crit) <= 1.25):
      regime_theta = 'S2A'
   else:
      regime_theta = 'atomization'
   
   #print tan_theta / tan_theta_crit, tan_theta, tan_theta_crit, regime_photo, regime_theta
   
   if regime_photo != regime_theta:
      print photo_filename, regime_photo, regime_theta

# TODO: Fix this hack that gets citations working in the legends.
os.system("cd ../outputs/figures/ ; for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")

macros_reg.close()
