#!/usr/bin/env python
# -*- coding: utf-8 -*-

from jetbreakup import *

#revno = lastchangedrevision[0:6]
revno = None

macros_reg = open('../outputs/macros/regression.tex', 'w')

# Following mcdermott_quality_2011: "Parameter space" is defined as the entire range of all dependent variables of interest. The "application space" is the subset useful for a particular application.

# TODO: Add error in coefficients and exponents.
# LATER: Change \theta correlation so that \lim_{\rho_g/\rho_l \rightarrow 0} \theta does not equal zero. It currently does, which is not true.

print

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

##########
# photos #
##########

#photo_df = df_jet_breakup
#photo_df = photo_df[photo_df['photo filename'].notnull()]
#parameter_space_plots(photo_df, 'photo filename', revno=revno)

#print len(photo_df), 'photographs.\n'

# TODO: Check if pipe thickness divided by pipe inner diameter (t_0/d_0) correlates with anything as suggested by kerstein_parameter_2017 p. 811 R5-5. Likely more important for higher gas densities.
# TODO: Test if there are any differences due to orientation and/or Froude number.
# TODO: Use nozzle material to see which data points are most likely to not be smooth.
# TODO: Variation of \overline{D} with \nu?
# TODO: Find coefficients for I_0 first with only studies with widely varying TI, then redo correlation taking that as a certainty?

# TODO: add error bars to the regression coefficients

############################
# breakup length (L_b/d_0) #
############################

# TODO: Instead of giving photographic measurement zero weight, give it a small weight?

L_bs_df = df_jet_breakup[df_jet_breakup['L_b/d_0'].notnull()]
L_bs_df = L_bs_df[L_bs_df['regime turb'] == 'turbulent']
L_bs_df = L_bs_df[L_bs_df['Ma_g'] < 0.3]
L_bs_df = L_bs_df[L_bs_df['regime L_b'] == 'second wind-induced']
L_bs_df_more = L_bs_df
parameter_space_plots(L_bs_df_more, 'L_b/d_0', revno=revno)
L_bs_df = L_bs_df[L_bs_df['L_b method'] == 'electrical']
latex_summary_table(L_bs_df_more, 'L_b/d_0')
#parameter_space_plots(L_bs_df, 'L_b/d_0', revno=revno)

key_array = []
for key in L_bs_df_more['key']:
   if not(key in key_array):
      key_array.append(key)

#At = (L_bs_df['rho_s'] - 1) / (L_bs_df['rho_s'] + 1)
A = np.column_stack([np.ones(len(L_bs_df)), log(L_bs_df['We_l0']), log(L_bs_df['Re_l0']), log(L_bs_df['I_0']), log(L_bs_df['rho_s'])])
B = log(L_bs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0, C_Re_l0, C_I_0, C_rho_s = result

#C_L_b = exp(a)

#print 'C_L_b   = ', C_L_b
#print 'C_We_l0 = ', C_We_l0
#print 'C_Re_l0 = ', C_Re_l0
#print 'C_I_0   = ', C_I_0
#print 'C_rho_s = ', C_rho_s
##print 'C_nu_s = ', C_nu_s

macros_reg.write(r'\newcommand{\xbavgReexp}{'+roundstr(C_Re_l0)+'}\n')
macros_reg.write(r'\newcommand{\xbavgrhosexp}{'+roundstr(C_rho_s)+'}\n')
macros_reg.write(r'\newcommand{\xbavgfullWeloexp}{'+roundstr(C_We_l0)+'}\n')

A = np.column_stack([np.ones(len(L_bs_df)), log(L_bs_df['We_l0']), log(L_bs_df['I_0'])])
B = log(L_bs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0, C_I_0 = result

C_xbavg = exp(a)

C_Re_l0 = 0
C_rho_s = 0

print 'C_xbavg = ', C_xbavg
print 'C_We_l0 = ', C_We_l0
print 'C_Re_l0 = ', C_Re_l0
print 'C_I_0   = ', C_I_0
print 'C_rho_s = ', C_rho_s

with open('../outputs/data/TSB_xbavg.pickle', 'w') as f:
   pickle.dump([C_xbavg, C_We_l0, C_Re_l0, C_I_0, C_rho_s], f)

print 'sigma log <x_b>/d_0 =', sigma_regression(a + C_We_l0 * log(L_bs_df['We_l0']) + C_Re_l0 * log(L_bs_df['Re_l0']) + C_I_0 * log(L_bs_df['I_0']) + C_rho_s * log(L_bs_df['rho_s']), log(L_bs_df['L_b/d_0']))

L_bs_df_more['L_b/d_0 predicted'] = C_xbavg * L_bs_df_more['We_l0']**C_We_l0 * L_bs_df_more['Re_l0']**C_Re_l0 * L_bs_df_more['I_0']**C_I_0 * L_bs_df_more['rho_s']**C_rho_s
plot_with_keys(L_bs_df_more, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno)
print 'R^2 =', coeff_of_determination(L_bs_df_more['L_b/d_0 predicted'], L_bs_df_more['L_b/d_0'])

macros_reg.write(r'\newcommand{\xbavgreg}{\frac{\xbavg}{d_0} = '+roundstr(C_xbavg)+r' \Tubarexp{'+roundstr(C_I_0)+'} \Welo^{'+roundstr(C_We_l0)+'}}\n')
macros_reg.write(r'\newcommand{\xbavgregrsquared}{'+roundstr(coeff_of_determination(L_bs_df_more['L_b/d_0 predicted'], L_bs_df_more['L_b/d_0']))+'}\n')
macros_reg.write(r'\newcommand{\xbavgregN}{\num{'+str(len(L_bs_df_more))+'}}\n')
macros_reg.write(r'\newcommand{\xbavgregNelectrical}{\num{'+str(len(L_bs_df))+'}}\n')
macros_reg.write(r'\newcommand{\xbavgregNstudies}{\num{'+str(len(key_array))+'}}\n')
macros_reg.write(r'\newcommand{\xbavgTurange}{$\num{'+str(round(100.*np.amin(L_bs_df_more['I_0']), 1))+r'}\% \leq \Tubar_0 \leq \num{'+str(round(100.*np.amax(L_bs_df_more['I_0']), 1))+'}\%$}\n')
macros_reg.write('\n')

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

if os.path.exists('../outputs/data/atomization_boundary.pickle'):
   df_other_comparison['We_l0_crit'] = We_l0_crit(df_other_comparison['I_0'], df_other_comparison['rho_s'])
   df_other_comparison = df_other_comparison[df_other_comparison['We_l0'] < df_other_comparison['We_l0_crit']] # Exclude atomization points.
else:
   print 'Run atomization.py to get the atomization regime boundary first.'
   exit(-1)

#df_other_comparison = df_other_comparison[df_other_comparison['Re_l0'] > 1.e5] # Use only likely turbulent points.
#df_other_comparison = df_other_comparison[df_other_comparison['I_0'] >= 0.01] # Use only likely turbulent points.
#df_other_comparison['We_l0_crit_TR'] = We_l0_crit_TR(df_other_comparison['I_0'])
#df_other_comparison = df_other_comparison[df_other_comparison['We_l0'] > df_other_comparison['We_l0_crit_TR']]

df_other_comparison['L_b/d_0 predicted'] = C_xbavg * df_other_comparison['We_l0']**C_We_l0 * df_other_comparison['Re_l0']**C_Re_l0 * df_other_comparison['I_0']**C_I_0 * df_other_comparison['rho_s']**C_rho_s
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

print '\n##########\n'

#print df_other_comparison['Re_l0']

#exit()

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

macros_reg.write(r'\newcommand{\thetaikeysall}{\citet{')
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

theta_df['We_l0_crit'] = We_l0_crit(theta_df['I_0'], theta_df['rho_s'])
theta_df = theta_df[theta_df['We_l0'] < theta_df['We_l0_crit']] # Exclude atomization points.

parameter_space_plots(theta_df, 'theta', revno=revno)

## # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
## Hypothesis 1prime: Power law with no density variation #
### # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#latex_summary_table(theta_df, 'theta')
#A = np.column_stack([np.ones(len(theta_df)), log(theta_df['We_l0']), log(theta_df['I_0'])]) # Power law excluding Re.
#tan_theta = np.tan(theta_df['theta']/2.)
#B = log(tan_theta)

#result, _, _, _ = np.linalg.lstsq(A, B)
#a, C_We_l0, C_I_0 = result

#C_theta = exp(a)

#print 'C_theta = ', C_theta
#print 'C_We_l0 = ', C_We_l0
#print 'C_I_0   = ', C_I_0

#theta_var = theta_df['We_l0']**(-1./5.) * theta_df['I_0']**(3./5.)
#C_v_d2 = theta_var.dot(tan_theta) / theta_var.dot(theta_var)

#print 'C_v_d2      = ', C_v_d2

#theta_df['theta predicted'] = 2. * np.arctan(C_theta * theta_df['We_l0']**C_We_l0 * theta_df['I_0']**C_I_0)
#tan_theta_predicted = np.tan(theta_df['theta predicted'] / 2.)
#print 'H1prime (Power law with no density variation) tan theta/2 R^2 =', coeff_of_determination(tan_theta_predicted, tan_theta)
#plot_with_keys(theta_df, 'correlation', 'theta predicted', 'theta', plot_type='linear', add_line=True, revno=revno)
#print 'H1prime (Power law with no density variation) theta R^2 =', coeff_of_determination(theta_df['theta predicted'], theta_df['theta'])
#print

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Weber number only (smooth pipes only), TSB regime #
## # # # # # # # # # # # # # # # # # # # # # # # # ##

latex_summary_table(theta_df, 'theta')
A = np.column_stack([np.ones(len(theta_df)), log(theta_df['We_l0'])])
tan_theta = np.tan(theta_df['theta']/2.)
B = log(tan_theta)

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0 = result

C_theta = exp(a)

print 'C_theta = ', C_theta
print 'C_We_l0 = ', C_We_l0

theta_df['theta predicted'] = 2. * np.arctan(C_theta * theta_df['We_l0']**C_We_l0)
tan_theta_predicted = np.tan(theta_df['theta predicted'] / 2.)
print 'tan theta/2 R^2 =', coeff_of_determination(tan_theta_predicted, tan_theta)
plot_with_keys(theta_df, 'correlation', 'theta predicted', 'theta', plot_type='linear', add_line=True, revno=revno)
print 'theta R^2 =', coeff_of_determination(theta_df['theta predicted'], theta_df['theta'])
print

with open('../outputs/data/TSB_thetai.pickle', 'w') as f:
   pickle.dump([C_theta, C_We_l0], f)

macros_reg.write(r'\newcommand{\thetaikeysused}{\citet{')
key_array = []
for key in theta_df['key']:
   if not(key in key_array):
      key_array.append(key)
      if len(key_array) == 1:
         macros_reg.write(key)
      else:
         macros_reg.write(','+key)
macros_reg.write('}}\n')

#macros_reg.write(r'\newcommand{\thetaireg}{\tan\left(\frac{\thetai}{2}\right) = '+roundstr(C_theta)+r' \Tubarexp{'+roundstr(C_I_0)+'} \Welo^{'+roundstr(C_We_l0)+'}}\n')
macros_reg.write(r'\newcommand{\thetaireg}{\tan\left(\frac{\thetai}{2}\right) = '+roundstr(C_theta)+r' \Welo^{'+roundstr(C_We_l0)+'}}\n')
macros_reg.write(r'\newcommand{\thetairegrsquared}{'+roundstr(coeff_of_determination(theta_df['theta predicted'], theta_df['theta']))+'}\n')
macros_reg.write(r'\newcommand{\thetairegN}{\num{'+str(len(theta_df))+'}}\n')
macros_reg.write('\n')

print '\n##########\n'

# # # # # # # # # # # # # # # # # # # # ##
# Comparison against data not from pipes #
## # # # # # # # # # # # # # # # # # # # #

# reitz_atomization_1978 w/ kent_nozzle_1983 (only use non-cavitating points because kent_nozzle_1983 was for air), ervine_behaviour_1987?, balewski_experimental_2008

######################################
# initial breakup location (x_i/d_0) #
######################################

# TODO: Add comparisons with non-pipe data: reitz_atomization_1978, wu_effects_1995
# TODO: Check against wu_onset_1995's non FDT data.
# TODO: zajac_correlation_1970 has data for variation of x_i with TI? Mentioned on p. 8.
# TODO: We*I_0^3 vs. x_i/d_0 plot

x_is_df = df_jet_breakup[df_jet_breakup['x_i/d_0'].notnull()]
x_is_df = x_is_df[x_is_df['regime turb'] == 'turbulent']
x_is_df = x_is_df[x_is_df['Re_l0'] > 5000]
x_is_df = x_is_df[x_is_df['Ma_g'] < 0.4]
#x_is_df = x_is_df[x_is_df['x_i/d_0'] > 4] # TODO: Why did I look at higher <x_i>/d_0? I recall thinking of a reason why lower <x_i>/d_0 should plateau, and I think I wrote it down somewhere, but I can't recall it now (2018-02-21).
x_is_df = x_is_df[x_is_df['key'] != 'reitz_atomization_1978'] # Theirs is from steady photos that likely won't resolve this right. That's what it seems off.
#x_is_df = x_is_df[x_is_df['key'] != 'kim_investigation_1983']
parameter_space_plots(x_is_df, 'x_i/d_0', revno=revno)
latex_summary_table(x_is_df, 'x_i/d_0')

A = np.column_stack([np.ones(len(x_is_df)), log(x_is_df['We_l0'] * x_is_df['I_0']**3)])#, log(x_is_df['rho_s'])])
B = log(x_is_df['x_i/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
#a, C_We_l0I_03, C_rho_s = result
a, C_We_l0I_03 = result

C_x_i = exp(a)

print 'C_x_i       = ', C_x_i
print 'C_We_l0I_03 = ', C_We_l0I_03
#print 'C_Re_l0     = ', C_Re_l0
#print 'C_I_0       = ', C_I_0
#print 'C_rho_s     = ', C_rho_s

x_is_df['x_i/d_0 predicted'] = C_x_i * (x_is_df['We_l0'] * x_is_df['I_0']**3)**C_We_l0I_03# * x_is_df['rho_s']**C_rho_s
plot_with_keys(x_is_df, 'correlation', 'x_i/d_0 predicted', 'x_i/d_0', plot_type='loglog', add_line=True, revno=revno)
print 'R^2 =', coeff_of_determination(x_is_df['x_i/d_0 predicted'], x_is_df['x_i/d_0'])
print '(log) R^2 =', coeff_of_determination(log(x_is_df['x_i/d_0 predicted']), log(x_is_df['x_i/d_0']))
print

macros_reg.write(r'\newcommand{\xiavgreg}{\frac{\xiavg}{d_0} = '+roundstr(C_x_i)+r' \left(\Tubarexp{3} \Welo\right)^{'+roundstr(C_We_l0I_03)+'}}\n')
macros_reg.write(r'\newcommand{\xiavgregrsquared}{'+roundstr(coeff_of_determination(x_is_df['x_i/d_0 predicted'], x_is_df['x_i/d_0']))+'}\n')
macros_reg.write(r'\newcommand{\xiavglogregrsquared}{'+roundstr(coeff_of_determination(log(x_is_df['x_i/d_0 predicted']), log(x_is_df['x_i/d_0'])))+'}\n')
macros_reg.write(r'\newcommand{\xiavgregN}{\num{'+str(len(x_is_df))+'}}\n')
macros_reg.write('\n')

C_We_l0I_03 = -1.
We_l0I_03 = (x_is_df['We_l0'] * x_is_df['I_0']**3.)**C_We_l0I_03
C_x_i = We_l0I_03.dot(x_is_df['x_i/d_0']) / We_l0I_03.dot(We_l0I_03)
print 'C_x_i       = ', C_x_i
print 'C_We_l0I_03 = ', C_We_l0I_03
x_is_df['x_i/d_0 predicted'] = C_x_i * (x_is_df['We_l0'] * x_is_df['I_0']**3.)**C_We_l0I_03
print '(alt) R^2 =', coeff_of_determination(x_is_df['x_i/d_0 predicted'], x_is_df['x_i/d_0'])
print '(alt log) R^2 =', coeff_of_determination(log(x_is_df['x_i/d_0 predicted']), log(x_is_df['x_i/d_0']))

macros_reg.write(r'\newcommand{\xiavgaltreg}{\frac{\xiavg}{d_0} = '+roundstr(C_x_i)+r' \left(\Tubarexp{3} \Welo\right)^{-1}}'+'\n')
macros_reg.write(r'\newcommand{\xiavgaltregrsquared}{'+roundstr(coeff_of_determination(x_is_df['x_i/d_0 predicted'], x_is_df['x_i/d_0']))+'}\n')
macros_reg.write(r'\newcommand{\xiavgaltlogregrsquared}{'+roundstr(coeff_of_determination(log(x_is_df['x_i/d_0 predicted']), log(x_is_df['x_i/d_0'])))+'}\n')
macros_reg.write('\n')

#x_is_df2 = x_is_df
#x_is_df  = x_is_df[x_is_df['x_i/d_0'] > 0.3]

#A = np.column_stack([np.zeros(len(x_is_df)), 1. / (x_is_df['We_l0'] * x_is_df['I_0']**3)])
#B = x_is_df['x_i/d_0']

#result, _, _, _ = np.linalg.lstsq(A, B)
#C_delta, C_We_l0I_03 = result

#print 'C_delta     = ', C_delta
#print 'C_We_l0I_03 = ', C_We_l0I_03

#x_is_df['x_i/d_0 predicted'] = C_delta + C_We_l0I_03 / (x_is_df['We_l0'] * x_is_df['I_0']**3)
#plot_with_keys(x_is_df, 'correlation', 'x_i/d_0 predicted', 'x_i/d_0', plot_type='loglog', add_line=True, revno=revno)
#print 'R^2 =', coeff_of_determination(x_is_df['x_i/d_0 predicted'], x_is_df['x_i/d_0'])
#print

########################
# droplet size (D/d_0) #
########################

# # # # # # # # # # # # # # # # # # # # ##
# Comparison against data not from pipes #
## # # # # # # # # # # # # # # # # # # # #

# ruff_structure_1991, dumouchel_role_2005, balewski_experimental_2008

# zajac_correlation_1970's data is for jets impinging on solid surfaces, unfortunately.

# # # # # # # # # # # # # # # # # # ##
# Comparison against data from pipes #
## # # # # # # # # # # # # # # # # # #

D_32_is_df = df_jet_breakup[df_jet_breakup['D_32/d_0'].notnull()]
D_32_is_df = D_32_is_df[D_32_is_df['regime turb'] == 'turbulent']
D_32_is_df = D_32_is_df[D_32_is_df['Re_l0'] > 5000]
D_32_is_df = D_32_is_df[D_32_is_df['Ma_g'] < 0.4]
D_32_is_df = D_32_is_df[D_32_is_df['rho_s'] > 500]
parameter_space_plots(D_32_is_df, 'D_32/d_0', revno=revno)
latex_summary_table(D_32_is_df, 'D_32/d_0')

A = np.column_stack([np.ones(len(D_32_is_df)), log(D_32_is_df['We_l0'] * D_32_is_df['I_0']**2.)])
B = log(D_32_is_df['D_32/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0I_02 = result

C_D_32 = exp(a)

print 'C_D_32      = ', C_D_32
print 'C_We_l0I_02 = ', C_We_l0I_02

D_32_is_df['D_32/d_0 predicted'] = C_D_32 * (D_32_is_df['We_l0'] * D_32_is_df['I_0']**2.)**C_We_l0I_02
plot_with_keys(D_32_is_df, 'correlation', 'D_32/d_0 predicted', 'D_32/d_0', plot_type='loglog', add_line=True, revno=revno)
print 'R^2 =', coeff_of_determination(D_32_is_df['D_32/d_0 predicted'], D_32_is_df['D_32/d_0'])
print

macros_reg.write(r'\newcommand{\SMDreg}{\frac{D_{32}}{d_0} = '+roundstr(C_D_32)+r' \left(\Tubarexp{2} \Welo\right)^{'+roundstr(C_We_l0I_02)+'}}\n')
macros_reg.write(r'\newcommand{\SMDregrsquared}{'+roundstr(coeff_of_determination(D_32_is_df['D_32/d_0 predicted'], D_32_is_df['D_32/d_0']))+'}\n')
macros_reg.write(r'\newcommand{\SMDregN}{\num{'+str(len(D_32_is_df))+'}}\n')
macros_reg.write('\n')

# set the exponent manually to get C_{v_\text{d}} consistent with theory

C_We_l0I_02 = -3./5.
We_l0I_02 = (D_32_is_df['We_l0'] * D_32_is_df['I_0']**2.)**C_We_l0I_02
C_D_32 = We_l0I_02.dot(D_32_is_df['D_32/d_0']) / We_l0I_02.dot(We_l0I_02)

print 'C_D_32      = ', C_D_32
print 'C_We_l0I_02 = ', C_We_l0I_02

D_32_is_df['D_32/d_0 predicted'] = C_D_32 * (D_32_is_df['We_l0'] * D_32_is_df['I_0']**2.)**C_We_l0I_02
print '(alt) R^2 =', coeff_of_determination(D_32_is_df['D_32/d_0 predicted'], D_32_is_df['D_32/d_0'])
print

#################################
# droplet velocity (v_d_bar/vp) #
#################################

v_d_bar_is = df_jet_breakup[df_jet_breakup['v_d_bar/vp'].notnull()]
v_d_bar_is = v_d_bar_is[v_d_bar_is['regime turb'] == 'turbulent']
v_d_bar_is = v_d_bar_is[v_d_bar_is['Re_l0'] > 5000]
v_d_bar_is = v_d_bar_is[v_d_bar_is['Ma_g'] < 0.4]
v_d_bar_is = v_d_bar_is[v_d_bar_is['rho_s'] > 500]
parameter_space_plots(v_d_bar_is, 'v_d_bar/vp', revno=revno)
latex_summary_table(v_d_bar_is, 'v_d_bar/vp')

A = np.column_stack([np.ones(len(v_d_bar_is)), log(v_d_bar_is['We_l0'] * v_d_bar_is['I_0']**2.)])
B = log(v_d_bar_is['v_d_bar/vp'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0I_02 = result

C_v_d_bar = exp(a)

print 'C_v_d_bar   = ', C_v_d_bar
print 'C_We_l0I_02 = ', C_We_l0I_02

v_d_bar_is['v_d_bar/vp predicted'] = C_v_d_bar * (v_d_bar_is['We_l0'] * v_d_bar_is['I_0']**2.)**C_We_l0I_02
plot_with_keys(v_d_bar_is, 'correlation', 'v_d_bar/vp predicted', 'v_d_bar/vp', plot_type='linear', add_line=True, revno=revno)
print 'R^2 =', coeff_of_determination(v_d_bar_is['v_d_bar/vp predicted'], v_d_bar_is['v_d_bar/vp'])
print

macros_reg.write(r'\newcommand{\vdavgreg}{\frac{\vdavg}{\vprimebar{0}} = '+roundstr(C_v_d_bar)+r' \left(\Tubarexp{2} \Welo\right)^{'+roundstr(C_We_l0I_02)+'}}\n')
macros_reg.write(r'\newcommand{\vdavgregrsquared}{'+roundstr(coeff_of_determination(v_d_bar_is['v_d_bar/vp predicted'], v_d_bar_is['v_d_bar/vp']))+'}\n')
macros_reg.write(r'\newcommand{\vdbarregN}{\num{'+str(len(v_d_bar_is))+'}}\n')
macros_reg.write('\n')

# set the exponent manually to get C_{v_\text{d}} consistent with theory

C_We_l0I_02 = -1./5.
We_l0I_02 = (v_d_bar_is['We_l0'] * v_d_bar_is['I_0']**2.)**C_We_l0I_02
C_v_d = We_l0I_02.dot(v_d_bar_is['v_d_bar/vp']) / We_l0I_02.dot(We_l0I_02)

print 'C_v_d       = ', C_v_d
print 'C_We_l0I_02 = ', C_We_l0I_02

v_d_bar_is['v_d_bar/vp predicted'] = C_v_d_bar * (v_d_bar_is['We_l0'] * v_d_bar_is['I_0']**2.)**C_We_l0I_02
#plot_with_keys(v_d_bar_is, 'correlation', 'v_d_bar/vp predicted', 'v_d_bar/vp', plot_type='loglog', add_line=True, revno=revno)
#print '(alt) R^2 =', coeff_of_determination(v_d_bar_is['v_d_bar/vp predicted'], v_d_bar_is['v_d_bar/vp']) # TODO: Figure out why this fails.
print

# TODO: Fix this hack that gets citations working in the legends.
os.system("cd ../outputs/figures/ ; for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")

macros_reg.close()
