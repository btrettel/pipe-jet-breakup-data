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

macros_model_comparisons = open('../outputs/macros/model_comparisons.tex', 'w')

print

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

#####################################
# Sauter mean diameter (D_{32}/d_0) #
#####################################

print 'Sauter mean diameter:'
print

D_32_is_df = df_jet_breakup[df_jet_breakup['D_32/d_0'].notnull()]
D_32_is_df = D_32_is_df[D_32_is_df['regime turb'] == 'turbulent']
D_32_is_df = D_32_is_df[D_32_is_df['Re_l0'] > 5000]
D_32_is_df = D_32_is_df[D_32_is_df['Ma_g'] < 0.4]

#Huh_func = np.ones(len(D_32_is_df['D_32/d_0']))
Huh_func = D_32_is_df['D_32/d_0'] / D_32_is_df['D_32/d_0']
C_D_32_Huh = Huh_func.dot(D_32_is_df['D_32/d_0']) / Huh_func.dot(Huh_func)
D_32_is_predicted = C_D_32_Huh * Huh_func
print 'C_D_32_Huh =', C_D_32_Huh
D_32_Huh_R2 = coeff_of_determination(D_32_is_predicted, D_32_is_df['D_32/d_0'])
print 'R^2 =', D_32_Huh_R2
print

Oh_l0 = np.sqrt(D_32_is_df['We_l0']) / D_32_is_df['Re_l0']
lambda_ms = (18.04 * (1. + 0.54 * np.sqrt(Oh_l0)) * (1. + 0.4 * ((D_32_is_df['We_l0'] / D_32_is_df['Re_l0']) / np.sqrt(D_32_is_df['rho_s']))**0.7)) / (1. + 0.27 * (D_32_is_df['We_l0'] / D_32_is_df['rho_s'])**1.67)**0.6
KHRT_func = lambda_ms / 2.
C_D_32_KHRT = KHRT_func.dot(D_32_is_df['D_32/d_0']) / KHRT_func.dot(KHRT_func)
D_32_is_predicted = C_D_32_KHRT * KHRT_func
print 'C_D_32_KHRT =', C_D_32_KHRT
D_32_KHRT_R2 = coeff_of_determination(D_32_is_predicted, D_32_is_df['D_32/d_0'])
print 'R^2 =', D_32_KHRT_R2
print

CDRSV_func = D_32_is_df['I_0']**(-6./5.)*D_32_is_df['We_l0']**(-3./5.)
C_D_32_CDRSV = CDRSV_func.dot(D_32_is_df['D_32/d_0']) / CDRSV_func.dot(CDRSV_func)
D_32_is_predicted = C_D_32_CDRSV * CDRSV_func
print 'C_D_32_CDRSV =', C_D_32_CDRSV, '(same for Faeth)'
D_32_CDRSV_R2 = coeff_of_determination(D_32_is_predicted, D_32_is_df['D_32/d_0'])
print 'R^2 =', D_32_CDRSV_R2
print

macros_model_comparisons.write(r'\newcommand{\CSMD}{'+roundstr(C_D_32_CDRSV)+'}\n')
macros_model_comparisons.write(r'\newcommand{\CSMDRsquared}{'+roundstr(D_32_CDRSV_R2)+'}\n\n')

D_32_is_df['D_32/d_0 predicted'] = D_32_is_predicted
plot_with_keys(D_32_is_df, 'CDRSV', 'D_32/d_0 predicted', 'D_32/d_0', plot_type='loglog', add_line=True, revno=revno)

with open('../outputs/data/TSB_D_32_is.pickle') as f:
   C_D_32, C_We_l0I_02 = pickle.load(f)
D_32_is_predicted = C_D_32*(D_32_is_df['I_0']**(2.)*D_32_is_df['We_l0'])**(C_We_l0I_02)
D_32_reg_R2 = coeff_of_determination(D_32_is_predicted, D_32_is_df['D_32/d_0'])
print 'regression:'
print 'R^2 =', D_32_reg_R2
print
print

#########################################
# droplet velocity ($\vdavg/\vprime_0$) #
#########################################

print 'Droplet velocity:'
print

v_d_bar_is_df = df_jet_breakup[df_jet_breakup['v_d_bar/vp'].notnull()]
v_d_bar_is_df = v_d_bar_is_df[v_d_bar_is_df['regime turb'] == 'turbulent']
v_d_bar_is_df = v_d_bar_is_df[v_d_bar_is_df['Re_l0'] > 5000]
v_d_bar_is_df = v_d_bar_is_df[v_d_bar_is_df['Ma_g'] < 0.4]

CDRSV_func = v_d_bar_is_df['I_0']**(-2./5.)*v_d_bar_is_df['We_l0']**(-1./5.)
C_v_d_bar_CDRSV = CDRSV_func.dot(v_d_bar_is_df['v_d_bar/vp']) / CDRSV_func.dot(CDRSV_func)
v_d_bar_is_predicted = C_v_d_bar_CDRSV * CDRSV_func
print 'C_v_d_bar_CDRSV =', C_v_d_bar_CDRSV, '(same for Faeth)'
v_d_bar_CDRSV_R2 = coeff_of_determination(v_d_bar_is_predicted, v_d_bar_is_df['v_d_bar/vp'])
print 'R^2 =', v_d_bar_CDRSV_R2
print

macros_model_comparisons.write(r'\newcommand{\Cvdavg}{'+roundstr(C_v_d_bar_CDRSV)+'}\n')
macros_model_comparisons.write(r'\newcommand{\CvdavgRsquared}{'+roundstr(v_d_bar_CDRSV_R2)+'}\n\n')

#######################################
# breakup onset location (\xiavg/d_0) #
#######################################

print 'Breakup onset length:'
print

xiavgs_df = df_jet_breakup[df_jet_breakup['x_i/d_0'].notnull()]
xiavgs_df = xiavgs_df[xiavgs_df['regime turb'] == 'turbulent']
xiavgs_df = xiavgs_df[xiavgs_df['Re_l0'] > 5000]
xiavgs_df = xiavgs_df[xiavgs_df['Ma_g'] < 0.4]
xiavgs_df = xiavgs_df[xiavgs_df['key'] != 'reitz_atomization_1978'] # Theirs is from steady photos that likely won't resolve this right. That's what it seems off.
#xiavgs_df = xiavgs_df[xiavgs_df['key'] != 'kim_investigation_1983']

Faeth_func = xiavgs_df['I_0']**(-9./5.) * xiavgs_df['We_l0']**(-2./5.)
C_xiavg_Faeth = Faeth_func.dot(xiavgs_df['x_i/d_0']) / Faeth_func.dot(Faeth_func)
xiavgs_predicted = C_xiavg_Faeth * Faeth_func
print 'C_xiavg_Faeth =', C_xiavg_Faeth
xiavgs_Faeth_R2 = coeff_of_determination(xiavgs_predicted, xiavgs_df['x_i/d_0'])
print 'R^2 =', xiavgs_Faeth_R2
print

CDRSV_func = xiavgs_df['I_0']**(-3.)*xiavgs_df['We_l0']**(-1.)
C_xiavg_CDRSV = CDRSV_func.dot(xiavgs_df['x_i/d_0']) / CDRSV_func.dot(CDRSV_func)
xiavgs_predicted = C_xiavg_CDRSV * CDRSV_func
print 'C_xiavg_CDRSV =', C_xiavg_CDRSV
xiavgs_CDRSV_R2 = coeff_of_determination(xiavgs_predicted, xiavgs_df['x_i/d_0'])
print 'R^2 =', xiavgs_CDRSV_R2
print

macros_model_comparisons.write(r'\newcommand{\Cxiavg}{'+roundstr(C_xiavg_CDRSV)+'}\n')
macros_model_comparisons.write(r'\newcommand{\CxiavgRsquared}{'+roundstr(xiavgs_CDRSV_R2)+'}\n\n')

xiavgs_df['x_i/d_0 predicted'] = xiavgs_predicted
plot_with_keys(xiavgs_df, 'CDRSV', 'x_i/d_0 predicted', 'x_i/d_0', plot_type='loglog', add_line=True, revno=revno)

with open('../outputs/data/TSB_xiavg.pickle') as f:
   C_x_i, C_We_l0I_03 = pickle.load(f)
xiavgs_predicted = C_x_i*(xiavgs_df['I_0']**(3.)*xiavgs_df['We_l0'])**(C_We_l0I_03)
xiavgs_reg_R2 = coeff_of_determination(xiavgs_predicted, xiavgs_df['x_i/d_0'])
print 'regression:'
print 'R^2 =', xiavgs_reg_R2
print
print

###############################
# breakup length (\xbavg/d_0) #
###############################

# TODO: Computer R^2 for full data set including photographic data.

print 'Breakup length:'
print

xbavgs_df = df_jet_breakup[df_jet_breakup['L_b/d_0'].notnull()]
xbavgs_df = xbavgs_df[xbavgs_df['regime turb'] == 'turbulent']
xbavgs_df = xbavgs_df[xbavgs_df['Ma_g'] < 0.3]
xbavgs_df = xbavgs_df[xbavgs_df['regime L_b'] == 'second wind-induced']
xbavgs_df_more = xbavgs_df
xbavgs_df = xbavgs_df[xbavgs_df['L_b method'] == 'electrical']

Faeth_func = xbavgs_df['We_l0']**(1./2.)
C_xbavg_Faeth = Faeth_func.dot(xbavgs_df['L_b/d_0']) / Faeth_func.dot(Faeth_func)
Faeth_func = xbavgs_df_more['We_l0']**(1./2.)
xbavgs_predicted = C_xbavg_Faeth * Faeth_func
print 'C_xbavg_Faeth =', C_xbavg_Faeth
xbavg_Faeth_R2 = coeff_of_determination(xbavgs_predicted, xbavgs_df_more['L_b/d_0'])
print 'R^2 =', xbavg_Faeth_R2
print

KHRT_func = 0.5 * xbavgs_df['rho_s']**(1./2.)
C_xbavg_KHRT = KHRT_func.dot(xbavgs_df['L_b/d_0']) / KHRT_func.dot(KHRT_func)
KHRT_func = 0.5 * xbavgs_df_more['rho_s']**(1./2.)
xbavgs_predicted = C_xbavg_KHRT * KHRT_func
print 'C_xbavg_KHRT =', C_xbavg_KHRT
xbavg_KHRT_R2 = coeff_of_determination(xbavgs_predicted, xbavgs_df_more['L_b/d_0'])
print 'R^2 =', xbavg_KHRT_R2
print

CDRSV_func = xbavgs_df['I_0']**(-3./5.)*xbavgs_df['We_l0']**(1./5.)
C_xbavg_CDRSV = CDRSV_func.dot(xbavgs_df['L_b/d_0']) / CDRSV_func.dot(CDRSV_func)
CDRSV_func = xbavgs_df_more['I_0']**(-3./5.)*xbavgs_df_more['We_l0']**(1./5.)
xbavgs_predicted = C_xbavg_CDRSV * CDRSV_func
print 'C_xbavg_CDRSV =', C_xbavg_CDRSV
xbavg_CDRSV_R2 = coeff_of_determination(xbavgs_predicted, xbavgs_df_more['L_b/d_0'])
print 'R^2 =', xbavg_CDRSV_R2
print

macros_model_comparisons.write(r'\newcommand{\Cxbavg}{'+roundstr(C_xbavg_CDRSV)+'}\n')
macros_model_comparisons.write(r'\newcommand{\CxbavgRsquared}{'+roundstr(xbavg_CDRSV_R2)+'}\n\n')

xbavgs_df_more['L_b/d_0 predicted'] = xbavgs_predicted
plot_with_keys(xbavgs_df_more, 'CDRSV', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno)

with open('../outputs/data/TSB_xbavg.pickle') as f:
   C_TSB, alpha_We, alpha_Re_l0_2WI, alpha_Tu_2WI, alpha_rho_s_2WI = pickle.load(f)
xbavgs_predicted = C_TSB*xbavgs_df_more['I_0']**alpha_Tu_2WI*xbavgs_df_more['We_l0']**alpha_We
xbavg_reg_R2 = coeff_of_determination(xbavgs_predicted, xbavgs_df_more['L_b/d_0'])
print 'regression:'
print 'R^2 =', xbavg_reg_R2
print
print

#########################
# spray angle (\thetai) #
#########################

thetai_df = df_jet_breakup[df_jet_breakup['theta'].notnull()]
thetai_df = thetai_df[thetai_df['regime turb'] == 'turbulent']
thetai_df = thetai_df[thetai_df['Ma_g'] < 0.3]
thetai_df = thetai_df[thetai_df['key'] != 'reitz_atomization_1978']
thetai_df = thetai_df[thetai_df['key'] != 'arai_break-up_1985']
thetai_df = thetai_df[thetai_df['key'] != 'rupe_dynamic_1962']
thetai_df = thetai_df[thetai_df['key'] != 'grant_newtonian_1965']
thetai_df = thetai_df[thetai_df['key'] != 'hoyt_pipe-exit_1980']
thetai_df = thetai_df[thetai_df['key'] != 'hoyt_effect_1985']
thetai_df = thetai_df[thetai_df['key'] != 'wu_liquid_1992']
thetai_df = thetai_df[thetai_df['key'] != 'skrebkov_turbulentnyye_1963']
thetai_df['We_l0_crit'] = We_l0_crit(thetai_df['I_0'], thetai_df['rho_s'])
thetai_df = thetai_df[thetai_df['We_l0'] < thetai_df['We_l0_crit']] # Exclude atomization points.

print 'Spray angle:'
print

Huh_func = thetai_df['I_0']
C_thetai_Huh = Huh_func.dot(np.tan(thetai_df['theta'] / 2.)) / Huh_func.dot(Huh_func)
thetai_predicted = 2. * np.arctan(C_thetai_Huh * Huh_func)
print 'C_thetai_Huh =', C_thetai_Huh
thetai_Huh_R2 = coeff_of_determination(thetai_predicted, thetai_df['theta'])
print 'R^2 =', thetai_Huh_R2
print

Oh_l0 = np.sqrt(thetai_df['We_l0']) / thetai_df['Re_l0']
lambda_ms = (18.04 * (1. + 0.54 * np.sqrt(Oh_l0)) * (1. + 0.4 * ((thetai_df['We_l0'] / thetai_df['Re_l0']) / np.sqrt(thetai_df['rho_s']))**0.7)) / (1. + 0.27 * (thetai_df['We_l0'] / thetai_df['rho_s'])**1.67)**0.6
omega_ms  = (0.96 + 0.38 * (thetai_df['We_l0'] / thetai_df['rho_s'])**1.5) / ((1. + np.sqrt(2.) * Oh_l0) * (1. + 1.4 * ((thetai_df['We_l0'] / thetai_df['Re_l0']) / thetai_df['rho_s'])**0.6))
KHRT_func = lambda_ms * omega_ms / thetai_df['We_l0']
C_thetai_KHRT = KHRT_func.dot(np.tan(thetai_df['theta'] / 2.)) / KHRT_func.dot(KHRT_func)
thetai_predicted = 2. * np.arctan(C_thetai_KHRT * KHRT_func)
print 'C_thetai_KHRT =', C_thetai_KHRT
thetai_KHRT_R2 = coeff_of_determination(thetai_predicted, thetai_df['theta'])
print 'R^2 =', thetai_KHRT_R2
print

CDRSV_func = thetai_df['I_0']**(3./5.)*thetai_df['We_l0']**(-1./5.)
C_thetai_CDRSV = CDRSV_func.dot(np.tan(thetai_df['theta'] / 2.)) / CDRSV_func.dot(CDRSV_func)
thetai_predicted = 2. * np.arctan(C_thetai_CDRSV * CDRSV_func)
print 'C_thetai_CDRSV =', C_thetai_CDRSV
thetai_CDRSV_R2 = coeff_of_determination(thetai_predicted, thetai_df['theta'])
print 'R^2 =', thetai_CDRSV_R2
print

macros_model_comparisons.write(r'\newcommand{\Cthetai}{'+roundstr(C_thetai_CDRSV)+'}\n')
macros_model_comparisons.write(r'\newcommand{\CthetaiRsquared}{'+roundstr(thetai_CDRSV_R2)+'}\n')

with open('../outputs/data/TSB_thetai.pickle') as f:
   C_theta, alpha_We_l0_theta, alpha_Tubar_0_theta = pickle.load(f)
thetai_predicted = 2. * np.arctan(C_theta*thetai_df['I_0']**alpha_Tubar_0_theta*thetai_df['We_l0']**alpha_We_l0_theta)
thetai_reg_R2 = coeff_of_determination(thetai_predicted, thetai_df['theta'])
print 'regression:'
print 'R^2 =', thetai_reg_R2
print
print

f = open('../outputs/tables/model_coefficient_and_r2_table.tex', 'w')
f.write(r'\begin{table}'+'\n')
f.write(r'\centering'+'\n')
f.write(r'\begin{tabular}{r|cccc|cccc}'+'\n')
f.write(r' & \multicolumn{4}{c}{coefficients} & \multicolumn{4}{|c}{$R^2$} \\'+'\n')
f.write(r' & $D_{ij}$ & $\xiavg$ & $\xbavg$ & $\thetai$ & $D_{ij}$ & $\xiavg$ & $\xbavg$ & $\thetai$ \\'+'\n')
f.write(r'\hline'+'\n')
f.write(r'Faeth & '+roundstr(C_D_32_CDRSV)+' & '+roundstr(C_xiavg_Faeth)+' & '+roundstr(C_xbavg_Faeth)+' & '+'---'+r' & '+roundstr(D_32_CDRSV_R2)+' & '+roundstr(xiavgs_Faeth_R2)+' & '+roundstr(xbavg_Faeth_R2)+' & '+'---'+r' \\'+'\n')
f.write(r'Huh & '+roundstr(C_D_32_Huh)+' & '+'---'+' & '+'---'+' & '+roundstr(C_thetai_Huh)+r' & '+roundstr(D_32_Huh_R2)+' & '+'---'+' & '+'---'+' & '+roundstr(thetai_Huh_R2)+r' \\'+'\n')
f.write(r'KH-RT & '+roundstr(C_D_32_KHRT)+' & '+'---'+' & '+roundstr(C_xbavg_KHRT)+' & '+roundstr(C_thetai_KHRT)+r' & '+roundstr(D_32_KHRT_R2)+' & '+'---'' & '+roundstr(xbavg_KHRT_R2)+' & '+roundstr(thetai_KHRT_R2)+r' \\'+'\n')
f.write(r'CDRSV & '+roundstr(C_D_32_CDRSV)+' & '+roundstr(C_xiavg_CDRSV)+' & '+roundstr(C_xbavg_CDRSV)+' & '+roundstr(C_thetai_CDRSV)+r' & '+roundstr(D_32_CDRSV_R2)+' & '+roundstr(xiavgs_CDRSV_R2)+' & '+roundstr(xbavg_CDRSV_R2)+' & '+roundstr(thetai_CDRSV_R2)+r' \\'+'\n')
f.write(r'regression & --- & --- & --- & --- & '+roundstr(D_32_reg_R2)+' & '+roundstr(xiavgs_reg_R2)+' & '+roundstr(xbavg_reg_R2)+' & '+roundstr(thetai_reg_R2)+'\n')
f.write(r'\end{tabular}'+'\n')
f.write(r'\caption{Calibrated model coefficients for multiple models and associated coefficients of determination ($R^2$).}'+'\n')
f.write(r'\label{tab:model-coefficient-table}'+'\n')
f.write(r'\end{table}'+'\n')
f.close()

# close macros file
macros_model_comparisons.close()

# TODO: Fix this hack that gets citations working in the legends.
os.system("cd ../outputs/figures/ ; for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")
