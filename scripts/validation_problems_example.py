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
import scipy.stats

#revno = lastchangedrevision[0:6]
revno = None

#def roundstr(number):
   #if abs(number) < 1e-8:
      ##return roundstre(number)
      #return r'\num{0}'
   #else:
      #number = "%.3f" % round(number, 3)
      #return r'\num{'+number+r'}'

def roundstrzero(number):
   if number < 0.:
      return '$<0$'
   else:
      number = "%.3f" % round(number, 3)
      return r'\num{'+number+r'}'

#def roundstre(number):
   #number = '%.2E' % Decimal(number)
   #return r'\num{'+number+r'}'

with open('../outputs/data/'+data_file+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

# Full case
print '--------------'
print 'Full data case'
print '--------------'

#L_bs_df = df_jet_breakup[df_jet_breakup['theta'].notnull()]
#L_bs_df = L_bs_df[L_bs_df['regime turb'] == 'turbulent']
#L_bs_df = L_bs_df[L_bs_df['Re_l0'] > 5000.]
#L_bs_df = L_bs_df[L_bs_df['Ma_g'] < 0.4]
#L_bs_df = L_bs_df[L_bs_df['rho_s'] > 500.]
#summary_table(L_bs_df, 'theta')

L_bs_df = df_jet_breakup[df_jet_breakup['L_b/d_0'].notnull()]
L_bs_df = L_bs_df[L_bs_df['regime turb'] == 'turbulent']
#L_bs_df = L_bs_df[L_bs_df['rho_s'] > 500]
L_bs_df = L_bs_df[L_bs_df['Ma_g'] < 0.3]
L_bs_df = L_bs_df[L_bs_df['regime L_b'] == 'second wind-induced']
#L_bs_df = L_bs_df[L_bs_df['key'] != 'phinney_breakup_1973'] # Removed until I can understand why it seems wrong.
L_bs_df_more = L_bs_df
L_bs_df = L_bs_df[L_bs_df['L_b method'] == 'electrical']
summary_table(L_bs_df, 'L_b/d_0')

print 'We only:'
A = np.column_stack([np.ones(len(L_bs_df)), log(L_bs_df['We_l0'])])
B = log(L_bs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0 = result

C_b = exp(a)

print 'C_b =', C_b
print 'C_We_l0 =', C_We_l0

L_bs_df['L_b/d_0 predicted'] = C_b * L_bs_df['We_l0']**C_We_l0
Rsquared = coeff_of_determination(L_bs_df['L_b/d_0 predicted'], L_bs_df['L_b/d_0'])
parameter_space_plots(L_bs_df, 'L_b/d_0', revno=revno, filename_extra='_We')
plot_with_keys(L_bs_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_We')
print 'R^2 =', Rsquared
print

C_b_2      = C_b
C_Tu_2     = 0.
C_We_2     = C_We_l0
C_Re_2     = 0.
Rsquared_2 = Rsquared

# Re, We, and Tu
print 'Re, We, and Tu:'
A = np.column_stack([np.ones(len(L_bs_df)), log(L_bs_df['Re_l0']), log(L_bs_df['We_l0']), log(L_bs_df['I_0'])])
B = log(L_bs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_Re_l0, C_We_l0, C_I_0 = result

C_b = exp(a)

L_bs_df['L_b/d_0 predicted'] = C_b * L_bs_df['Re_l0']**C_Re_l0 * L_bs_df['We_l0']**C_We_l0 * L_bs_df['I_0']**C_I_0
Rsquared = coeff_of_determination(L_bs_df['L_b/d_0 predicted'], L_bs_df['L_b/d_0'])

delta_C_Re_l0 = scipy.stats.t.ppf(1 - (1 - interval_probability_level) / 2, len(L_bs_df) - 1) * (C_Re_l0 / sqrt(len(L_bs_df) - 2)) * sqrt(1 / Rsquared - 1)

print 'C_b =', C_b
print 'C_Re_l0 =', C_Re_l0, '+-', delta_C_Re_l0
print 'C_We_l0 =', C_We_l0
print 'C_I_0   =', C_I_0

macros_validation = open('../outputs/macros/validation.tex', 'w')
macros_validation.write(r'\newcommand{\CRelowithuncertainty}{'+roundstr(C_Re_l0)+r' \pm \num{'+str(round(10000. * delta_C_Re_l0) / 10000.)+'}}\n')

parameter_space_plots(L_bs_df, 'L_b/d_0', revno=revno, filename_extra='_all')
plot_with_keys(L_bs_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_all')
print 'R^2 =', Rsquared
print

C_b_1      = C_b
C_Tu_1     = C_I_0
C_We_1     = C_We_l0
C_Re_1     = C_Re_l0
Rsquared_1 = Rsquared

# We and Tu only
print 'We and Tu only:'
A = np.column_stack([np.ones(len(L_bs_df)), log(L_bs_df['We_l0']), log(L_bs_df['I_0'])])
B = log(L_bs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
#a, C_We_l0, C_Re_l0, C_I_0 = result
a, C_We_l0, C_I_0 = result

C_b = exp(a)

print 'C_b =', C_b
print 'C_We_l0 =', C_We_l0
#print 'C_Re_l0 =', C_Re_l0
print 'C_I_0   =', C_I_0

L_bs_df['L_b/d_0 predicted'] = C_b * L_bs_df['We_l0']**C_We_l0 * L_bs_df['I_0']**C_I_0
Rsquared = coeff_of_determination(L_bs_df['L_b/d_0 predicted'], L_bs_df['L_b/d_0'])
parameter_space_plots(L_bs_df, 'L_b/d_0', revno=revno, filename_extra='_WeTu')
plot_with_keys(L_bs_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_WeTu')
print 'R^2 =', Rsquared
print

C_b_3      = C_b
C_Tu_3     = C_I_0
C_We_3     = C_We_l0
C_Re_3     = 0.
Rsquared_3 = Rsquared

# Re only
print 'Re:'
A = np.column_stack([np.ones(len(L_bs_df)), log(L_bs_df['Re_l0'])])
B = log(L_bs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_Re_l0 = result

C_b = exp(a)

print 'C_b =', C_b
print 'C_Re_l0 =', C_Re_l0

L_bs_df['L_b/d_0 predicted'] = C_b * L_bs_df['Re_l0']**C_Re_l0
Rsquared = coeff_of_determination(L_bs_df['L_b/d_0 predicted'], L_bs_df['L_b/d_0'])
parameter_space_plots(L_bs_df, 'L_b/d_0', revno=revno, filename_extra='_Re')
plot_with_keys(L_bs_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_Re')
print 'R^2 =', Rsquared
print

C_b_4      = C_b
C_Tu_4     = 0.
C_We_4     = 0.
C_Re_4     = C_Re_l0
Rsquared_4 = Rsquared

print 'Writing LaTeX table file...'
f = open('../outputs/tables/example_table_all.tex', 'w')

f.write(r'\begin{table}'+'\n')
f.write(r'\centering'+'\n')
f.write(r'\begin{tabular}{r|ccccc}'+'\n')
f.write(r' & all & $\Welo$ & $\Welo$, $\Tubar_0$ & $\Relo$ \\'+'\n')
f.write(r'\hline'+'\n')
f.write(r'$C_\text{b}$ & '+roundstr(C_b_1)+' & '+roundstr(C_b_2)+' & '+roundstr(C_b_3)+' & '+roundstr(C_b_4)+r' \\'+'\n')
f.write(r'$C_\Tu$ & '+roundstr(C_Tu_1)+' & '+roundstr(C_Tu_2)+' & '+roundstr(C_Tu_3)+' & '+roundstr(C_Tu_4)+r' \\'+'\n')
f.write(r'$C_\We$ & '+roundstr(C_We_1)+' & '+roundstr(C_We_2)+' & '+roundstr(C_We_3)+' & '+roundstr(C_We_4)+r' \\'+'\n')
f.write(r'$C_\Re$ & '+roundstr(C_Re_1)+' & '+roundstr(C_Re_2)+' & '+roundstr(C_Re_3)+' & '+roundstr(C_Re_4)+r' \\'+'\n')
f.write(r'\hline'+'\n')
f.write(r'$R^2$ & '+roundstrzero(Rsquared_1)+' & '+roundstrzero(Rsquared_2)+' & '+roundstrzero(Rsquared_3)+' & '+roundstrzero(Rsquared_4)+'\n')
f.write(r'\end{tabular}'+'\n')
f.write(r'\caption{Table showing the effects of using different variables in a regression analysis for the breakup length, $\xbavg$. All the available electrical conductivity data fitting the data quality guidelines is used. The regression equation is the same as in \tabref{confounded-example-table}. Regression equation: $\xbavg/d_0 = C_\text{b} \Tubarexp{C_\Tu} \Welo^{C_\We} \Relo^{C_\Re}$. Total number of data points: '+str(len(L_bs_df_more))+'.\label{tab:full-example-table}}\n')
f.write(r'\end{table}'+'\n')

f.close()

macros_validation.write(r'\newcommand{\xbavgregrsquareda}{'+roundstr(Rsquared_1)+'}\n')
macros_validation.write(r'\newcommand{\xbavgregrsquaredb}{'+roundstr(Rsquared_2)+'}\n')
macros_validation.write(r'\newcommand{\xbavgregrsquaredc}{'+roundstr(Rsquared_3)+'}\n')
macros_validation.write(r'\newcommand{\xbavgregrsquaredd}{'+roundstr(Rsquared_4)+'}\n')
macros_validation.close()

# Confounded case
print '---------------'
print 'Confounded case'
print '---------------'

#L_bs_df_less = L_bs_df
L_bs_df_less = L_bs_df_more

#L_bs_df_less = L_bs_df_less[L_bs_df_less['key'] == 'chen_mechanics_1962']

#L_bs_df_less = L_bs_df_less[L_bs_df_less['key'] == 'arai_break-up_1985']

L_bs_df_less = L_bs_df_less[L_bs_df_less['key'] == 'kusui_liquid_1969']
L_bs_df_less = L_bs_df_less[L_bs_df_less['roughness'] == 'smooth']
#L_bs_df_less = L_bs_df_less[L_bs_df_less['Re_l0'] < 1.e5]

summary_table(L_bs_df_less, 'L_b/d_0')

print 'We only:'
A = np.column_stack([np.ones(len(L_bs_df_less)), log(L_bs_df_less['We_l0'])])
B = log(L_bs_df_less['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_We_l0 = result

C_b = exp(a)

print 'C_b =', C_b
print 'C_We_l0 =', C_We_l0

L_bs_df_less['L_b/d_0 predicted'] = C_b * L_bs_df_less['We_l0']**C_We_l0
Rsquared = coeff_of_determination(L_bs_df_less['L_b/d_0 predicted'], L_bs_df_less['L_b/d_0'])
parameter_space_plots(L_bs_df_less, 'L_b/d_0', revno=revno, filename_extra='_confounded_We')
plot_with_keys(L_bs_df_less, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_confounded_We')
print 'R^2 =', Rsquared
print

C_b_2      = C_b
C_Tu_2     = 0.
C_We_2     = C_We_l0
C_Re_2     = 0.
Rsquared_2 = Rsquared

# Re, We, and Tu
print 'Re, We, and Tu:'
A = np.column_stack([np.ones(len(L_bs_df_less)), log(L_bs_df_less['Re_l0']), log(L_bs_df_less['We_l0']), log(L_bs_df_less['I_0'])])
B = log(L_bs_df_less['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_Re_l0, C_We_l0, C_I_0 = result

C_b = exp(a)

L_bs_df_less['L_b/d_0 predicted'] = C_b * L_bs_df_less['Re_l0']**C_Re_l0 * L_bs_df_less['We_l0']**C_We_l0 * L_bs_df_less['I_0']**C_I_0
Rsquared = coeff_of_determination(L_bs_df_less['L_b/d_0 predicted'], L_bs_df_less['L_b/d_0'])

delta_C_Re_l0 = scipy.stats.t.ppf(1 - (1 - interval_probability_level) / 2, len(L_bs_df_less) - 1) * (C_Re_l0 / sqrt(len(L_bs_df_less) - 2)) * sqrt(1 / Rsquared - 1)

# Uncertainty estimates: Toolkit_10

print 'C_b =', C_b
print 'C_Re_l0 =', C_Re_l0, '+-', delta_C_Re_l0
print 'C_We_l0 =', C_We_l0
print 'C_I_0   =', C_I_0

parameter_space_plots(L_bs_df_less, 'L_b/d_0', revno=revno, filename_extra='_confounded_all')
plot_with_keys(L_bs_df_less, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_confounded_all')
print 'R^2 =', Rsquared
print

C_b_1      = C_b
C_Tu_1     = C_I_0
C_We_1     = C_We_l0
C_Re_1     = C_Re_l0
Rsquared_1 = Rsquared

# We and Tu only
print 'We and Tu only:'
A = np.column_stack([np.ones(len(L_bs_df_less)), log(L_bs_df_less['We_l0']), log(L_bs_df_less['I_0'])])
B = log(L_bs_df_less['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
#a, C_We_l0, C_Re_l0, C_I_0 = result
a, C_We_l0, C_I_0 = result

C_b = exp(a)

print 'C_b =', C_b
print 'C_We_l0 =', C_We_l0
#print 'C_Re_l0 =', C_Re_l0
print 'C_I_0   =', C_I_0

L_bs_df_less['L_b/d_0 predicted'] = C_b * L_bs_df_less['We_l0']**C_We_l0 * L_bs_df_less['I_0']**C_I_0
Rsquared = coeff_of_determination(L_bs_df_less['L_b/d_0 predicted'], L_bs_df_less['L_b/d_0'])
parameter_space_plots(L_bs_df_less, 'L_b/d_0', revno=revno, filename_extra='_confounded_WeTu')
plot_with_keys(L_bs_df_less, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_confounded_WeTu')
print 'R^2 =', Rsquared
print

C_b_3      = C_b
C_Tu_3     = C_I_0
C_We_3     = C_We_l0
C_Re_3     = 0.
Rsquared_3 = Rsquared

# Re only
print 'Re:'
A = np.column_stack([np.ones(len(L_bs_df_less)), log(L_bs_df_less['Re_l0'])])
B = log(L_bs_df_less['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, C_Re_l0 = result

C_b = exp(a)

print 'C_b =', C_b
print 'C_Re_l0 =', C_Re_l0

L_bs_df_less['L_b/d_0 predicted'] = C_b * L_bs_df_less['Re_l0']**C_Re_l0
Rsquared = coeff_of_determination(L_bs_df_less['L_b/d_0 predicted'], L_bs_df_less['L_b/d_0'])
parameter_space_plots(L_bs_df_less, 'L_b/d_0', revno=revno, filename_extra='_confounded_Re')
plot_with_keys(L_bs_df_less, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_confounded_Re')
print 'R^2 =', Rsquared
print

C_b_4      = C_b
C_Tu_4     = 0.
C_We_4     = 0.
C_Re_4     = C_Re_l0
Rsquared_4 = Rsquared

print 'Writing LaTeX table file...'
f = open('../outputs/tables/example_table_confounded.tex', 'w')

f.write(r'\begin{table}'+'\n')
f.write(r'\centering'+'\n')
f.write(r'\begin{tabular}{r|ccccc}'+'\n')
f.write(r' & all & $\Welo$ & $\Welo$, $\Tubar_0$ & $\Relo$ \\'+'\n')
f.write(r'\hline'+'\n')
f.write(r'$C_\text{b}$ & '+roundstr(C_b_1)+' & '+roundstr(C_b_2)+' & '+roundstr(C_b_3)+' & '+roundstr(C_b_4)+r' \\'+'\n')
f.write(r'$C_\Tu$ & '+roundstr(C_Tu_1)+' & '+roundstr(C_Tu_2)+' & '+roundstr(C_Tu_3)+' & '+roundstr(C_Tu_4)+r' \\'+'\n')
f.write(r'$C_\We$ & '+roundstr(C_We_1)+' & '+roundstr(C_We_2)+' & '+roundstr(C_We_3)+' & '+roundstr(C_We_4)+r' \\'+'\n')
f.write(r'$C_\Re$ & '+roundstr(C_Re_1)+' & '+roundstr(C_Re_2)+' & '+roundstr(C_Re_3)+' & '+roundstr(C_Re_4)+r' \\'+'\n')
f.write(r'\hline'+'\n')
f.write(r'$R^2$ & '+roundstrzero(Rsquared_1)+' & '+roundstrzero(Rsquared_2)+' & '+roundstrzero(Rsquared_3)+' & '+roundstrzero(Rsquared_4)+'\n')
f.write(r'\end{tabular}'+'\n')
f.write(r'\caption{Table showing the effects of using different variables in a regression analysis for the breakup length, $\xbavg$. Only data from \citet{kusui_liquid_1969} is used. Total number of data points: '+str(len(L_bs_df_less))+'.\label{tab:confounded-example-table}}\n')
f.write(r'\end{table}'+'\n')

# TODO: Add page number for arai_break-up_1985.

f.close()

# TODO: Fix this hack that gets citations working in the legends.
os.system("cd ../outputs/figures/ ; for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")
