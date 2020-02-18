#!/usr/bin/env python
# -*- coding: utf-8 -*-

from jetbreakup import *
from scipy.optimize import fsolve
#import scipy.stats

lastchangedby       = '$LastChangedBy: $'
lastchangedrevision = '$LastChangedRevision: $'
lastchangeddate     = '$LastChangedDate: $'
headurl             = '$HeadURL: $'

#revno = lastchangedrevision.replace('$LastChangedRevision: ', '').replace('$', '')
revno = None

macros_atomization = open('macros_atomization.tex', 'w')

input_filename = 'pipe-jet-breakup-data'

print

with open(input_filename+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

xbavgs_df = df_jet_breakup
xbavgs_df = xbavgs_df[xbavgs_df['regime L_b'] == 'atomization']
xbavgs_df = xbavgs_df[xbavgs_df['regime turb'] == 'turbulent']
xbavgs_df = xbavgs_df[xbavgs_df['Re_l0'] > Re_turb]

#summary_table(xbavgs_df)

## 2019-07-24:
## K_c, We_l0, Re_l0 do not seem to help.
## If there was cavitation, I suppose it would be at the beginning of the nozzle. As the pipe is long, I suspect the flow may be merely flashing, not cavitating, and the influence of strong fluctuations at the inlet are minimized. So perhaps K_c won't have an effect here? It could also be that none of these cases are cavitating. No one mentioned that they were. Though cavitation might explain why Arai's L_0/d_0 = 40 differs from L_0/d_0 = 50.
## Ma_g helps, so I am including it. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html

#gamma = 1.4 # all gases are diatomic (air or nitrogen) as of 2019-07-24

#A = np.column_stack([np.ones(len(xbavgs_df)), log(xbavgs_df['I_0']), log(xbavgs_df['rho_s']), log(1. + ((gamma - 1.) / 2.) * xbavgs_df['Ma_g'])])
#B = log(xbavgs_df['L_b/d_0'])

#result, _, _, _ = np.linalg.lstsq(A, B)
#a, alpha_Tu_A, alpha_rho_s, alpha_Ma_g = result

#C_A = exp(a)

#print 'C_A         =', C_A
#print 'alpha_Tu_A  =', alpha_Tu_A
#print 'alpha_rho_s =', alpha_rho_s
#print 'alpha_Ma_g  =', alpha_Ma_g

#xbavgs_df['L_b/d_0 predicted'] = C_A * xbavgs_df['I_0']**alpha_Tu_A * xbavgs_df['rho_s']**alpha_rho_s * (1. + ((gamma - 1.) / 2.) * xbavgs_df['Ma_g'])**alpha_Ma_g
#plot_with_keys(xbavgs_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_TA')
#print 'R^2 =', coeff_of_determination(xbavgs_df['L_b/d_0 predicted'], xbavgs_df['L_b/d_0'])

##print xbavgs_df['gas']
##print xbavgs_df['Ma_g']

#print
#print

# Low-Mach regression
xbavgs_df = xbavgs_df[xbavgs_df['Ma_g'] < 0.3] # 0.2 basically eliminates the variation in the density ratio
summary_table(xbavgs_df)

macros_atomization.write(r'\newcommand{\atomrhosrange}{\numrange{'+roundstr(np.amin(xbavgs_df['rho_s']), num=False)+'}{'+roundstr(np.amax(xbavgs_df['rho_s']), num=False)+'}}\n')

A = np.column_stack([np.ones(len(xbavgs_df)), log(xbavgs_df['I_0']), log(xbavgs_df['rho_s'])])
B = log(xbavgs_df['L_b/d_0'])

result, _, _, _ = np.linalg.lstsq(A, B)
a, alpha_Tu_A, alpha_rho_s = result

C_A = exp(a)

print 'C_A         =', C_A
print 'alpha_Tu_A  =', alpha_Tu_A
print 'alpha_rho_s =', alpha_rho_s

with open('atomization_xbavg.pickle', 'w') as f:
   pickle.dump([C_A, alpha_Tu_A, alpha_rho_s], f)

macros_atomization.write(r'\newcommand{\xbavgatomreg}{\frac{\xbavg}{d_0} = '+roundstr(C_A)+r' \Tubarexp{'+roundstr(alpha_Tu_A)+r'} \left(\frac{\rhol}{\rhog}\right)^{'+roundstr(alpha_rho_s)+'}}\n')

xbavgs_df['L_b/d_0 predicted'] = C_A * xbavgs_df['I_0']**alpha_Tu_A * xbavgs_df['rho_s']**alpha_rho_s
plot_with_keys(xbavgs_df, 'correlation', 'L_b/d_0 predicted', 'L_b/d_0', plot_type='linear', add_line=True, revno=revno, filename_extra='_atomization_low_Mach')
print 'R^2 =', coeff_of_determination(xbavgs_df['L_b/d_0 predicted'], xbavgs_df['L_b/d_0'])
macros_atomization.write(r'\newcommand{\xbavgatomregrsquared}{'+roundstr(coeff_of_determination(xbavgs_df['L_b/d_0 predicted'], xbavgs_df['L_b/d_0']))+'}\n')
macros_atomization.write(r'\newcommand{\xbavgatomregN}{\num{'+str(len(xbavgs_df))+'}}\n')

with open('TSB_xbavg.pickle') as f:
   C_TSB, alpha_We, alpha_Re_l0_2WI, alpha_Tu_2WI, alpha_rho_s_2WI = pickle.load(f)

C_TSBtoA           = (C_A / C_TSB)**(1./alpha_We)
alpha_Tu_2WItoA    = (alpha_Tu_A - alpha_Tu_2WI) / alpha_We
alpha_rho_s_2WItoA = alpha_rho_s / alpha_We

print
print 'C_TSBtoA           =', C_TSBtoA
print 'alpha_Tu_2WItoA    =', alpha_Tu_2WItoA
print 'alpha_rho_s_2WItoA =', alpha_rho_s_2WItoA
print

with open('atomization_boundary.pickle', 'w') as f:
   pickle.dump([C_TSBtoA, alpha_Tu_2WItoA, alpha_rho_s_2WItoA], f)

macros_atomization.write(r'\newcommand{\Wecritatom}{\Welocrit = '+roundstr(C_TSBtoA)+r' \Tubarexp{'+roundstr(alpha_Tu_2WItoA)+r'} \left(\frac{\rhol}{\rhog}\right)^{'+roundstr(alpha_rho_s_2WItoA)+'}}\n')
#macros_atomization.write(r'\newcommand{\Wegcritatom}{\Wegocrit \approx 2 \Tubarexp{-1}}'+'\n')
macros_atomization.write(r'\newcommand{\Wegcritatom}{\Wegocrit = '+roundstr(C_TSBtoA)+r' \Tubarexp{'+roundstr(alpha_Tu_2WItoA)+r'}}'+'\n')

We_g_crit = C_TSBtoA * (0.05**alpha_Tu_2WItoA) * ((1000./1.2)**(alpha_rho_s_2WItoA - 1.))
assert(We_g_crit < 100.)
assert(We_g_crit > 10.)
#We_g_crit = C_TSBtoA * (0.05**round(alpha_Tu_2WItoA, 3)) * ((1000./1.2)**(round(alpha_rho_s_2WItoA, 2) - 1.))
print 'We_g(Tu_0 = 5%) =', We_g_crit
macros_atomization.write(r'\newcommand{\WegcritReitz}{'+roundstr(We_g_crit)+'}\n')

# TODO: Analyze data from ahn_effects_2006 to determine K_\text{c,crit} for sharp-edged inlets.
# LATER: Try a step function in K_c for atomization breakup length. Different critical K_c values for rough and smooth inlets?
# LATER: Check Kusui's outlier in atomization breakup length data

macros_atomization.close()

# TODO: Fix this hack that gets citations working in the legends.
os.system("for i in *.pgf; do sed -i 's/TdTEi/\_/g' $i; done")
