#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import matplotlib as mpl
import os
import pandas as pd
import pickle
import re
import scipy
import sys
import socket
import warnings
import math

from decimal import Decimal
from numpy import sqrt, log, exp, std, mean
from math import erf, floor
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, DrawingArea, HPacker
from scipy import interpolate
from string import lower
from scipy.special import lambertw
from scipy.optimize import fsolve
from scipy.integrate import quad
from git import Repo

# Configuration
data_file = 'pipe-jet-breakup-data'

# determine root_dir
repo = Repo(str(__file__), search_parent_directories=True)
root_dir = repo.working_dir+'/'

#print 'Root directory:', root_dir

if not(os.path.exists(root_dir+'outputs/')):
   os.mkdir(root_dir+'outputs/')
if not(os.path.exists(root_dir+'outputs/data/')):
   os.mkdir(root_dir+'outputs/data/')
if not(os.path.exists(root_dir+'outputs/figures/')):
   os.mkdir(root_dir+'outputs/figures/')
if not(os.path.exists(root_dir+'outputs/macros/')):
   os.mkdir(root_dir+'outputs/macros/')
if not(os.path.exists(root_dir+'outputs/tables/')):
   os.mkdir(root_dir+'outputs/tables/')

headcommit = repo.head.commit

lastchangedby       = headcommit.committer.name
lastchangedrevision = headcommit.hexsha
lastchangeddate     = headcommit.authored_date

# https://tex.stackexchange.com/a/391078
mpl.use('pgf')
pgf_with_pdflatex_report = {
   "pgf.texsystem": "pdflatex",
   "text.usetex": True,         # use LaTeX to write all text
   "font.family": "serif",
   "font.serif": [],            # blank entries should cause plots 
   "font.sans-serif": [],       # to inherit fonts from the document
   "font.monospace": [],
   "font.size": 12,
   "xtick.direction": 'in',
   "ytick.direction": 'in',
   "grid.color": 'k',
   "grid.linewidth": 0.5,
   #"grid.linestyle": '--',
   #"grid.linestyle": 'dotted',
   "grid.linestyle": ':',
   "xtick.bottom": True,
   "xtick.top": True,
   "ytick.left": True,
   "ytick.right": True,
   "legend.framealpha": 1.,
   "legend.fancybox": False,
   "legend.shadow": False,
   "legend.edgecolor": 'k',
   "patch.linewidth": 0.5, # legend border thickness, https://stackoverflow.com/a/38379176/1124489
   "pgf.preamble": [
      r"\usepackage[T1]{fontenc}",
      r"\usepackage{amstext}",
      r"\usepackage{amsmath}",
      r"\usepackage{xcolor}",
      r"\renewcommand{\cite}[1]{#1}",
      r"\newcommand{\citet}[1]{#1}",
      r"\DeclareMathOperator{\LambertW}{W}",
      #r"\newcommand{\LambertW}{\mathrm{W}}",
      ],
   "legend.numpoints" : 1       # the number of points in the legend line, https://stackoverflow.com/a/6146871
}
pgf_with_pdflatex_talk = {
   "pgf.texsystem": "pdflatex",
   "text.usetex": True,         # use LaTeX to write all text
   "font.family": "sans-serif",
   "font.serif": [],            # blank entries should cause plots 
   "font.sans-serif": [],       # to inherit fonts from the document
   "font.monospace": [],
   "pgf.preamble": [
      r"\usepackage[T1]{fontenc}",
      r"\usepackage{amstext}",
      r"\usepackage{color}",
      r"\renewcommand{\cite}[1]{#1}",
      r"\newcommand{\citet}[1]{#1}",
      ],
}
pgf_with_pdflatex_talk = pgf_with_pdflatex_report
mpl.rcParams.update(pgf_with_pdflatex_report)
import matplotlib.pyplot as plt

water_df = pd.read_csv(root_dir+'data/fluid-properties/water.csv', sep=',', header=0)
water_df['Viscosity (m2/s)'] = water_df['Viscosity (Pa*s)'] / water_df['Density (kg/m3)']

g      = 9.80665 # m/s
MW_air = 28.97   # kg/kmol
MW_N2  = 2*14.007  # kg/kmol
P_atm  = 101325.  # Pa
Rbar  = 8.314e3  # J/(kmol*K)
T_zero = 273.15  # K

interval_probability_level = 0.95
Re_trans                   = 2300.
Re_turb                    = 4000.
T_std                      = 25.   # C, standard assumed gas temperature
f_ratio_cutoff             = 2.    # f / f_laminar of this or higher is classified as turbulent
N_pts                      = 1e3  # number of points to use when plotting smooth curves

# https://en.wikipedia.org/wiki/Viscosity#Effect_of_temperature_on_the_viscosity_of_a_gas
mu_0_air = 18.27e-6
T_0_air  = 291.15 # K
C_air    = 120 # K

# TODO: Add conversion between mm Hg and absolute atmospheric pressure

def rho_ideal_gas(pressure, temperature, MW): # kg/m^3
   # temperature in C
   return (pressure * MW) / (Rbar * (temperature + T_zero))

def mu_ideal_gas(temperature, mu_0, T_0, C): # Pa*s
   # temperature in C
   # Sutherland's formula
   # https://en.wikipedia.org/wiki/Viscosity#Effect_of_temperature_on_the_viscosity_of_a_gas
   return mu_0 * ((T_0 + C) / (temperature + T_zero + C)) * ((temperature + T_zero) / T_0)**(3/2)

def c_ideal_gas(gamma, T, MW): # m/s
   assert(T > 200.) # This is absolute temperature.
   return sqrt(gamma * Rbar * T / MW) 

def gamma_air(T):
   return np.interp(T, [-40., -20., 0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100., 200., 300., 400., 500., 1000.], [1.401, 1.401, 1.401, 1.401, 1.401, 1.401, 1.401, 1.401, 1.400, 1.400, 1.400, 1.399, 1.399, 1.399, 1.398, 1.397, 1.390, 1.379, 1.368, 1.357, 1.321])

rho_water = interpolate.interp1d(water_df['Temperature (C)'], water_df['Density (kg/m3)'])

mu_water = interpolate.interp1d(water_df['Temperature (C)'], water_df['Viscosity (Pa*s)'])

nu_water = interpolate.interp1d(water_df['Temperature (C)'], water_df['Viscosity (m2/s)'])

sigma_water = interpolate.interp1d(water_df['Temperature (C)'], water_df['Surf. Tension (N/m)'])

P_v_water = interpolate.interp1d(water_df['Temperature (C)'], water_df['Vapor Pressure (Mpa)'] * 1e6) # Pa, vapor pressure of water

temp_from_nu_water = interpolate.interp1d(water_df['Viscosity (m2/s)'][::-1], water_df['Temperature (C)'][::-1])

def friction_factor_smooth(Re, Re_crit=2300.):
   #assert Re > 0
   
   # TODO: This correlation produces a noticeable kink in the Re vs I plot. Use one without this problem?
   # From incropera_fundamentals_2007 p. 490.
   
   if Re > 0.:
      if Re < Re_crit:
         return friction_factor_laminar(Re)
      else:
         return (0.790 * log(Re) - 1.64)**(-2.)
         if Re > 5.e6:
            warnings.warn("Re above upper limit of regression: Re =", Re)
      #elif ((Re > Re_crit) and (Re < 2e4)):
         #return 0.316 * Re**(-0.25)
      #else:
         #return 0.184 * Re**(-0.2)
   else:
      return np.nan

def friction_factor_smooth_array(Re_array, Re_crit=2300.):
   if isinstance(Re_crit, int):
      Re_crit = float(Re_crit)
   if isinstance(Re_crit, float):
      Re_crit = Re_crit * np.ones(len(Re_array))
   f_array = np.zeros(len(Re_array))
   i = 0
   for Re in Re_array:
      f_array[i] = friction_factor_smooth(Re, Re_crit[i])
      i = i + 1
   
   return f_array

#def friction_factor_rough(Re, eps_o_D):
   #assert Re > 0
   #assert eps_o_D > 0
   
   ## TODO: This has not yet been tested. Test it.
   
   ## From more_analytical_2006.
   
   #a = eps_o_D / 3.7;
   #b = 1.257 / Re;
   #c = 1.7372;
   
   #x = exp(a / (b * c)) / (b * c);
   
   #f = 1 ./ (c * lambertw(x, 0) - a / b)**2; # Fanning friction factor
   #f = 4 * f; # Darcy friction factor

def friction_factor_laminar(Re):
   return 64. / Re

def I_fully_developed(f):
   return 0.36552*f**0.45867

def I_fully_developed_array(Re_array, Re_turb=4000.):
   if isinstance(Re_turb, int):
      Re_turb = float(Re_turb)
   if isinstance(Re_turb, float):
      Re_turb = Re_turb * np.ones(len(Re_array))
   i = 0
   I_array = np.zeros(len(Re_array))
   for Re in Re_array:
      if Re > Re_turb[i]:
         I_array[i] = I_fully_developed(friction_factor_smooth(Re))
      else:
         I_array[i] = np.nan
      
      i = i + 1
   
   return I_array

def Lambda_s(roughness, f, direction):
   # Data from powe_turbulence_1970.
   
   # TODO: Average with robertson_study_1965?
   if roughness == 'smooth':
      if direction == 'u':
         Lambda_s = 0.0791970187645
      elif direction == 'v':
         Lambda_s = 0.04558141583441
      elif direction == 'w':
         Lambda_s = 0.09033242035238
      else:
         sys.exit('Invalid direction.')
      
      assert(Lambda_s > 0)
   elif roughness == 'rough':
      if direction == 'u':
         Lambda_s = np.interp(f, [0.021259, 0.062848], [0.08037865532689, 0.05220755185866])
      elif direction == 'v':
         Lambda_s = np.interp(f, [0.021259, 0.062848], [0.07741487586879, 0.02030618096583])
      elif direction == 'w':
         Lambda_s = np.interp(f, [0.021259, 0.062848], [0.07178763776754, 0.02316331301628])
      else:
         sys.exit('Invalid direction.')
      
      assert(Lambda_s > 0)
   else:
      #sys.exit('Invalid roughness.')
      Lambda_s = np.nan
   
   return Lambda_s

def Lambda_s_array(roughness_array, f_array, direction):
   if isinstance(roughness_array, str):
      roughness = roughness_array
      roughness_array = []
      i = 0
      for f in f_array:
         roughness_array.append(roughness)
         i = i + 1
   
   i = 0
   Lambda_s_array = np.zeros(len(roughness_array))
   for roughness in roughness_array:
      Lambda_s(roughness, f_array[i], direction)
      
      i = i + 1
   
   return Lambda_s_array

def sigma_regression(predicted, actual):
   actualbar = np.average(actual)
   SS_tot    = np.sum((actual - actualbar)**2.)
   
   return sqrt((1. - coeff_of_determination(predicted, actual)**2.) * SS_tot / len(predicted))

def coeff_of_determination(predicted, actual):
   predicted = predicted.values
   actual    = actual.values
   
   actualbar = np.average(actual)
   #actualbar = np.sum(actual)/len(actual)
   SS_res    = np.sum((actual - predicted)**2.)
   #SS_reg    = np.sum((predicted - actualbar)**2.)
   SS_tot    = np.sum((actual - actualbar)**2.)
   
   #print actualbar
   #print (actual - predicted)**2.
   #print (actual - actualbar)**2.
   #print SS_res, SS_res, SS_tot, SS_res + SS_res - SS_tot
   #print predicted
   
   rsquared = 1 - SS_res / SS_tot
   #rsquared = SS_reg / SS_tot
   #rsquared = np.corrcoef(predicted, actual)[0, 1]**2.
   #print rsquared
   
   #assert((SS_res + SS_reg - SS_tot) < 1e-6)
   # https://en.wikipedia.org/wiki/Coefficient_of_determination
   # > In some cases the total sum of squares equals the sum of the two other sums of squares defined above
   
   #assert(rsquared >= 0.)
   # http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
   # > Note that it is possible to get a negative R-square for equations that do not contain a constant term.
   # https://stats.stackexchange.com/questions/12900/when-is-r-squared-negative
   
   assert(rsquared <= 1.)
   
   return rsquared

def summary_table(df_key, variable=None):
   # Also make this function do some error checking. This will ensure that you narrow down the data points responsible for problems more easily.
   # TODO: Find errors in the data by checking redundant data, e.g., comparing a calculated Re against the given Re.
   # TODO: Assert that turbulent, laminar can't appear in regime photo and regime L_b, Assert that breakup regimes can't appear in regime turb.
   # TODO: Check that I is only assigned for turbulent cases.
   
   key_array = []
   for key in df_key['key']:
      if not(key in key_array):
         key_array.append(key)
   
   if variable is None:
      if len(key_array) > 1:
         print str(len(key_array))+' studies ('+str(len(df_key))+' data points):'
         print key_array
      else:
         print df_key['key'][0]+' ('+str(len(df_key))+' data points):'
      
      i = 0
      total_regimes = 0
      combined_regime_array = []
      for regime_photo, regime_turb in zip(df_key['regime photo'], df_key['regime turb']):
         if isinstance(regime_photo, basestring):
            total_regimes = total_regimes + 1
            #print i, regime_photo
            #if regime_photo == 'bursting':
               #regime_photo = 'F2S'
            
            if (regime_photo == 'Rayleigh') and (regime_turb == 'laminar'):
               regime_photo = 'laminar Rayleigh'
            
            if (regime_photo == 'Rayleigh') and (regime_turb == 'transitional'):
               regime_photo = 'transitional Rayleigh'
            
            if (regime_photo == 'Rayleigh') and (regime_turb == 'turbulent'):
               regime_photo = 'turbulent Rayleigh'
            
            combined_regime_array.append(regime_photo)
         else:
            combined_regime_array.append(np.nan)
         
         i = i + 1
      
      i = 0
      for regime_L_b, regime_turb in zip(df_key['regime L_b'], df_key['regime turb']):
         if isinstance(regime_L_b, basestring):
            if (regime_L_b == 'Rayleigh') and (regime_turb == 'laminar'):
               regime_L_b = 'laminar Rayleigh'
            
            if (regime_L_b == 'Rayleigh') and (regime_turb == 'transitional'):
               regime_L_b = 'transitional Rayleigh'
               
            if (regime_L_b == 'Rayleigh') and (regime_turb == 'turbulent'):
               regime_L_b = 'turbulent Rayleigh'
            if isinstance(combined_regime_array[i], basestring):
               if not(combined_regime_array[i] == regime_L_b):
                  raise ValueError('Inconsistent regimes:', i, regime_L_b, combined_regime_array[i])
            else:
               combined_regime_array[i] = regime_L_b
               total_regimes = total_regimes + 1
         
         i = i + 1
      
      regime_list = []
      for regime in combined_regime_array:
         if not(regime in regime_list):
            regime_list.append(regime)
      
      print len(df_key[df_key['L_b/d_0'].notnull()]), 'xbavg/d_0 data points'
      print len(df_key[df_key['theta'].notnull()]), 'theta data points'
      print len(df_key[df_key['D_10/d_0'].notnull()]) + len(df_key[df_key['D_30/d_0'].notnull()]) + len(df_key[df_key['D_32/d_0'].notnull()]), 'D data points'
      print len(df_key[df_key['x_i/d_0'].notnull()]), 'x_i/d_0 data points'
      print len(df_key[df_key['v_d_bar/vp'].notnull()]), 'v_d_bar/vp data points'
      print len(df_key[df_key['x_trans/d_0'].notnull()]), 'x_trans/d_0 data points'
      print len(df_key[df_key['x_e/d_0'].notnull()]), 'x_e/d_0 data points'
      print len(df_key[df_key['photo filename'].notnull()]), 'photographs'
      print len(df_key[df_key['regime photo'].notnull()]), 'photographic regimes'
      print len(df_key[df_key['regime L_b'].notnull()]), 'breakup length regimes'
      print total_regimes, 'regimes in total'
      print 'regimes present:', regime_list
      print len(df_key[df_key['regime turb'] == 'turbulent']), 'turbulent data points'
      #print len(df_key[df_key['regime L_b'] == 'second wind-induced']) / len(df_key), 'second wind-induced'
      #print len(df_key[df_key['regime L_b'] == 'atomization']) / len(df_key), 'atomization'
   else:
      key_array = []
      for key in df_key['key']:
         if not(key in key_array):
            key_array.append(key)
      print 'Num. of applicable studies with '+variable+':', key_array
      print 'Number of data points:', len(df_key)
   num_measured_f = 0.
   for est_f in df_key['est f']:
      if not(est_f):
         num_measured_f = num_measured_f + 1.
   
   print str(int(1000. * num_measured_f / len(df_key)) / 10.) + '% measured friction factors'
   
   num_measured_turb = 0.
   for est_turb in df_key['est turb']:
      if not(est_turb):
         num_measured_turb = num_measured_turb + 1.
   
   print str(int(1000. * num_measured_turb / len(df_key)) / 10.) + '% measured hydrodynamic regimes'
   
   num_measured_nu = 0.
   for est_nu in df_key['est nu_s']:
      if not(est_nu):
         num_measured_nu = num_measured_nu + 1.
   
   print str(int(1000. * num_measured_nu / len(df_key)) / 10.) + '% measured viscosities'
   
   num_rough_tubes = 0.
   for roughness in df_key['roughness']:
      if roughness == 'rough':
         num_rough_tubes = num_rough_tubes + 1.
   
   df_key_turb = df_key[df_key['regime turb'] == 'turbulent']
   
   print str(int(1000. * num_rough_tubes / len(df_key)) / 10.) + '% rough pipes'
   print 'We_l0     = [', '%.1E' % Decimal(np.amin(df_key['We_l0'])), '--', '%.1E' % Decimal(np.amax(df_key['We_l0'])), ']'
   print 'Re_l0     = [', '%.1E' % Decimal(np.amin(df_key['Re_l0'])), '--', '%.1E' % Decimal(np.amax(df_key['Re_l0'])), ']'
   print 'Tu_0      = [', str(round(np.amin(df_key_turb['I_0']), 3)), '--', str(round(np.amax(df_key_turb['I_0']), 3)), ']'
   #print 'Tu_0      = [', np.amin(df_key['I_0']), '--', np.amax(df_key['I_0']),     ']'
   print 'Fr_0      = [', '%.1E' % Decimal(np.amin(df_key['Fr_0'])), '--', '%.1E' % Decimal(np.amax(df_key['Fr_0'])), ']'
   print 'rho_s     = [', str(round(np.amin(df_key['rho_s']), 1)), '--', str(round(np.amax(df_key['rho_s']), 1)), ']'
   print 'nu_s      = [', str(round(np.amin(df_key['nu_s']), 3)), '--', str(round(np.amax(df_key['nu_s']), 3)), ']'
   print 'K_c       = [', str(round(np.amin(df_key['K_c']), 3)), '--', str(round(np.amax(df_key['K_c']), 3)), ']'
   print 'Ma_g      = [', str(round(np.amin(df_key['Ma_g']), 3)), '--', str(round(np.amax(df_key['Ma_g']), 3)), ']'
   print 'L_0/d_0   = [', str(round(np.amin(df_key['L_0/d_0']), 1)), '--', str(round(np.amax(df_key['L_0/d_0']), 1)), ']'
   print 'c         = [', str(round(np.amin(df_key['c']), 1)), '--', str(round(np.amax(df_key['c']), 1)), ']'
   if variable is None:
      print 
   else:
      if variable[0:5] == 'theta': # plot in degrees for understandability
         theta_deg = df_key[variable] * (180 / np.pi)
         print variable+' = [', np.amin(theta_deg), ' ', np.amax(theta_deg), '] (deg.)\n'
      else:
         print variable+' = [', np.amin(df_key[variable]), ' ', np.amax(df_key[variable]), ']\n'
   
   assert(np.amin(df_key['We_l0']) > 0.)
   assert(np.amin(df_key['Re_l0']) > 0.)
   assert(np.amin(df_key['Fr_0']) > 0.)
   assert(np.amin(df_key['rho_s']) > 0.)
   assert(np.amin(df_key['nu_s']) > 0.)
   assert(np.amin(df_key['Ma_g']) > 0.)
   
   assert(np.amax(df_key['Ma_g']) > 0.)
   assert(np.amax(df_key['Ma_g']) < 1.0)
   
   if len(df_key[df_key['theta'].notnull()]) > 0:
      assert(np.amax(df_key['theta']) > 0.)
      assert(np.amax(df_key['theta']) < 1.)
   
   regime_array = []
   for regime in df_key['regime L_b']:
      if not(regime in regime_array):
         regime_array.append(regime)
         
         if isinstance(regime, basestring):
            if regime != 'Rayleigh':
               if regime != 'R2F':
                  if regime != 'first wind-induced':
                     if regime != 'F2S':
                        if regime != 'second wind-induced':
                           if regime != 'S2A':
                              if regime != 'atomization':
                                 if regime != 'bursting':
                                    if regime != 'R2B':
                                       if regime != 'R2S':
                                          if regime != 'F2R':
                                             raise ValueError('Invalid regime for '+df_key['key'][0]+':', regime)
   
   regime_array = []
   for regime in df_key['regime photo']:
      if not(regime in regime_array):
         regime_array.append(regime)
         
         if isinstance(regime, basestring):
            if regime != 'Rayleigh':
               if regime != 'R2F':
                  if regime != 'first wind-induced':
                     if regime != 'F2S':
                        if regime != 'second wind-induced':
                           if regime != 'S2A':
                              if regime != 'atomization':
                                 if regime != 'bursting':
                                    if regime != 'R2B':
                                       if regime != 'R2S':
                                          if regime != 'F2R':
                                             raise ValueError('Invalid regime for '+df_key['key'][0]+':', regime)
   
   #if key != 'eisenklam_flow_1958':
      #df_key_2 = df_key[df_key['regime L_b'].notnull()]
      #df_key_2 = df_key_2[df_key_2['Re_l0'] > 1e4]
      
      #for regime_L_b in df_key_2['regime L_b']:
         #if regime_L_b == 'first wind-induced':
            #raise ValueError('High Re regime F.')
      
      #df_key_2 = df_key[df_key['regime photo'].notnull()]
      #df_key_2 = df_key_2[df_key_2['Re_l0'] > 1e4]
      
      #for regime_photo in df_key_2['regime photo']:
         #if regime_photo == 'first wind-induced':
            #raise ValueError('High Re regime F.')
   
   df_key_2 = df_key[df_key['regime turb'] == 'turbulent']
   est_f_array = []
   for est_f in df_key_2['est f']:
      est_f_array.append(est_f)
   Re_l0_array = []
   for Re_l0 in df_key_2['Re_l0']:
      Re_l0_array.append(Re_l0)
   i = 0
   for I_0 in df_key_2['I_0']:
      if not(np.isnan(I_0)):
         if I_0 == 0:
            print df_key_2['I_0']
            raise ValueError('Zero turbulence intensity for turbulent case.')
         elif I_0 < I_fully_developed(friction_factor_laminar(Re_l0_array[i])):
            raise ValueError('Turbulent case has turbulence intensity lower than that of a smooth pipe in laminar flow.')
         elif I_0 < I_fully_developed(friction_factor_smooth(Re_l0_array[i])):
            if est_f_array[i]:
               raise ValueError('Turbulent case has turbulence intensity lower than that of a smooth pipe.')
      
      i = i + 1
   
   #L_b
   #\theta
   #D
   #x_i
   #photos

def latex_summary_table(df_key, variable):
   # DONE: Make LaTeX table generated automatically by the processing script.
   # WON'T: Use \uparrow, \downarrow, and \rightarrow for orientation in table.
   # WON'T: If writing out page figs, add slashes after periods to make sure the periods are not treated as the ends of sentences.
   
   table_filename = filename_escape(root_dir+'outputs/tables/summary_table_'+variable+'.tex')
   
   print 'Writing '+table_filename+'...'
   
   f = open(table_filename, 'w')
   
   f.write(r'\begin{table}'+'\n')
   f.write(r'\centering'+'\n')
   f.write(r'\begin{tabular}{rcc}'+'\n')
   f.write(r' & min. & max. \\'+'\n')
   f.write(r'\hline'+'\n')
   if variable == 'theta': # plot in degrees for understandability
      theta_deg = df_key[variable] * (180 / np.pi)
      f.write(convert_variable_to_latex(variable)+' & '+str(round(np.amin(theta_deg), 3))+' & '+str(round(np.amax(theta_deg), 3))+r' deg. \\'+'\n')
   else:
      if variable == 'D_32/d_0':
         precision = 3
      elif variable == 'L_b/d_0':
         precision = 1
      elif variable == 'v_d_bar/vp':
         precision = 3
      elif variable == 'x_i/d_0':
         precision = 1
      else:
         raise ValueError('Unable to determine precision for variable:', variable)
         
      
      f.write(convert_variable_to_latex(variable)+' & '+str(round(np.amin(df_key[variable]), precision))+' & '+str(round(np.amax(df_key[variable]), precision))+r' \\'+'\n')
   
   f.write(convert_variable_to_latex('We_l0')+' & '+convert_number_to_latex(np.amin(df_key['We_l0']))+' & '+convert_number_to_latex(np.amax(df_key['We_l0']))+r' \\'+'\n')
   f.write(convert_variable_to_latex('Re_l0')+' & '+convert_number_to_latex(np.amin(df_key['Re_l0']))+' & '+convert_number_to_latex(np.amax(df_key['Re_l0']))+r' \\'+'\n')
   f.write(convert_variable_to_latex('I_0')+' & '+str(round(np.amin(100. * df_key['I_0']), 1))+r'\% & '+str(round(100. * np.amax(df_key['I_0']), 1))+r'\% \\'+'\n')
   f.write(convert_variable_to_latex('rho_s')+' & '+str(round(np.amin(df_key['rho_s']), 1))+' & '+str(round(np.amax(df_key['rho_s']), 1))+r' \\'+'\n')
   f.write(convert_variable_to_latex('nu_s')+' & '+str(round(np.amin(df_key['nu_s']), 3))+' & '+str(round(np.amax(df_key['nu_s']), 3))+r' \\'+'\n')
   f.write(convert_variable_to_latex('L_0/d_0')+' & '+str(round(np.amin(df_key['L_0/d_0']), 1))+' & '+str(round(np.amax(df_key['L_0/d_0']), 1))+r' \\'+'\n')
   
   f.write(r'\end{tabular}'+'\n')
   f.write(r'\caption{Ranges of dependent and independent variables of compiled data for '+convert_variable_to_latex(variable)+'.}\n')
   
   key_array = []
   for key in df_key['key']:
      if not(key in key_array):
         key_array.append(key)
   
   f.write('Total data points: '+str(len(df_key))+r'\\ Studies considered: ')
   i = 0
   for key in key_array:
      i = i + 1
      
      num_points = len(df_key[df_key['key'] == key])
      
      # cite papers not dissertations
      if key == 'grant_newtonian_1965':
         key = 'grant_newtonian_1966'
      elif key == 'chen_mechanics_1962':
         key = 'chen_disintegration_1964'
      elif key == 'skrebkov_turbulentnyye_1963':
         key = 'skrebkov_turbulent_1966'
      if i == len(key_array):
         f.write(r'\citet{'+key+'} ('+str(num_points)+' data points)')
      else:
         f.write(r'\citet{'+key+'} ('+str(num_points)+' data points), ')
   f.write(r' \\'+'\n')
   
   num_measured_f = 0.
   for est_f in df_key['est f']:
      if not(est_f):
         num_measured_f = num_measured_f + 1.
   
   f.write(str(round(100. * num_measured_f / len(df_key), 1)) + r'\% measured friction factors \\'+'\n')
   
   num_measured_turb = 0.
   for est_turb in df_key['est turb']:
      if not(est_turb):
         num_measured_turb = num_measured_turb + 1.
   
   f.write(str(round(100. * num_measured_turb / len(df_key), 1)) + r'\% measured hydrodynamic regimes \\'+'\n')
   
   num_measured_nu = 0.
   for est_nu in df_key['est nu_s']:
      if not(est_nu):
         num_measured_nu = num_measured_nu + 1.
   
   f.write(str(round(100. * num_measured_nu / len(df_key), 1)) + r'\% measured viscosities \\'+'\n')
   
   num_rough_tubes = 0.
   for roughness in df_key['roughness']:
      if roughness == 'rough':
         num_rough_tubes = num_rough_tubes + 1.
   
   f.write(str(round(100. * num_rough_tubes / len(df_key), 1)) + r'\% rough pipes \\'+'\n')
   
   f.write(r'\label{tab:summary-table-'+variable+'}\n')
   f.write(r'\end{table}')
   
   f.close()

def plot_with_keys(df_variable, filename_header, input_variable_1, input_variable_2, plot_type='loglog', add_line=False, revno=None, actual=False, label_with_nozzle_type=False, fds_vers=None, filename_extra='', output_dir=root_dir):
   # TODO: Label different studies differently in these plots.
   # https://stats.stackexchange.com/questions/104622/what-does-an-actual-vs-fitted-graph-tell-us#comment202248_104622
   # > the convention of plotting the values that are fixed (conditional on predictors) on the x-axis and the values that are random on the y-axis
   # http://www.sciencedirect.com/science/article/pii/S0304380008002305
   # > Observed (in the y-axis) vs. predicted (in the x-axis) (OP) regressions should be used
   
   if input_variable_1[0:5] == 'theta': # plot in degrees for understandability
      df_variable[input_variable_1] = df_variable[input_variable_1] * (180. / np.pi)
   if input_variable_2[0:5] == 'theta':
      df_variable[input_variable_2] = df_variable[input_variable_2] * (180. / np.pi)
   
   #https://matplotlib.org/api/markers_api.html
   #marker_array = ["x", "+", ",", "v", "D", "^", "<", ">", "s", "p", "P", "*", "h", "H", "X", "d", "|", "_", "o", ".", "1", "2", "3", "4", "8"]
   
   # https://matplotlib.org/examples/lines_bars_and_markers/marker_reference.html
   marker_array = ['o', 's', 'D', 'x', '+', '>', '<', '^', 'v', '.']
   
   #plt.loglog(df_variable['We_l0'], df_variable['Re_l0'], 'x')
   if not(label_with_nozzle_type):
      extra_arr = df_variable['key']
   else:
      extra_arr = df_variable['nozzle']
   i = 0
   key_array = []
   for key, extra in zip(df_variable['key'], extra_arr):
      if not(label_with_nozzle_type):
         key_with_extra = key
      else:
         key_with_extra = key + extra
      
      if not(key_with_extra in key_array):
         #print key_with_extra
         key_array.append(key_with_extra)
         df_key = df_variable[df_variable['key'] == key]
         if label_with_nozzle_type:
            df_key = df_key[df_key['nozzle'] == extra]
         
         # cite papers not dissertations
         if key == 'grant_newtonian_1965':
            key = 'grant_newtonian_1966'
         elif key == 'chen_mechanics_1962':
            key = 'chen_disintegration_1964'
         elif key == 'skrebkov_turbulentnyye_1963':
            key = 'skrebkov_turbulent_1966'
         
         if not(label_with_nozzle_type):
            if len(df_key) == 1:
               key_label = r'\citet{'+escape_key(key)+r'} ('+str(len(df_key))+r' pt.)'
            else:
               key_label = r'\citet{'+escape_key(key)+r'} ('+str(len(df_key))+r' pts.)'
         else:
            try:
               int(extra[0])
            except: # if the nozzle name is not a number
               if len(df_key) == 1:
                  key_label = r'\citet{'+escape_key(key)+r'}, '+extra+' nozzle ('+str(len(df_key))+r' pt.)'
               else:
                  key_label = r'\citet{'+escape_key(key)+r'}, '+extra+' nozzle ('+str(len(df_key))+r' pts.)'
            else: # if the nozzle name is a number
               if len(df_key) == 1:
                  key_label = r'\citet{'+escape_key(key)+r'}, nozzle '+extra+' ('+str(len(df_key))+r' pt.)'
               else:
                  key_label = r'\citet{'+escape_key(key)+r'}, nozzle '+extra+' ('+str(len(df_key))+r' pts.)'
            
         #key_label = 'test'
         #print key_label
         if filename_header == 'correlation':
            error_2 = np.zeros(len(df_key[input_variable_2]))
            quantified_error = False
            len_stat = len(df_key[df_key[input_variable_2+' stat error'].notnull()])
            len_res  = len(df_key[df_key[input_variable_2+' resolution'].notnull()])
            if len_stat > 0:
               error_2 = error_2 + df_key[input_variable_2+' stat error']
               quantified_error = True
            if len_res > 0:
               error_2 = error_2 + df_key[input_variable_2+' resolution']
               quantified_error = True
            
            #if quantified_error:
               #print error_2, key
            
            if input_variable_2[0:5] == 'theta':
               error_2 = error_2 * (180 / np.pi)
         
         if plot_type == 'loglog':
            if 'quantified_error' in locals():
               if not(quantified_error):
                  plt.loglog(df_key[input_variable_1], df_key[input_variable_2], linestyle='None', marker=marker_array[i], label=key_label)
               else:
                  # https://matplotlib.org/examples/pylab_examples/log_demo.html
                  #print key
                  axes = plt.gca()
                  axes.set_xscale("log", nonposx='clip')
                  axes.set_yscale("log", nonposy='clip')
                  plt.errorbar(df_key[input_variable_1], df_key[input_variable_2], yerr=error_2, linestyle='None', marker=marker_array[i], label=key_label)
            else:
               plt.loglog(df_key[input_variable_1], df_key[input_variable_2], linestyle='None', marker=marker_array[i], label=key_label)
         elif plot_type == 'semilogx':
            plt.semilogx(df_key[input_variable_1], df_key[input_variable_2], linestyle='None', marker=marker_array[i], label=key_label)
         elif plot_type == 'semilogy':
            plt.semilogy(df_key[input_variable_1], df_key[input_variable_2], linestyle='None', marker=marker_array[i], label=key_label)
         elif plot_type == 'linear':
            if 'quantified_error' in locals():
               if not(quantified_error):
                  plt.plot(df_key[input_variable_1], df_key[input_variable_2], linestyle='None', marker=marker_array[i], label=key_label)
               else:
                  #print key
                  plt.errorbar(df_key[input_variable_1], df_key[input_variable_2], yerr=error_2, linestyle='None', marker=marker_array[i], label=key_label)
            else:
               plt.plot(df_key[input_variable_1], df_key[input_variable_2], linestyle='None', marker=marker_array[i], label=key_label)
         else:
            raise ValueError('Invalid plot_type:', plot_type)
         i = i + 1
   if input_variable_2 == 'rho_s':
      #min_var_1 = np.min(df_variable[input_variable_1])
      #max_var_1 = np.max(df_variable[input_variable_1])
      axes = plt.gca()
      min_var_1 = axes.get_xlim()[0]
      max_var_1 = axes.get_xlim()[1]
      rho_s_stp = rho_water(T_std) / rho_ideal_gas(P_atm, T_std, MW_air)
      if plot_type == 'loglog':
         plt.loglog([min_var_1, max_var_1], [rho_s_stp, rho_s_stp], label='water in air @ STP', ls='dashed')
      elif plot_type == 'linear':
         plt.plot([min_var_1, max_var_1], [rho_s_stp, rho_s_stp], label='water in air @ STP', ls='dashed')
      else:
         raise ValueError('Invalid plot_type:', plot_type)
   if input_variable_2 == 'I_0':
      axes = plt.gca()
      min_var_1 = axes.get_xlim()[0]
      max_var_1 = axes.get_xlim()[1]
      if plot_type == 'loglog':
         Re_l0_array = np.logspace(log(min_var_1) / log(10), log(max_var_1) / log(10), N_pts)
         I_0_array   = I_fully_developed_array(Re_l0_array, Re_turb)
         plt.loglog(Re_l0_array, I_0_array, label='FD flow, smooth pipe', ls='dashed')
      elif plot_type == 'semilogx':
         Re_l0_array = np.logspace(log(min_var_1) / log(10), log(max_var_1) / log(10), N_pts)
         I_0_array   = I_fully_developed_array(Re_l0_array, Re_turb)
         plt.semilogx(Re_l0_array, I_0_array, label='FD flow, smooth pipe', ls='dashed')
      elif plot_type == 'linear':
         Re_l0_array = np.linspace(min_var_1, max_var_1, N_pts)
         I_0_array   = I_fully_developed_array(Re_l0_array, Re_turb)
         plt.plot(Re_l0_array, I_0_array, label='FD flow, smooth pipe', ls='dashed')
      else:
         raise ValueError('Invalid plot_type:', plot_type)
   if input_variable_1 == 'Re_l0' and not(label_with_nozzle_type):
      plt.axvline(Re_trans, color='b', linestyle=None, marker=None, label=r'approx.\ onset of turbulence')
      plt.axvline(Re_turb, color='r', linestyle=None, marker=None, label=r'approx.\ fully turbulent')
   if (input_variable_2 == 'Re_l0') and not(label_with_nozzle_type):
      plt.axhline(Re_trans, color='b', linestyle=None, marker=None, label=r'approx.\ onset of turbulence')
      plt.axhline(Re_trans, color='r', linestyle=None, marker=None, label=r'approx.\ fully turbulent')
   
   # https://matplotlib.org/users/recipes.html#transparent-fancy-legends
   # http://matplotlib.1069221.n5.nabble.com/legend-symbols-is-duplicated-td10125.html
   plt.legend(loc='best', fancybox=True, framealpha=0.5, numpoints=1)
   
   extra = ''
   if add_line:
      # TODO: Change so that this always computes R^2 in the original coordinates, not log coordinates.
      #print 'R^2 =', coeff_of_determination(df_variable[input_variable_1], df_variable[input_variable_2])
      
      axes = plt.gca()
      box1 = TextArea('$R^2$ = '+str(int(1.e5*coeff_of_determination(df_variable[input_variable_1], df_variable[input_variable_2])) / 1.e5), textprops=dict(color="k"))
      anchored_box = AnchoredOffsetbox(loc=2,
                                       child=box1,
                                       frameon=False,
                                       bbox_to_anchor=(0., 1.),
                                       bbox_transform=axes.transAxes)
      #axes.add_artist(anchored_text)
      # loc=2, bbox_to_anchor=(0., -0.05),
      axes.add_artist(anchored_box)
      
      #extra = '\n\n$R^2$ = '+str(int(1.e5*coeff_of_determination(df_variable[input_variable_1], df_variable[input_variable_2])) / 1.e5)
      
      #min_var = min(np.min(df_variable[input_variable_1]), np.min(df_variable[input_variable_2]))
      #max_var = max(np.max(df_variable[input_variable_1]), np.max(df_variable[input_variable_2]))
      if ('eta_R_max' in input_variable_1) and ('eta_R_max' in input_variable_2):
         min_var = 0.
         max_var = 1.
      else:
         axes = plt.gca()
         min_var = min(axes.get_xlim()[0], axes.get_ylim()[0])
         max_var = max(axes.get_xlim()[1], axes.get_ylim()[1])
      
      plt.plot([min_var, max_var], [min_var, max_var], zorder=-1, linestyle='solid', linewidth=0.4, color='k')
      axes.set_xlim([min_var, max_var])
      axes.set_ylim([min_var, max_var])
   
   if not(revno is None):
      axes = plt.gca()
      #anchored_text = AnchoredText('Rev. '+revno, loc=4) #, frameon=False)
      box1 = TextArea('Rev.\ '+revno, textprops=dict(color="k"))
      if (filename_header == 'correlation') or (filename_header == 'other comparison') or actual:
         anchored_box = AnchoredOffsetbox(loc=4,
                                          child=box1,
                                          frameon=False,
                                          bbox_to_anchor=(1., 0.),
                                          bbox_transform=axes.transAxes)
         #axes.add_artist(anchored_text)
         # old: loc=4, bbox_to_anchor=(1., -0.20),
         axes.add_artist(anchored_box)
      else:
         anchored_box = AnchoredOffsetbox(loc=2,
                                          child=box1,
                                          frameon=False,
                                          bbox_to_anchor=(0., 1.),
                                          bbox_transform=axes.transAxes)
         axes.add_artist(anchored_box)
      
      #if extra == '':
         #extra = '\n\nSVN rev.\ '+revno
      #else:
         #extra = extra + ', SVN rev.\ '+revno
   
   if not(fds_vers is None):
      #axes = plt.gca()
      #box2 = TextArea(fds_vers, textprops=dict(color="k"))
      #anchored_box = AnchoredOffsetbox(loc=9,
                                       #child=box2,
                                       #frameon=False,
                                       #bbox_transform=axes.transAxes)
      ##axes.add_artist(anchored_text)
      ## old: loc=2, bbox_to_anchor=(0., 1.),
      #axes.add_artist(anchored_box)
      
      plt.title(fds_vers)
      
      #if extra == '':
         #extra = '\n\n' + fds_vers
      #else:
         #extra = extra + '\n' + fds_vers
   
   if input_variable_2[0:5] == 'theta':
      if input_variable_1[0:5] == 'theta':
         plt.xlabel(convert_variable_to_latex(input_variable_1) + ' (deg.)' + extra)
      else:
         plt.xlabel(convert_variable_to_latex(input_variable_1) + extra)
      
      if filename_header == 'correlation':
         plt.ylabel(convert_variable_to_latex(input_variable_2) + ' actual (deg.)')
      else:
         plt.ylabel(convert_variable_to_latex(input_variable_2) + ' (deg.)')
   else:
      plt.xlabel(convert_variable_to_latex(input_variable_1) + extra)
      if (filename_header == 'correlation') or (filename_header == 'other comparison') or actual:
         plt.ylabel(convert_variable_to_latex(input_variable_2) + ' actual')
      else:
         plt.ylabel(convert_variable_to_latex(input_variable_2))
   plt.grid()
   #plt.show()
   # http://matplotlib.org/users/annotations_intro.html
   # http://stackoverflow.com/a/33417697
   
   if input_variable_2 == 'rho_s':
      axes = plt.gca()
      #axes.set_ylim([1, axes.get_ylim()[1]])
      axes.set_ylim([0.9 * axes.get_ylim()[0], axes.get_ylim()[1]])
   if 'eta_R_max' in input_variable_1:
      axes.set_xlim([0., 1.])
   if 'eta_R_max' in input_variable_2:
      axes.set_ylim([0., 1.])
   #plt.tight_layout() # TODO: comment out for talk
   fig = plt.gcf()
   plt.savefig(output_dir+'outputs/figures/'+filename_escape(filename_header+'_'+input_variable_1+'_vs_'+input_variable_2+filename_extra+'.png'))
   #fig.set_size_inches(12., 8., forward=True) # poster
   fig.set_size_inches(6., 4., forward=True) # report
   mpl.rcParams.update(pgf_with_pdflatex_report)
   
   # Move legend out of the plot to have extra space.
   # https://stackoverflow.com/a/43439132
   plt.legend(bbox_to_anchor=(0., -0.15), loc='upper left', frameon=False, fontsize=11)
   
   #mpl.rcParams.update({"font.family": "serif"})
   plt.savefig(output_dir+'outputs/figures/'+filename_escape(filename_header+'_'+input_variable_1+'_vs_'+input_variable_2+filename_extra+'_report.pgf'), bbox_inches="tight")
   
   fig.set_size_inches(5.5, 4., forward=True) # paper
   plt.savefig(output_dir+'outputs/figures/'+filename_escape(filename_header+'_'+input_variable_1+'_vs_'+input_variable_2+filename_extra+'_paper.pgf'), bbox_inches="tight")
   
   fig.set_size_inches(4., 2.5, forward=True) # talk
   mpl.rcParams.update(pgf_with_pdflatex_talk)
   axes = plt.gca()
   axes.legend_.remove()
   plt.savefig(output_dir+'outputs/figures/'+filename_escape(filename_header+'_'+input_variable_1+'_vs_'+input_variable_2+filename_extra+'_talk.pgf'), bbox_inches="tight")
   #mpl.rcParams.update(pgf_with_pdflatex_report)
   plt.close()
   
   print 'Saved plot for', filename_header, input_variable_1, 'vs.', input_variable_2
   
   # Convert back to radians.
   if input_variable_1[0:5] == 'theta':
      df_variable[input_variable_1] = df_variable[input_variable_1] * (np.pi / 180.)
   if input_variable_2[0:5] == 'theta':
      df_variable[input_variable_2] = df_variable[input_variable_2] * (np.pi / 180.)

def parameter_space_plots(df_variable, variable, revno=None, trajectory=False, filename_extra=''):
   if not(trajectory):
      summary_table(df_variable, variable)
      plot_with_keys(df_variable, 'parameter_space_'+variable+filename_extra, 'We_l0', 'Re_l0', revno=revno)
      plot_with_keys(df_variable, 'parameter_space_'+variable+filename_extra, 'We_l0', 'rho_s', revno=revno)
      plot_with_keys(df_variable, 'parameter_space_'+variable+filename_extra, 'Re_l0', 'I_0', plot_type='semilogx', revno=revno)
   else:
      plot_with_keys(df_variable, 'parameter_space_'+variable+filename_extra, 'We_l0', 'Re_l0', revno=revno, label_with_nozzle_type=True)
      plot_with_keys(df_variable, 'parameter_space_'+variable+filename_extra, 'Fr_0', 'avg_x_b/d_0', plot_type='linear', revno=revno, label_with_nozzle_type=True)
      plot_with_keys(df_variable, 'parameter_space_'+variable+filename_extra, 'Fr_h_0', 'theta_0', plot_type='semilogx', revno=revno, label_with_nozzle_type=True)

def escape_key(key):
   return key.replace('_', 'TdTEi')
   #return key

def filename_escape(filename):
   return filename.replace('/d_0', 's').replace('/vp', 's').replace(' ', '_')

def convert_number_to_latex(number):
   number = '%.1E' % Decimal(number)
   number = number.replace('+0', '+').replace('-0', '-').replace('E+', r' \cdot 10^{').replace('E-', r' \cdot 10^{-')
   return '$'+number+'}$'

def convert_variable_to_latex(variable):
   if variable.split(' ')[-1] == 'predicted':
      extra = r' predicted'
   elif variable.split(' ')[-1] == 'actual':
      extra = r' actual'
   else:
      extra = r''
   
   variable = variable.split(' ')[0]
   
   # TODO: Change the variable names to match those in the text. E.g., L_\text{b} instead of L_b, etc. You get a compilation error if you do this right now.
   
   if variable == 'L_b/d_0':
      return r'$\langle x_\text{b} \rangle/d_0$'+extra
   elif variable == 'theta':
      return r'$\theta_\text{i}$'+extra
   elif variable == 'Re_l0':
      return r'$\text{Re}_{\ell0}$'+extra
   elif variable == 'We_l0':
      return r'$\text{We}_{\ell0}$'+extra
   elif variable == 'We_g0':
      return r'$\text{We}_\text{g0}$'+extra
   elif variable == 'rho_s':
      return r'$\rho_\ell/\rho_\text{g}$'+extra
   elif variable == 'nu_s':
      return r'$\nu_\ell/\nu_\text{g}$'+extra
   elif variable == 'I_0':
      return r'$\overline{\text{Tu}}_0$'+extra
   elif variable == 'x_i/d_0':
      return r'$\langle x_\text{i} \rangle/d_0$'+extra
   elif variable == 'D_32/d_0':
      return r'$D_{32}/d_0$'+extra
   elif variable == 'v_d_bar/vp':
      return r'$\langle v_\text{d} \rangle / \overline{v^\prime_0}$'+extra
   elif variable == 'eta_R_max':
      return r'$\eta_\text{max}$'+extra
   elif variable == 'theta_0':
      return r'$\theta_0$'+extra
   elif variable == 'xbavg/d_0':
      return r'$\langle x_\text{b} \rangle/d_0$'+extra # TODO: convert L_b/d_0 to this later
   elif variable == 'D_max/d_0':
      return r'$D_\text{max}/d_0$'+extra
   elif variable == 'Fr_0':
      return r'$\text{Fr}_0$'+extra
   elif variable == 'Fr_h_0':
      return r'$\text{Fr}_{h_0}$'+extra
   elif variable == 'L_0/d_0':
      return r'$L_0/d_0$'+extra
   else:
      raise ValueError('Unable to convert variable name to TeX:', variable)

#def round_to_nearest_in_list(x, x_list):
   
#def dissipation():

def spray_angle(Tu_0, We_l0):
   #return 2. * np.arctan(0.001007 * Tu_0**0.3469 * We_l0**0.4024)
   
   #with open(root_dir+'outputs/data/TSB_thetai.pickle') as f:
      #C_theta, C_We_l0, C_I_0 = pickle.load(f)
      
      #return 2. * np.arctan(C_theta * Tu_0**C_I_0 * We_l0**C_We_l0)
   
   with open(root_dir+'outputs/data/TSB_thetai.pickle') as f:
      C_theta, C_We_l0 = pickle.load(f)
      
      return 2. * np.arctan(C_theta * We_l0**C_We_l0)

def breakup_length(Tu_0, We_l0):
   #return 3.5453 * Tu_0**(-0.2623) * We_l0**0.3387
   
   with open(root_dir+'outputs/data/TSB_xbavg.pickle') as f:
      C_TSB, alpha_We, alpha_Re_l0_2WI, alpha_Tu_2WI, alpha_rho_s_2WI = pickle.load(f)
      
      return C_TSB * Tu_0**alpha_Tu_2WI * We_l0**alpha_We

def breakup_length_atomization(Tu_0, rho_s):
   #return 5.30525284643 * Tu_0**(-0.567852186033) * rho_s**(0.334893454839)
   
   with open(root_dir+'outputs/data/atomization_xbavg.pickle') as f:
      C_A, alpha_Tu_A, alpha_rho_s = pickle.load(f)
      
      return C_A * Tu_0**alpha_Tu_A * rho_s**alpha_rho_s

def atoi(text):
   return int(text) if text.isdigit() else text

def natural_keys(text):
   '''
   alist.sort(key=natural_keys) sorts in human order
   http://nedbatchelder.com/blog/200712/human_sorting.html
   (See Toothy's implementation in the comments)
   '''
   return [ atoi(c) for c in re.split('(\d+)', text) ]

def We_l0_crit(Tu_0, rho_s):
   # second-wind-induced-atomization boundary
   #return 1.97458901974 * Tu_0**(-1.02350307758) * rho_s**1.02781138701
   #return 3.17250862504 * Tu_0**(-0.876275337814) * rho_s**1.00350442958
   
   with open(root_dir+'outputs/data/atomization_boundary.pickle') as f:
      C_TSBtoA, alpha_Tu_2WItoA, alpha_rho_s_2WItoA = pickle.load(f)
      
      return C_TSBtoA * Tu_0**alpha_Tu_2WItoA * rho_s**alpha_rho_s_2WItoA

def We_l0_crit_TR(Tu_0):
   We_T_crit = 8.
   return We_T_crit * Tu_0**(-2.)

def Re_l0_crit_DT(We_l0, Re_x_trans=3.e5, C_LR=13.4):
   return (Re_x_trans - 3. * C_LR * We_l0) / (C_LR * We_l0**(1./2.))

def x_b_s(Tu_0_arr, We_l0_arr, rho_s):
   x_b_s_arr = np.array([])
   for Tu_0, We_l0 in zip(Tu_0_arr, We_l0_arr):
      if We_l0 < We_l0_crit(Tu_0, rho_s):
         x_b_s_arr = np.append(x_b_s_arr, breakup_length(Tu_0, We_l0))
      else:
         x_b_s_arr = np.append(x_b_s_arr, breakup_length_atomization(Tu_0, rho_s))
   
   return x_b_s_arr

def D_32_s(x_s, We_l0):
   # from wu_primary_1992 p. 312, modified to remove the assumption about the integral scale (see p. 305 for that).
   # see 2018-08-30 handwritten notes for the derivation
   return 0.54 * (x_s / We_l0**0.54)**0.57

def eta_max_func(C_d, rho_s, Fr_0, D_max_s, avg_x_b_s, theta_0, Fr_h_0, alpha):
   C_d_s = (3. / 2.) * (C_d / rho_s) * (Fr_0 / D_max_s) * (1. - alpha)**2.
   
   a = 1. + ((avg_x_b_s / Fr_0) * np.sin(theta_0) - 0.5 * (avg_x_b_s / Fr_0)**2. + 1. / Fr_h_0) * (C_d_s * np.cos(theta_0))**2.
   b = 1. + C_d_s * (np.sin(theta_0) - avg_x_b_s / Fr_0) * np.cos(theta_0)
   
   eta_max_max_s_analytical = (-a/b - np.real(lambertw(-np.exp(-a/b)/b, -1))) / C_d_s + avg_x_b_s * np.cos(theta_0) / Fr_0;
   eta_max_max_analytical   = eta_max_max_s_analytical / np.sqrt(1. + 2. / Fr_h_0)
   
   return eta_max_max_analytical

def m_to_in(d):
   return d / 25.4e-3

def in_to_m(d):
   return d * 25.4e-3

def ft_to_m(L):
   return L * 0.3048

# Munson, 2005, p. 438
K_L_reentrant       = 0.8
K_L_sharpedged      = 0.5
K_L_slightlyrounded = 0.2
K_L_wellrounded     = 0.04
K_L_guess           = K_L_sharpedged

# DONE: Change cavitation number to be calculated by a function, so you can try different versions of it.
def K_c(P_out, P_v, rho_l, Ubar_0, K_L, f, L_0s):
   # (P_in - P_v) / (0.5 * rho_l * Ubar_0**2)
   # P_in is rarely provided so instead I estimate it.
   return (P_out - P_v) / (0.5 * rho_l * Ubar_0**2) + 1 + K_L + f * L_0s

def f_FDS(D, D_FDS, gamma, sigma=None):
   if D <= D_FDS:
      if sigma is None:
         sigma = 2. / (sqrt(2. * np.pi) * log(2.) * gamma)
      f = (1. / (D * sigma * sqrt(2. * np.pi))) * exp(-0.5 * (log(D / D_FDS) / sigma)**2.)
   else:
      f = 0.693 * gamma * D_FDS**(-gamma) * D**(gamma - 1.) * exp(-0.693 * (D / D_FDS)**gamma)
   
   return f

def integrand(D, p, D_FDS, gamma, sigma=None):
   return D**p * f_FDS(D, D_FDS, gamma, sigma=sigma)

def charD(p, q, D_FDS, gamma, D_max=np.inf, sigma=None):
   numerator   = quad(integrand, 0., D_max, args=(p, D_FDS, gamma, sigma))[0]
   denominator = quad(integrand, 0., D_max, args=(q, D_FDS, gamma, sigma))[0]
   return (numerator / denominator)**(1./(p - q))

def f_lognormal(D, D_FDS, sigma):
   f = (1. / (D * sigma * sqrt(2. * np.pi))) * exp(-0.5 * (log(D / D_FDS) / sigma)**2.)
   
   return f

def integrand_lognormal(D, p, D_FDS, sigma):
   return D**p * f_lognormal(D, D_FDS, sigma)

def charD_lognormal(p, q, D_FDS, sigma, D_max=np.inf):
   numerator   = quad(integrand_lognormal, 1.e-3, D_max, args=(p, D_FDS, sigma))[0]
   denominator = quad(integrand_lognormal, 1.e-3, D_max, args=(q, D_FDS, sigma))[0]
   #print denominator
   return (numerator / denominator)**(1./(p - q))

def Wu_D32s(xs, We_l0):
   return 0.54 * (xs / (We_l0**(0.54)))**(0.57)

def stability_curve(d_0, rho_l, nu_l, sigma, filename, output_dir=root_dir):
   We_corner_FtoS = 1.e2
   N_Re = 1e3
   Re_arr = np.logspace(-1., 7., N_Re)
   We_arr = (rho_l * nu_l**2.) / (sigma * d_0) * Re_arr**2.
   xbavgs_arr = np.array([])
   Re_LNR = np.nan
   We_LNR = np.nan
   #D_bool   = False
   LR_bool  = False
   LNR_bool = False
   TR_bool  = False
   TSB_bool = False
   A_bool   = False
   prev_regime = None
   prev_We     = 1.
   xbavg_center = 4.e1
   
   K = 0.37 # clanet_transition_1999
   t = 0.   # tube thickness
   
   def dripping_We_func(We_lo):
      d_0     = rho_l * (Re * nu_l)**2. / (sigma * We_lo)
      Bo      = np.sqrt(rho_l * g * d_0**2. / (2. * sigma))
      Bo_outer  = np.sqrt(rho_l * g * (d_0 + t)**2. / (2. * sigma))
      We_crit = 4. * (Bo_outer / Bo) * (1. + K * Bo * Bo_outer - ((1. + K * Bo * Bo_outer)**2. - 1)**(1./2.))**2.
      return We_crit - We_lo
   
   # approximate center for xbavg, to place the text
   #xbavg_center = 
   
   We_guess = 4.
   for Re, We in zip(Re_arr, We_arr):
      if Re > 500:
         We_crit_dripping = 0.
      else:
         We_crit_dripping = fsolve(dripping_We_func, We_guess)[0]
      We_guess = We_crit_dripping
      
      if We < We_crit_dripping:
         # dripping
         xbavgs = 2.0 * ((6. * sigma) / (rho_l * g * d_0**2.))**(1./3.)
         prev_regime = 'dripping'
      else:
         if Re < Re_turb:
            if Re < Re_l0_crit_DT(We):
               # laminar Rayleigh
               xbavgs = 8.5 * (We**0.5 + 3. * We / Re)
               if not(LR_bool):
                  LR_bool = True
                  if (log(We) + log(prev_We)) > 0.:
                     plt.axvline(We, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
                     if We > 2:
                        plt.text(np.exp(0.5 * (log(We) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
                     else:
                        plt.text(np.exp(0.5 * (log(We) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center', fontsize='small', bbox=dict(boxstyle='square,pad=0.0',fc='white', ec='none'))
                  prev_We = We
               prev_regime = 'laminar Rayleigh'
            else:
               # laminar non-Rayleigh
               if (Re_LNR is np.nan) and (We_LNR is np.nan):
                  Re_LNR = Re
                  We_LNR = We
                  xbavgs_endLR = xbavgs
               #xbavg = 12. * (We_LNR**(3./2.) * 3. * We_LNR**2. / Re_LNR) / We # not sure why this doesn't work
               xbavgs = xbavgs_endLR * (We_LNR / We)
               if not(LNR_bool):
                  LNR_bool = True
                  plt.axvline(We, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
                  plt.text(np.exp(0.5 * (log(We) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
                  prev_We = We
               prev_regime = 'downstream transition'
         else:
            Tu = I_fully_developed(friction_factor_smooth(Re))
            
            if We < We_l0_crit_TR(Tu):
               # turbulent Rayleigh
               xbavgs = 3.27 * (We**0.5 + 3. * We / Re)
               if not(TR_bool):
                  TR_bool = True
                  plt.axvline(We, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
                  plt.text(np.exp(0.5 * (log(We) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
                  prev_We = We
               prev_regime = 'turbulent Rayleigh'
            elif We < We_l0_crit(Tu, 1000./1.2):
               # turbulent surface breakup
               xbavgs = breakup_length(Tu, We)
               if not(TSB_bool):
                  TSB_bool = True
                  plt.axvline(We, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
                  plt.text(np.exp(0.5 * (log(We) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
                  prev_We = We
               prev_regime = 'turbulent surface breakup'
            else:
               # atomization
               xbavgs = breakup_length_atomization(Tu, 1000./1.2)
               if not(A_bool):
                  A_bool = True
                  plt.axvline(We, marker=None, color='k', zorder=4, linewidth=0.8, linestyle='--')
                  plt.text(np.exp(0.5 * (log(We) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
                  prev_We = We
               prev_regime = 'atomization'
      
      xbavgs_arr = np.append(xbavgs_arr, xbavgs)
   
   plt.text(np.exp(0.5 * (log(1.e6) + log(prev_We))), xbavg_center, prev_regime, backgroundcolor='w', rotation=90, verticalalignment='center', horizontalalignment='center')
   
   ## apply some smoothing
   #N_smooth = int(round(0.05 * N_Re))
   #xbavg_arr2 = xbavg_arr
   #xbavg_arr = np.array([])
   #for i in range(0, len(xbavg_arr2)):
      #if i < N_smooth:
         #xbavg_arr = np.append(xbavg_arr, xbavg_arr2[i])
      #elif (i > N_smooth) and (i < (len(xbavg_arr2) - N_smooth)):
         #xbavg_arr = np.append(xbavg_arr, np.average(xbavg_arr2[i - N_smooth:i + N_smooth]))
      #else:
         #xbavg_arr = np.append(xbavg_arr, xbavg_arr2[i])
   
   plt.loglog(We_arr, xbavgs_arr, linestyle='-', marker=None)
   plt.xlim([1.e0, 1.e6])
   plt.ylim([1.e0, 1.e3])
   plt.xlabel(r'$\mathrm{We}_{\ell0}$')
   plt.ylabel(r'$\langle x_\text{b} \rangle / d_0$')
   plt.grid()
   
   fig = plt.gcf()
   fig.set_size_inches(6., 4., forward=True) # report
   plt.savefig(output_dir+'outputs/figures/stability_curve_'+filename+'.png')
   plt.savefig(output_dir+'outputs/figures/stability_curve_'+filename+'.pgf', bbox_inches="tight")
   fig.set_size_inches(5., 3., forward=True) # paper
   plt.savefig(output_dir+'outputs/figures/stability_curve_'+filename+'_paper.pgf', bbox_inches="tight")
   plt.close()

# TODO: Round to 3 sig. figs.
#def roundstr(number):
   #if abs(number) < 1e-8:
      ##return roundstre(number)
      #return r'\num{0}'
   #else:
      #number = "%.3f" % round(number, 3)
      #return r'\num{'+number+'}'

# http://randlet.com/blog/python-significant-figures-format/
def to_precision(x, p):
   """
   returns a string representation of x formatted with a precision of p

   Based on the webkit javascript implementation taken from here:
   https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
   """
   x = float(x)

   if x == 0.:
      return "0." + "0"*(p-1)

   out = []

   if x < 0:
      out.append("-")
      x = -x

   e = int(math.log10(x))
   tens = math.pow(10, e - p + 1)
   n = math.floor(x/tens)

   if n < math.pow(10, p - 1):
      e = e -1
      tens = math.pow(10, e - p+1)
      n = math.floor(x / tens)

   if abs((n + 1.) * tens - x) <= abs(n * tens -x):
      n = n + 1

   if n >= math.pow(10,p):
      n = n / 10.
      e = e + 1

   m = "%.*g" % (p, n)

   if e < -2 or e >= p:
      out.append(m[0])
      if p > 1:
         out.append(".")
         out.extend(m[1:p])
      out.append('e')
      if e > 0:
         out.append("+")
      out.append(str(e))
   elif e == (p -1):
      out.append(m)
   elif e >= 0:
      out.append(m[:e+1])
      if e+1 < len(m):
         out.append(".")
         out.extend(m[e+1:])
   else:
      out.append("0.")
      out.extend(["0"]*-(e+1))
      out.append(m)

   return "".join(out)

def roundstr(x, p=3, num=True):
   if num:
      return r'\num{'+to_precision(x, p)+'}'
   else:
      return to_precision(x, p)
