#!/usr/bin/env python
# -*- coding: utf-8 -*-

from jetbreakup import *
from string import join
import scipy.stats

metadata = [lastchangedby, lastchangedrevision, lastchangeddate]

print

# TODO: Generate "breakdown" table in validation paper automatically from this. Add extra fields to the database for the data quality issues mentioned so that the table can be generated automatically.
# TODO: Split into multiple files, one per study.

# TODO: Obtain more data for liquid-liquid systems, liquid jets in cross flow, liquid jets with co-flows
# TODO: chow_experimental_1980 doesn't have many data points but discusses many things done to have a high quality experiment. Add some more database fields inspired by this.
# TODO: Reject all laminar data that isn't fully developed.
# TODO: Can K_c be negative?
# TODO: Check that all basic quantities (e.g., Re, Fr, etc.) are defined

# WON'T: Change c to area contraction ratio rather than diameter contraction ratio (both \geq 0 here).
# TODO: check if need confidence interval or prediction interval in database
# TODO: \langle x_\text{b} \rangle ==> confidence interval for average
# TODO: \langle x_\text{i} \rangle ==> confidence interval for average

# TODO: Change variable names to be more consistent with later work. L_b ==> x_b, I_0 ==> Tu_0

# TODO: duffie_factors_1953/duffie_factors_1951
# TODO: richardson_mechanism_1954
# TODO: keith_liquid-liquid_1955/keith_drop_1951
# TODO: christiansen_breakup_1957/christiansen_influence_1957
# TODO: ryley_construction_1963
# TODO: smirnov_effect_1965
# TODO: meister_formation_1966/scheele_drop_1968/scheele_drop_1968-1/meister_drop_1969/meister_prediction_1969
# TODO: rao_drop_1966 ("capillary" but what is its length?)
# TODO: dotson_study_1967
# TODO: newman_preliminary_1967
# TODO: kroesser_stability_1968, kroesser_viscoelastic_1969
# TODO: takahashi_effect_1969, takahashi_laminar_1971/kitamura_stability_1986, takahashi_stability_1972/takahashi_breakup_1972/takahashi_laminar_1971, kitamura_stability_1982
# TODO: Add fenn_newtonian_1969/fenn_ambient_1968 ==> Should be useful to examine density ratio effects on LR-DT boundary.
# TODO: rutland_theoretical_1970 (see p. 1692R: 4 mm diameter, 30 cm length nozzle)
# TODO: skelland_dispersed_1971/minhas_dispersed_1969
# TODO: gordon_non-newtonian_1973 (both Newtonian and non-Newtonian fluids)
# TODO: parkin_production_1973 (p. 67: nozzle not long enough)
# TODO: skelland_jet_1974/johnson_jet_1973
# TODO: van_de_sande_air_1974 fig III.13 (pdf p. 41), fig. III.15 (turbulent Rayleigh, pdf p. 43), fig. III.17 (can't disambiguate, pdf p. 46, also van_de_sande_surface_1973 fig. 3)
# TODO: lafrance_capillary_1975 (most experimental results are forced, but p. 78 has some earlier unforced data), lafrance_drop_1974
# TODO: skelland_dispersed_1977

# TODO: kitamura_influence_1978
# fig. 8: \rhol/\rhog \approx 1.3, highest Re \approx 3000
# fig. 9: \rhol/\rhog \approx 0.8, highest Re \approx 1500

# TODO: van_den_akker_spontaneous_1980 (long enough for laminar fully developed flow?)
# TODO: abbott_analysis_1982
# TODO: mayer_zur_1993/mayer_coaxial_1994/mayer_rocket_1995
# TODO: eroglu_coaxial_1991/eroglu_initial_1991/eroglu_wave_1991/farago_morphological_1992/farago_parametric_1990
# TODO: olinger_lock-states_1993 (forced dripping)
# TODO: hardalupas_characteristics_1994
# TODO: Add amagai_frequency_1997/arai_surface_1999
# TODO: leroux_break-up_1997
# TODO: rhim_measurement_1998
# TODO: Add tamaki_effects_1998
# TODO: Add clanet_transition_1999
# TODO: Add blaisot_determination_2000, godelle_phase_2000
# TODO: Add godelle_symbolic_2000
# TODO: Add malot_experimental_2001
# TODO: blaisot_instabilities_2003
# TODO: aalburg_primary_2005/aalburg_deformation_2002 (has data with no cross-flow: fig. 4, figs. 7--9)
# TODO: Add salyers_spray_2010
# TODO: Add walker_effect_2012

# TODO: Add Ruff and Tseng data, which apparently has low \rhol/\rhog photographic regimes according to wu_effects_1995 p. 189, fig. 7 (pdf p. 15) / faeth_structure_1995 p. 117, fig. 1.5.
# TODO: Based on the same figure as the previous line's TODO, so does wu_effects_1995.

# TODO: Add radial droplet velocity uncertainty. wu_liquid_1992 p. 129
# TODO: Add droplet size uncertainty from wu_liquid_1992 too.
# TODO: Add breakup lengths from sauerwein_theoretische_1992 pp. 121--122?
# TODO: engelbert_breakup_1995
# TODO: mun_effects_1998 (L_0/d_0 = 25.6. Not long enough for laminar jets?)
# TODO: Add breakup lengths from rhim_measurement_1998?
# TODO: tang_laminar_2003/tang_cylindrical_2004
# TODO: badens_laminar_2005
# TODO: Add spray angle from osta_effect_2010 pdf p. 70
# TODO: umemura_two-valued_2011
# TODO: moallemi_breakup_2016 (see pdf p. 10)
# TODO: omocea_breakup_2016 (unclear on length of nozzle)
# TODO: rajendran_experimental_2017/rajendran_experimental_2012
# TODO: kiaoulias_evaluation_2019 (has pressure drop so I can better estimate f!)
# TODO: patrascu_dominant-wavelength_2019 (unclear on length of nozzle)

# TODO: Add atomization regime markers from palmer_water_1962 fig. 2 (also see report by same author).
# TODO: Add motion compensation key for photos. thorne_effect_1978, taylor_water_1983 (I later read that Hoyt and Taylor had two different setups with different methods of motion compensation?)
# TODO: Add transcription error. Break into error from image resolution, error from not knowing the center of the point due to its odd shape, and error from imprecise placement of the points.
# TODO: Get more spray angles from photos.
# TODO: Add Knudsen number. https://en.wikipedia.org/wiki/Knudsen_number

# TODO: Add other sources of error for spray angle.
# inter-rater reliability for droplet size: wu_measurements_1986 p. 945R (they call it "person-to-person bias")
# differences in definitions of spray angle:
# ~/reference/Engineering/Fluid dynamics/Multiphase/Liquid jets/Breakup/Spray angle/Definition/
# balewski_experimental_2010-1 discusses the effects of different threshold values
# TODO: Give people photos and ask them to measure the spray angle. Use this as an estimate of the "inter-rater reliability".

# TODO: Data validation for incrementing variables like pressure, velocity, Reynolds number, Weber number, etc. Do this for certain series, resetting when a new series is encountered?

# TODO: Look into papers on scale effects in hydraulic structures for other ideas.
# TODO: Add whether the flow is cavitating or not.
# TODO: Add field to mark as true when data for the wrong fluid was used due to the correct data being unavailable or otherwise.
# TODO: Add duration of jet in time.
# TODO: Add Jacob number for evaporation/flashing. See Incropera p. 376, kitamura_critical_1986, cleary_flashing_2007. Alternative: how close the temperature is to boiling. nezgada_effect_1970, nezgada_effect_1970-1
# TODO: Check whether spray angles are for the near or far fields.
# TODO: Characterize velocity profile. logan_flow_1963
# TODO: Add index of refraction. corey_droplet_1969 pdf p. 37, balewski_experimental_2008

# TODO: Classify all photos as turbulent or laminar at the nozzle. Add assertion to check this.
# TODO: Find surface transition location in all photos.
# TODO: Double check all orientations.
# TODO: List limits and number of data points for each dependent variable. Use these for the correlations developed.

# WON'T: Write function to mostly automatically determine the regime from the breakup length. Notice when "major" (i.e., not including noise) changes in slope occur. The other regime column is for visual classification. L regime vs V regime.

# MAYBE: For Mansour's data, estimate the friction factor from the u' and v' measurements. Compare this against smooth tubes. You can interpolate his u' and v' data to get the right TKE. You'll need the w' correlation for fully developed flow, however. ==> He doesn't have plane averaged TKE?
# TODO: Add more data from Phinney for various atmospheric densities.

# TODO: Source of SMD/MMD = 1.2? How does this vary with We, Re, I, etc.?

# TODO: Does Pandas automatically treat the missing data in D_10 correctly?

# MAYBE: Add information about cavitation inception. History of water? keller_influence_1980, ooi_scale_1985 (dissolved air content in ppm), brennen_cavitation_2013 p. 22
# MAYBE: Add other Mach number: U_0/c_l.
# MAYBE: Add distance to flat plate divided by nozzle diameter for breakup measurements. See takahashi_stability_1972 p. 75. li_study_2007 seems to be on the same topic, but I don't know if they come to the same conclusion.
# WON'T: Add some measure of total vorticity? Luis seems interested in vorticity generation in the flow.
# MAYBE: birkhoff_jets_1957 p. 330: > Also, dust in the air, the chemistry of the surrounding gas, and electrification may affect the final behavior decisively, by determining whether colliding drops amalgamate or rebound.
# MAYBE: Electrically ground the nozzle?
# TODO: Add ECN conditions:
# https://ecn.sandia.gov/diesel-spray-combustion/target-condition/spray-ab/
# https://ecn.sandia.gov/diesel-spray-combustion/experimental-diagnostics/
# https://ecn.sandia.gov/workshop/ECN2/SprayDev&Vap%20-%20Abstract.pdf
# TODO: Range testing nuisance variables. http://www.artofsoaking.com/2017/04/28/thoughts-on-water-gun-water-blaster-range-testing-2017/ : > As noted above, beyond the wind variable, I became more acutely aware of two other variables that can adversely affect stream performance: nozzle jitters and water gas content.
# Other nuisance variables: alberdi_strategies_2011 ("the fluid storage and delivery apparatus should not be treated as an afterthought [Alberdi et al. [2011], Ramesh et al. [2004]" from rouly_design_2015), bardi_engine_2012, meijer_engine_2012, payri_engine_2014, pei_engine_2015

# TODO: Change x_low, etc. to x_s_low as this is actually x_low/d_0.

# Can't disambiguate variables.
# WON'T: See if you can convert van_de_sande_surface_1973 fig. 3 data to spray angle.
# WON'T: Add branam_injection_2002 fig. 25 p. 8 if this is possible. Might not be able to get unambiguous Re, We, etc. values.

#########
# Rules #
#########

# - Use a different row for each photo. Only put breakup, etc. data in the first row. The other rows have We, Re, etc., but now breakup data. Order the photos by location in x, so typically the first row would have the photo of the jet at the nozzle.
# - When naming photo files, if the figure contains multiple photos but not labels for each photo individually, distinguish each with a letter in alphabetical order from the first photo in the figure.
# - All jet images are made horizontal with velocity from left to right.
# - Include images with unknown Re, etc. if you are confident they are of pipe jets.
# - Typically, if something is not mentioned, put np.nan in that field. Exceptions include fields defined by whether something was mentioned, e.g., the end checked field, and the trip field. (No trip mentioned is sufficient to assume they didn't use one.)
# - Use the full spray angle. Accept units in radians, not degrees or tan(\theta / 2).

data_file             = 'pipe-jet-breakup-data'
z                           = scipy.stats.norm.ppf(1 - (1 - interval_probability_level) / 2) # 2.5% on either side for a 5% confidence interval # DONE: change to t and use t when the number of data points is small, particularly for Grant and Chen
breakup_length_sigmas = 0.1291 # TODO: Obtain this from automatic analysis of the data. Update the paper automatically too.
cols = ['key', 'alt key',
        
        'liquid', 'gas', 'degas', 'regulation',
        
        'bend', 'flow straightener', 'screen', 'contraction shape', 'c', 'trip', 'L_r/d_0', 'eps_r/d_0',
        
        'L_0/d_0', 'roughness', 'f page fig', 'pipe material', 'est f', 'check FD',
        
        't/d_0', 'L_tip/d_0', 'end checked',
        
        'orientation', 'vibration isolation', 'd_chamber/d_0',
        
        'We_l0', 'Re_l0', 'I_0', 'Fr_0', 'rho_s', 'nu_s', 'K_c', 'Ma_g', 'Lambda_r_s', 'est rho_s', 'est nu_s',
        
        'regime photo', 'regime L_b', 'regime turb', 'est turb',
        
        'L_b/d_0', 'L_b/d_0 stat error', 'L_b/d_0 resolution', 'L_b method', 'L_b/d_0 page fig', 'L_b/d_0 transcription method',
        
        'theta', 'theta stat error', 'theta resolution', 'theta page fig', 'theta transcription method',
        
        'D_10/d_0', 'D_10/d_0 stat error', 'D_10/d_0 resolution', 'D_30/d_0', 'D_30/d_0 stat error', 'D_30/d_0 resolution', 'D_32/d_0', 'D_32/d_0 stat error', 'D_32/d_0 resolution', 'droplet x/d_0', 'D/d_0 page fig', 'D/d_0 transcription method', 'D/d_0 measurement method', 'v_d_bar/vp', 'v_d_bar/vp stat error', 'v_d_bar/vp resolution', 'v_d_bar/vp page fig', 'v_d_bar/vp transcription method', 'e_p', 'e_p page fig', 'e_p transcription method',
        
        'x_trans/d_0', 'x_trans/d_0 page fig', 'x_trans/d_0 transcription method',
        
        'x_i/d_0', 'x_i/d_0 stat error', 'x_i/d_0 resolution', 'x_i/d_0 page fig', 'x_i/d_0 transcription method',
        
        'x_e/d_0', 'x_e/d_0 stat error', 'x_e/d_0 resolution', 'x_e/d_0 page fig', 'x_e/d_0 transcription method',
        
        'photo filename', 'exposure time', 'flash', 'x_low', 'x_mid', 'x_high', 'photo page fig', 'spectrum', 'lighting', 'background color', 'grid', 'distance', 'f-stop', 'focal length', 'camera model', 'sensor']

########################
# asset_hydraulic_1951 #
########################

# DONE: Add transitional regimes (R2F, F2S, S2A).

asset_hydraulic_1951_csv = pd.read_csv('../data/asset_hydraulic_1951/asset_hydraulic_1951.csv', sep=',', header=0)

fluid_asset_hydraulic_1951          = asset_hydraulic_1951_csv['liquid']
Re_l0_asset_hydraulic_1951          = asset_hydraulic_1951_csv['Re_l0']
We_l0_asset_hydraulic_1951          = asset_hydraulic_1951_csv['We_l0']
Fr_0_asset_hydraulic_1951           = asset_hydraulic_1951_csv['Fr_0']
ts_asset_hydraulic_1951             = asset_hydraulic_1951_csv['t/d_0']
L_bs_asset_hydraulic_1951           = asset_hydraulic_1951_csv['L_b/d_0']
sigma_L_bs_asset_hydraulic_1951     = asset_hydraulic_1951_csv['sigma L_b/d_0'] * L_bs_asset_hydraulic_1951 / 100
regime_photo_asset_hydraulic_1951   = asset_hydraulic_1951_csv['photo regime']
regime_turb_asset_hydraulic_1951    = asset_hydraulic_1951_csv['turb regime']
Ubar_0_asset_hydraulic_1951         = asset_hydraulic_1951_csv['Ubar_0']
d_0_asset_hydraulic_1951            = asset_hydraulic_1951_csv['d_0']
D_10_asset_hydraulic_1951           = asset_hydraulic_1951_csv['D_10']
sigma_D_10_asset_hydraulic_1951     = asset_hydraulic_1951_csv['sigma D_10'] * D_10_asset_hydraulic_1951 / 100
photo_page_fig_asset_hydraulic_1951 = asset_hydraulic_1951_csv['photo page fig']
photo_filename_asset_hydraulic_1951 = asset_hydraulic_1951_csv['photo filename']
x_lows_asset_hydraulic_1951         = asset_hydraulic_1951_csv['x_low'] / d_0_asset_hydraulic_1951
x_mids_asset_hydraulic_1951         = asset_hydraulic_1951_csv['x_mid/L_b'] * L_bs_asset_hydraulic_1951
I_0_asset_hydraulic_1951            = I_fully_developed_array(Re_l0_asset_hydraulic_1951, Re_turb)
Lambda_r_s_asset_hydraulic_1951     = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_asset_hydraulic_1951, Re_turb), 'v')

N = 100 # Number of data points. Worst case scenario; see p. 5: > Nevertheless for this study 100 to 300 pictures were used to determine one length [or?] diameter
e_L_bs_asset_hydraulic_1951 = z * sigma_L_bs_asset_hydraulic_1951 / sqrt(N)
e_D_10_asset_hydraulic_1951 = z * sigma_D_10_asset_hydraulic_1951 / sqrt(N)

# p. 5: > The length was computed from Schiller's (8) equation which gave a conservative value for the turbulent case.
# schiller_entwicklung_1922's correlation according to durst_development_2005: L/D = 0.0288 Re
# But p. 6 suggests they used a different correlation. That is below.
L_0s_asset_hydraulic_1951 = 0.07 * Re_l0_asset_hydraulic_1951

# Estimate the liquid density.
i = len(Re_l0_asset_hydraulic_1951) - 1
rho_s_asset_hydraulic_1951 = np.zeros(len(Re_l0_asset_hydraulic_1951))
nu_s_asset_hydraulic_1951  = np.zeros(len(Re_l0_asset_hydraulic_1951))
K_c_asset_hydraulic_1951   = np.zeros(len(Re_l0_asset_hydraulic_1951))
Ma_g_asset_hydraulic_1951  = np.zeros(len(Re_l0_asset_hydraulic_1951))
rho_l_bar = 0.
nu_g_bar  = 0.
for Re_l0 in Re_l0_asset_hydraulic_1951[::-1]:
   if np.isnan(Re_l0):
      rho_s_asset_hydraulic_1951[i] = np.nan
      nu_s_asset_hydraulic_1951[i]  = np.nan
      K_c_asset_hydraulic_1951[i]   = np.nan
      Ma_g_asset_hydraulic_1951[i]  = np.nan
   else:
      fluid  = fluid_asset_hydraulic_1951[i]
      Ubar_0 = Ubar_0_asset_hydraulic_1951[i]
      d_0    = d_0_asset_hydraulic_1951[i]
      We_l0  = We_l0_asset_hydraulic_1951[i]
      
      nu_l = (Ubar_0 * d_0) / Re_l0
      
      if fluid == 'benzyl alcohol':
         # https://pubchem.ncbi.nlm.nih.gov/compound/benzyl_alcohol#section=Density
         rho_l = 1041.9; # kg/m^3
         
         # https://pubchem.ncbi.nlm.nih.gov/compound/benzyl_alcohol#section=Viscosity
         T = np.interp(rho_l * nu_l, (np.array([5.474, 2.760, 1.618, 1.055]) * 1e-3)[::-1], (np.array([25., 50., 75., 100.])[::-1])) # C
         rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
         nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
         
         # https://pubchem.ncbi.nlm.nih.gov/compound/benzyl_alcohol#section=Surface-Tension
         sigma = np.interp(T, np.array([20., 80]), np.array([39., 33.])) * 1e-3 # N/m
         
         # https://pubchem.ncbi.nlm.nih.gov/compound/benzyl_alcohol#section=Vapor-Pressure
         P_v = np.interp(T, np.array([20., 25., 57.8]), np.array([13.2, 12.5, 133.3])) # Pa
      
      if fluid == 'tide solution':
         # 3%; see p. 5
         
         # Assume density equals that of water? Use the average temperature of the 3 water data points.
         rho_l = rho_l_bar
         
         # Assume that the air viscosity is the average of the 3 points for water.
         nu_g  = nu_g_bar
         
         # Calculate a surface tension consistent with the above choices
         sigma = rho_l * Ubar_0**2 * d_0 / We_l0
         
         # Assume the vapor pressure is that of water. This is the worst-case scenario.
         # TODO: Change this to be more accurate?
         P_v = P_v_water(T)
      
      if fluid == 'water':
         T     = temp_from_nu_water(nu_l)
         rho_l = rho_water(T)
         sigma = sigma_water(T)
         rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
         nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
         P_v   = P_v_water(T)
         
         rho_l_bar = rho_l_bar + rho_l / 3.
         nu_g_bar  = nu_g_bar + nu_g / 3.
      
      We_l0_check = rho_l * Ubar_0**2 * d_0 / sigma
      
      #print abs(We_l0_check - We_l0) / We_l0
      assert(abs(We_l0_check - We_l0) / We_l0 < 0.03)
      
      #print T, rho_l / rho_g, abs(We_l0_check - We_l0) / We_l0
      
      rho_s_asset_hydraulic_1951[i] = rho_l / rho_g
      nu_s_asset_hydraulic_1951[i]  = nu_l / nu_g
      K_c_asset_hydraulic_1951[i]   = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_guess, friction_factor_smooth(Re_l0), L_0s_asset_hydraulic_1951[i])
      Ma_g_asset_hydraulic_1951[i]  = Ubar_0 / c_ideal_gas(gamma_air(T), T + T_zero, MW_air)
      #print 0.5 * rho_l * Ubar_0**2, P_atm - P_v , K_c_asset_hydraulic_1951[i]
   
   if not(np.isnan(K_c_asset_hydraulic_1951[i])):
      assert(K_c_asset_hydraulic_1951[i] > 0)
   
   i = i - 1

df_asset_hydraulic_1951 = pd.DataFrame({'key': 'asset_hydraulic_1951',
                                        'liquid' : fluid_asset_hydraulic_1951})
df_asset_hydraulic_1951['alt key']    = np.nan
df_asset_hydraulic_1951['gas']        = 'air'
df_asset_hydraulic_1951['degas']      = np.nan
df_asset_hydraulic_1951['regulation'] = 'pressure regulator' # p. 5: > the air pressure above the liquid was kept uniform by a pressure regulator

df_asset_hydraulic_1951['bend']              = np.nan
df_asset_hydraulic_1951['flow straightener'] = np.nan # Was not mentioned on p. 5.
df_asset_hydraulic_1951['contraction shape'] = np.nan
df_asset_hydraulic_1951['screen']            = np.nan
df_asset_hydraulic_1951['c']                 = np.nan
df_asset_hydraulic_1951['trip']              = False
df_asset_hydraulic_1951['L_r/d_0']           = 0
df_asset_hydraulic_1951['eps_r/d_0']         = np.nan

df_asset_hydraulic_1951['L_0/d_0']       = L_0s_asset_hydraulic_1951 
df_asset_hydraulic_1951['roughness']     = 'smooth'
df_asset_hydraulic_1951['f page fig']    = np.nan
df_asset_hydraulic_1951['pipe material'] = 'glass' # p. 5: "B. Methods. \n The glass nozzles, mounted in a vertical position
df_asset_hydraulic_1951['est f']         = True
df_asset_hydraulic_1951['check FD']      = False

df_asset_hydraulic_1951['t/d_0']       = ts_asset_hydraulic_1951 # Estimated from the photos.
df_asset_hydraulic_1951['L_tip/d_0']   = np.nan
df_asset_hydraulic_1951['end checked'] = False

df_asset_hydraulic_1951['orientation']         = 'down' # p. 6 : "the nozzle faced downwards
df_asset_hydraulic_1951['vibration isolation'] = True # p. 5: > The tubing was firmly mounted by supports which were mechanically insulated from the floor by felt one-inch thick. It was found that ordinary mechanical disturbances due to walking, for example, did not affect the jets.
df_asset_hydraulic_1951['d_chamber/d_0']       = np.inf

df_asset_hydraulic_1951['We_l0']            = We_l0_asset_hydraulic_1951
df_asset_hydraulic_1951['We_l0 resolution'] = np.nan
df_asset_hydraulic_1951['Re_l0']            = Re_l0_asset_hydraulic_1951
df_asset_hydraulic_1951['Re_l0 resolution'] = np.nan
df_asset_hydraulic_1951['I_0']              = I_0_asset_hydraulic_1951
df_asset_hydraulic_1951['I_0 resolution']   = np.nan
df_asset_hydraulic_1951['Fr_0']             = Fr_0_asset_hydraulic_1951
df_asset_hydraulic_1951['Fr_0 resolution']  = np.nan
df_asset_hydraulic_1951['rho_s']            = rho_s_asset_hydraulic_1951
df_asset_hydraulic_1951['rho_s resolution'] = np.nan
df_asset_hydraulic_1951['nu_s']             = nu_s_asset_hydraulic_1951
df_asset_hydraulic_1951['nu_s resolution']  = np.nan
df_asset_hydraulic_1951['K_c']              = K_c_asset_hydraulic_1951
df_asset_hydraulic_1951['K_c resolution']   = np.nan
df_asset_hydraulic_1951['Ma_g']             = Ma_g_asset_hydraulic_1951
df_asset_hydraulic_1951['Ma_g resolution']  = np.nan
df_asset_hydraulic_1951['Lambda_r_s']       = Lambda_r_s_asset_hydraulic_1951
df_asset_hydraulic_1951['est rho_s']        = True
df_asset_hydraulic_1951['est nu_s']         = True

df_asset_hydraulic_1951['regime photo']                 = regime_photo_asset_hydraulic_1951
df_asset_hydraulic_1951['regime L_b']                   = np.nan
df_asset_hydraulic_1951['regime turb']                  = regime_turb_asset_hydraulic_1951
df_asset_hydraulic_1951['est turb']                     = False

df_asset_hydraulic_1951['L_b/d_0']                      = L_bs_asset_hydraulic_1951
df_asset_hydraulic_1951['L_b/d_0 stat error']           = e_L_bs_asset_hydraulic_1951
df_asset_hydraulic_1951['L_b/d_0 resolution']           = np.nan
df_asset_hydraulic_1951['L_b method']                   = 'photographic'
df_asset_hydraulic_1951['L_b/d_0 page fig']             = 'p. 11, tab. 2'
df_asset_hydraulic_1951['L_b/d_0 transcription method'] = 'table'

df_asset_hydraulic_1951['theta']                      = np.nan
df_asset_hydraulic_1951['theta stat error']           = np.nan
df_asset_hydraulic_1951['theta resolution']           = np.nan
df_asset_hydraulic_1951['theta page fig']             = np.nan
df_asset_hydraulic_1951['theta transcription method'] = np.nan

df_asset_hydraulic_1951['D_10/d_0']                   = D_10_asset_hydraulic_1951 / d_0_asset_hydraulic_1951
df_asset_hydraulic_1951['D_10/d_0 stat error']        = e_D_10_asset_hydraulic_1951
df_asset_hydraulic_1951['D_10/d_0 resolution']        = np.nan
df_asset_hydraulic_1951['D_30/d_0']                   = np.nan
df_asset_hydraulic_1951['D_30/d_0 stat error']        = np.nan
df_asset_hydraulic_1951['D_30/d_0 resolution']        = np.nan
df_asset_hydraulic_1951['D_32/d_0']                   = np.nan
df_asset_hydraulic_1951['D_32/d_0 stat error']        = np.nan
df_asset_hydraulic_1951['D_32/d_0 resolution']        = np.nan
df_asset_hydraulic_1951['droplet x/d_0']              = L_bs_asset_hydraulic_1951
df_asset_hydraulic_1951['D/d_0 page fig']             = 'p. 13, tab. 3'
df_asset_hydraulic_1951['D/d_0 transcription method'] = 'table'
df_asset_hydraulic_1951['D/d_0 measurement method']   = 'photographic'
df_asset_hydraulic_1951['v_d_bar/vp']                        = np.nan
df_asset_hydraulic_1951['v_d_bar/vp stat error']             = np.nan
df_asset_hydraulic_1951['v_d_bar/vp resolution']             = np.nan
df_asset_hydraulic_1951['v_d_bar/vp page fig']               = np.nan
df_asset_hydraulic_1951['v_d_bar/vp transcription method']   = np.nan
df_asset_hydraulic_1951['e_p']                        = np.nan
df_asset_hydraulic_1951['e_p page fig']               = np.nan
df_asset_hydraulic_1951['e_p transcription method']   = np.nan

df_asset_hydraulic_1951['x_trans/d_0']                      = np.nan
df_asset_hydraulic_1951['x_trans/d_0 page fig']             = np.nan
df_asset_hydraulic_1951['x_trans/d_0 transcription method'] = np.nan

df_asset_hydraulic_1951['x_i/d_0']                      = np.nan
df_asset_hydraulic_1951['x_i/d_0 stat error']           = np.nan
df_asset_hydraulic_1951['x_i/d_0 resolution']           = np.nan
df_asset_hydraulic_1951['x_i/d_0 page fig']             = np.nan
df_asset_hydraulic_1951['x_i/d_0 transcription method'] = np.nan

df_asset_hydraulic_1951['x_e/d_0']                      = np.nan
df_asset_hydraulic_1951['x_e/d_0 stat error']           = np.nan
df_asset_hydraulic_1951['x_e/d_0 resolution']           = np.nan
df_asset_hydraulic_1951['x_e/d_0 page fig']             = np.nan
df_asset_hydraulic_1951['x_e/d_0 transcription method'] = np.nan

df_asset_hydraulic_1951['photo filename']   = photo_filename_asset_hydraulic_1951
df_asset_hydraulic_1951['exposure time']    = 10e-6 # s, p. 5
df_asset_hydraulic_1951['flash']            = True # p. 5
df_asset_hydraulic_1951['x_low']            = x_lows_asset_hydraulic_1951
df_asset_hydraulic_1951['x_mid']            = x_mids_asset_hydraulic_1951
df_asset_hydraulic_1951['x_high']           = np.nan
df_asset_hydraulic_1951['photo page fig']   = photo_page_fig_asset_hydraulic_1951
df_asset_hydraulic_1951['spectrum']         = 'visible'
df_asset_hydraulic_1951['lighting']         = 'diffuse background' # p. 5
df_asset_hydraulic_1951['background color'] = 'light'
df_asset_hydraulic_1951['grid']             = False
df_asset_hydraulic_1951['distance']         = np.nan
df_asset_hydraulic_1951['f-stop']           = np.nan
df_asset_hydraulic_1951['focal length']     = np.nan
df_asset_hydraulic_1951['camera model']     = np.nan
df_asset_hydraulic_1951['sensor']           = '16 mm film' # p. 5

df_asset_hydraulic_1951 = df_asset_hydraulic_1951[cols]
summary_table(df_asset_hydraulic_1951)
df_jet_breakup = df_asset_hydraulic_1951

########################
# betchov_breakup_1955 #
########################

# Issues with this data:
# 1. r in table V should be 0.2 mm, not 0.02. 0.2 mm is what is written in table II. This was checked for the sigma = 74 case, which appears to be water. The Ohnesorge number matches.

# betchov_breakup_1955 fig. 19, pdf p. 59
# p. 18 (pdf p. 29) has some details.
# Fluid? 1 cs viscosity. Table I (p. 30, pdf p. 71)
# Definition of W? Based on air density
# Table V (p. 36, pdf p. 77)

betchov_breakup_1955_csv = pd.read_csv('../data/betchov_breakup_1955/betchov_breakup_1955.csv', sep=',', header=0)

liquid_betchov_breakup_1955         = betchov_breakup_1955_csv['fluid']
nu_l_betchov_breakup_1955           = betchov_breakup_1955_csv['nu_l (m^2/s)']
sigma_betchov_breakup_1955          = betchov_breakup_1955_csv['sigma (mN/mm)'] * 1e-3 # N/m
rho_s_betchov_breakup_1955          = 1. / betchov_breakup_1955_csv['rho_g/rho_l']
We_g0_betchov_breakup_1955          = 2. * betchov_breakup_1955_csv['We (r)']
Re_approx_betchov_breakup_1955      = 2. * betchov_breakup_1955_csv['Re (r)']
page_fig_betchov_breakup_1955       = betchov_breakup_1955_csv['page fig']
regime_photo_betchov_breakup_1955   = betchov_breakup_1955_csv['regime']
photo_filename_betchov_breakup_1955 = betchov_breakup_1955_csv['photo filename']
photo_page_fig_betchov_breakup_1955 = betchov_breakup_1955_csv['photo page fig']

# Table II, p. 31
L_0 = 2.e-2
d_0 = 0.4e-3

We_l0_betchov_breakup_1955 = We_g0_betchov_breakup_1955 * rho_s_betchov_breakup_1955

Re_l0_betchov_breakup_1955       = []
Fr_0_betchov_breakup_1955        = []
nu_s_betchov_breakup_1955        = []
K_c_betchov_breakup_1955         = []
Ma_g_betchov_breakup_1955        = []
regime_turb_betchov_breakup_1955 = []
for We_l0, liquid, Re_approx, rho_s, nu_l, sigma in zip(We_l0_betchov_breakup_1955, liquid_betchov_breakup_1955, Re_approx_betchov_breakup_1955, rho_s_betchov_breakup_1955, nu_l_betchov_breakup_1955, sigma_betchov_breakup_1955):
   if liquid == 'water':
      rho_l  = 1000 # kg/m^3
      Ubar_0 = np.sqrt((sigma * We_l0) / (rho_l * d_0))
      P_v    = P_v_water(T)
   elif liquid == 'Dow 1 cs':
      rho_l = 820 # kg/m^3
      Ubar_0 = np.sqrt((sigma * We_l0) / (rho_l * d_0))
      P_v = np.nan
   else:
      Ubar_0 = Re_approx * nu_l / d_0
      rho_l  = sigma * We_l0 / (Ubar_0**2. * d_0)
      P_v = np.nan
   
   rho_g = rho_l / rho_s
   nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
   
   Re_l0 = Ubar_0 * d_0 / nu_l
   
   Re_l0_betchov_breakup_1955.append(Re_l0)
   Fr_0_betchov_breakup_1955.append(Ubar_0**2 / (g * d_0))
   nu_s_betchov_breakup_1955.append(nu_l / nu_g)
   #K_c_betchov_breakup_1955.append((P_atm - P_v) / (0.5 * rho_l * Ubar_0**2))
   K_c_betchov_breakup_1955.append(K_c(P_atm, P_v, rho_l, Ubar_0, K_L_guess, friction_factor_smooth(Re_l0), L_0/d_0))
   Ma_g_betchov_breakup_1955.append(Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air))
   
   if Re_l0 > Re_turb:
      regime_turb_betchov_breakup_1955.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_betchov_breakup_1955.append('transitional')
   else:
      regime_turb_betchov_breakup_1955.append('laminar')

I_0_betchov_breakup_1955        = I_fully_developed_array(Re_l0_betchov_breakup_1955, Re_turb)
Lambda_r_s_betchov_breakup_1955 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_betchov_breakup_1955, Re_turb), 'v')

df_betchov_breakup_1955 = pd.DataFrame({'key': 'betchov_breakup_1955',
                                       'liquid': liquid_betchov_breakup_1955})
df_betchov_breakup_1955['alt key']    = np.nan
df_betchov_breakup_1955['gas']        = 'air'
df_betchov_breakup_1955['degas']      = np.nan
df_betchov_breakup_1955['regulation'] = np.nan

df_betchov_breakup_1955['bend']              = np.nan
df_betchov_breakup_1955['flow straightener'] = np.nan
df_betchov_breakup_1955['screen']            = np.nan
df_betchov_breakup_1955['contraction shape'] = np.nan
df_betchov_breakup_1955['c']                 = np.nan
df_betchov_breakup_1955['trip']              = np.nan
df_betchov_breakup_1955['L_r/d_0']           = 0
df_betchov_breakup_1955['eps_r/d_0']         = np.nan

df_betchov_breakup_1955['L_0/d_0']       = L_0/d_0
df_betchov_breakup_1955['roughness']     = 'smooth'
df_betchov_breakup_1955['f page fig']    = np.nan
df_betchov_breakup_1955['pipe material'] = 'hypodermic needle' # p. 31, tab. II
df_betchov_breakup_1955['est f']         = True
df_betchov_breakup_1955['check FD']      = np.nan

df_betchov_breakup_1955['t/d_0']       = np.nan
df_betchov_breakup_1955['L_tip/d_0']   = np.nan
df_betchov_breakup_1955['end checked'] = np.nan

df_betchov_breakup_1955['orientation']         = np.nan
df_betchov_breakup_1955['vibration isolation'] = np.nan
df_betchov_breakup_1955['d_chamber/d_0']       = np.nan

df_betchov_breakup_1955['We_l0']      = We_l0_betchov_breakup_1955
df_betchov_breakup_1955['Re_l0']      = Re_l0_betchov_breakup_1955
df_betchov_breakup_1955['I_0']        = I_0_betchov_breakup_1955
df_betchov_breakup_1955['Fr_0']       = Fr_0_betchov_breakup_1955
df_betchov_breakup_1955['rho_s']      = rho_s_betchov_breakup_1955
df_betchov_breakup_1955['nu_s']       = nu_s_betchov_breakup_1955
df_betchov_breakup_1955['K_c']        = K_c_betchov_breakup_1955
df_betchov_breakup_1955['Ma_g']       = Ma_g_betchov_breakup_1955
df_betchov_breakup_1955['Lambda_r_s'] = Lambda_r_s_betchov_breakup_1955
df_betchov_breakup_1955['est rho_s']  = False
df_betchov_breakup_1955['est nu_s']   = True

df_betchov_breakup_1955['regime photo'] = regime_photo_betchov_breakup_1955
df_betchov_breakup_1955['regime L_b']   = np.nan
df_betchov_breakup_1955['regime turb']  = regime_turb_betchov_breakup_1955
df_betchov_breakup_1955['est turb']     = True

df_betchov_breakup_1955['L_b/d_0']                      = np.nan
df_betchov_breakup_1955['L_b/d_0 stat error']           = np.nan
df_betchov_breakup_1955['L_b/d_0 resolution']           = np.nan
df_betchov_breakup_1955['L_b method']                   = np.nan
df_betchov_breakup_1955['L_b/d_0 page fig']             = np.nan
df_betchov_breakup_1955['L_b/d_0 transcription method'] = np.nan

df_betchov_breakup_1955['theta']                      = np.nan
df_betchov_breakup_1955['theta stat error']           = np.nan
df_betchov_breakup_1955['theta resolution']           = np.nan
df_betchov_breakup_1955['theta page fig']             = np.nan
df_betchov_breakup_1955['theta transcription method'] = np.nan

df_betchov_breakup_1955['D_10/d_0']                   = np.nan
df_betchov_breakup_1955['D_10/d_0 stat error']        = np.nan
df_betchov_breakup_1955['D_10/d_0 resolution']        = np.nan
df_betchov_breakup_1955['D_32/d_0']                   = np.nan
df_betchov_breakup_1955['D_32/d_0 stat error']        = np.nan
df_betchov_breakup_1955['D_32/d_0 resolution']        = np.nan
df_betchov_breakup_1955['D_30/d_0']                   = np.nan
df_betchov_breakup_1955['D_30/d_0 stat error']        = np.nan
df_betchov_breakup_1955['D_30/d_0 resolution']        = np.nan
df_betchov_breakup_1955['droplet x/d_0']              = np.nan
df_betchov_breakup_1955['D/d_0 page fig']             = np.nan
df_betchov_breakup_1955['D/d_0 transcription method'] = np.nan
df_betchov_breakup_1955['D/d_0 measurement method']   = np.nan
df_betchov_breakup_1955['v_d_bar/vp']                        = np.nan
df_betchov_breakup_1955['v_d_bar/vp stat error']             = np.nan
df_betchov_breakup_1955['v_d_bar/vp resolution']             = np.nan
df_betchov_breakup_1955['v_d_bar/vp page fig']               = np.nan
df_betchov_breakup_1955['v_d_bar/vp transcription method']   = np.nan
df_betchov_breakup_1955['e_p']                        = np.nan
df_betchov_breakup_1955['e_p page fig']               = np.nan
df_betchov_breakup_1955['e_p transcription method']   = np.nan

df_betchov_breakup_1955['x_trans/d_0']                      = np.nan
df_betchov_breakup_1955['x_trans/d_0 page fig']             = np.nan
df_betchov_breakup_1955['x_trans/d_0 transcription method'] = np.nan

df_betchov_breakup_1955['x_i/d_0']                      = np.nan
df_betchov_breakup_1955['x_i/d_0 stat error']           = np.nan
df_betchov_breakup_1955['x_i/d_0 resolution']           = np.nan
df_betchov_breakup_1955['x_i/d_0 page fig']             = np.nan
df_betchov_breakup_1955['x_i/d_0 transcription method'] = np.nan

df_betchov_breakup_1955['x_e/d_0']                      = np.nan
df_betchov_breakup_1955['x_e/d_0 stat error']           = np.nan
df_betchov_breakup_1955['x_e/d_0 resolution']           = np.nan
df_betchov_breakup_1955['x_e/d_0 page fig']             = np.nan
df_betchov_breakup_1955['x_e/d_0 transcription method'] = np.nan

df_betchov_breakup_1955['photo filename']   = photo_filename_betchov_breakup_1955
df_betchov_breakup_1955['exposure time']    = np.nan
df_betchov_breakup_1955['flash']            = np.nan
df_betchov_breakup_1955['x_low']            = np.nan
df_betchov_breakup_1955['x_mid']            = np.nan
df_betchov_breakup_1955['x_high']           = np.nan
df_betchov_breakup_1955['photo page fig']   = np.nan
df_betchov_breakup_1955['spectrum']         = np.nan
df_betchov_breakup_1955['lighting']         = np.nan
df_betchov_breakup_1955['background color'] = np.nan
df_betchov_breakup_1955['grid']             = np.nan
df_betchov_breakup_1955['distance']         = np.nan
df_betchov_breakup_1955['f-stop']           = np.nan
df_betchov_breakup_1955['focal length']     = np.nan
df_betchov_breakup_1955['camera model']     = np.nan
df_betchov_breakup_1955['sensor']           = np.nan

df_betchov_breakup_1955 = df_betchov_breakup_1955[cols]
summary_table(df_betchov_breakup_1955)
df_jet_breakup = pd.concat([df_jet_breakup, df_betchov_breakup_1955])

##################################################
# eisenklam_flow_1958 / hooper_experimental_1959 #
##################################################

# Issues with this data:
# 1. It is unclear what the length in table 2 is, and whether it differs from the length in table 3. When I transcribed the data I thought the table 2 length was a laminar transition length (like table 3), but it appears actually that is false because some of the cases where the length is nonzero appear turbulent at the exit. I suspect this is x_i now.
# 2. Is the friction factor in table 2 uncorrected? It seems so, as they are higher than expected for laminar flow.

T = None
d_0 = None
i = None
j = None

# Liquid: water
# Ambient conditions: atmospheric and vacuum (29.0" Hg; see fig. 9)
# Orientation: down (not written but suggested by figure 1)

# TODO: Change header of spreadsheet to be standard, i.e., put units in the first row.

# TODO: Use Ubar_0, d_0, etc., to check that Re_l0 is correct for this data.
# TODO: Check whether the flow is fully developed. The shortest tubes here might not be fully developed at the highest Reynolds numbers.
# TODO: Use friction factor to determine whether data not manually categorized is turbulent or not.

# DONE: Was there any discussion of how the breakup length was measured? Fig. 13 says "(glassy jet)" which might indicate that the breakup length given was actually a transition length. I assume the length given is the surface transition location.
# DONE: Transcribe data from fig. 13 (pdf p. 46). The atmospheric pressure case is in table 2, tube 2. I assume the friction factor for the vacuum case is the same as that of the atmospheric case.
# DONE: Evaluate t/d_0 for eisenklam_flow_1958 based on the photos of the tube.
# DONE: t/d_0 photos: error analysis

eisenklam_flow_1958_csv = pd.read_csv('../data/eisenklam_flow_1958/eisenklam_flow_1958.csv', sep=',', header=0, skiprows=[1])

page_fig_eisenklam_flow_1958                     = 'p. '+str(eisenklam_flow_1958_csv['page'])+', '+eisenklam_flow_1958_csv['page fig']
tube_eisenklam_flow_1958                         = eisenklam_flow_1958_csv['tube']
d_0_eisenklam_flow_1958                          = eisenklam_flow_1958_csv['d_0'] # m
L_0s_eisenklam_flow_1958                         = eisenklam_flow_1958_csv['L_0/d_0']
Re_l0_eisenklam_flow_1958                        = eisenklam_flow_1958_csv['Re_l0']
Ubar_0_eisenklam_flow_1958                       = eisenklam_flow_1958_csv['Ubar_0'] # m/s
f_eisenklam_flow_1958                            = eisenklam_flow_1958_csv['f']
trip_eisenklam_flow_1958                         = eisenklam_flow_1958_csv['trip']
x_i_eisenklam_flow_1958                          = eisenklam_flow_1958_csv['x_i'] # m
x_trans_eisenklam_flow_1958                      = eisenklam_flow_1958_csv['x_trans'] # m
x_trans_transcription_method_eisenklam_flow_1958 = eisenklam_flow_1958_csv['transcription method']
regime_turb_eisenklam_flow_1958_2                = eisenklam_flow_1958_csv['regime turb']
regime_photo_eisenklam_flow_1958                 = eisenklam_flow_1958_csv['regime photo']
photo_page_fig_eisenklam_flow_1958               = eisenklam_flow_1958_csv['photo page fig']
photo_filename_eisenklam_flow_1958               = eisenklam_flow_1958_csv['photo filename']

# Non-dimensionalize x_trans_eisenklam_flow_1958. Need to do loop due to nans in the raw data.
i = 0
x_transs_eisenklam_flow_1958    = np.zeros(len(x_trans_eisenklam_flow_1958))
x_is_eisenklam_flow_1958        = np.zeros(len(x_trans_eisenklam_flow_1958))
est_turb_eisenklam_flow_1958    = np.zeros(len(x_trans_eisenklam_flow_1958))
I_0_eisenklam_flow_1958         = np.zeros(len(f_eisenklam_flow_1958))
Lambda_r_s_eisenklam_flow_1958  = np.zeros(len(f_eisenklam_flow_1958))
est_f_eisenklam_flow_1958       = []
trip2_eisenklam_flow_1958       = []
regime_turb_eisenklam_flow_1958 = []
sumnu = 0
N_nu  = 0
for x_trans in x_trans_eisenklam_flow_1958:
   x_i         = x_i_eisenklam_flow_1958[i]
   regime_turb = regime_turb_eisenklam_flow_1958_2[i]
   tube        = tube_eisenklam_flow_1958[i]
   Re_l0       = Re_l0_eisenklam_flow_1958[i]
   d_0         = d_0_eisenklam_flow_1958[i]
   f           = f_eisenklam_flow_1958[i]
   
   if np.isnan(x_i):
      x_is_eisenklam_flow_1958[i] = np.nan
   else:
      x_is_eisenklam_flow_1958[i] = x_i / d_0
   
   if np.isnan(x_trans):
      x_transs_eisenklam_flow_1958[i] = np.nan
   else:
      x_transs_eisenklam_flow_1958[i] = x_trans / d_0
   
   nu = (Ubar_0_eisenklam_flow_1958[i] * d_0) / Re_l0
   if not(np.isnan(nu)):
      sumnu = sumnu + nu
      
      N_nu = N_nu + 1
   
   if trip_eisenklam_flow_1958[i]:
      trip2_eisenklam_flow_1958.append(True)
   else:
      trip2_eisenklam_flow_1958.append(False)
   
   # TODO: Add transitional Re.
   # p. 6 (pdf p. 11): tubes 1 and 3 transition at about Re = 9000
   # p. 10 (pdf p. 15): tube 3 turbulent at Re = 9200
   # p. 10 (pdf p. 15): tube 4 laminar at Re = 12500 (I put this in the spreadsheet.)
   # p. 11: transition Reynolds number of about 10000?
   # p. 16 (pdf p. 21): turbulent jet in fig. 15 has Re between 3000 and 10000
   
   #if ((regime_turb != 'turbulent') and (regime_turb != 'laminar')):
   f_laminar = friction_factor_laminar(Re_l0)
   if (f / f_laminar > f_ratio_cutoff):
      regime_turb = 'turbulent'
      print Re_l0, f, f_laminar
   else:
      regime_turb = 'laminar'
   
   #if isinstance(regime_turb, basestring):
   if ((regime_turb != 'turbulent') and (regime_turb != 'laminar')):
      if tube == 1:
         Re_trans_eisenklam_flow_1958 = 9000
      elif tube == 2:
         Re_trans_eisenklam_flow_1958 = 13000
      elif tube == 3:
         Re_trans_eisenklam_flow_1958 = 9200
      else: # tube 5 and others
         Re_trans_eisenklam_flow_1958 = Re_trans
      
      if Re_l0 < Re_trans_eisenklam_flow_1958:
         regime_turb_eisenklam_flow_1958.append('laminar')
         I_0_eisenklam_flow_1958[i] = 0
         Lambda_r_s_eisenklam_flow_1958[i] = np.nan
         est_f_eisenklam_flow_1958.append(False)
      else:
         regime_turb_eisenklam_flow_1958.append('turbulent')
         if np.isnan(f):
            #print Re_l0
            I_0_eisenklam_flow_1958[i] = I_fully_developed(friction_factor_smooth(Re_l0))
            Lambda_r_s_eisenklam_flow_1958[i] = Lambda_s('smooth', friction_factor_smooth(Re_l0), 'v')
            est_f_eisenklam_flow_1958.append(True)
         else:
            I_0_eisenklam_flow_1958[i] = I_fully_developed(f)
            Lambda_r_s_eisenklam_flow_1958[i] = Lambda_s('smooth', f, 'v')
            est_f_eisenklam_flow_1958.append(False)
      
      est_turb_eisenklam_flow_1958[i] = True
   else:
      regime_turb_eisenklam_flow_1958.append(regime_turb)
      est_turb_eisenklam_flow_1958[i] = False
      
      if regime_turb == 'turbulent':
         if np.isnan(f):
            #print Re_l0
            I_0_eisenklam_flow_1958[i] = I_fully_developed(friction_factor_smooth(Re_l0))
            Lambda_r_s_eisenklam_flow_1958[i] = Lambda_s('smooth', friction_factor_smooth(Re_l0), 'v')
            est_f_eisenklam_flow_1958.append(True)
         else:
            I_0_eisenklam_flow_1958[i] = I_fully_developed(f)
            Lambda_r_s_eisenklam_flow_1958[i] = Lambda_s('smooth', f, 'v')
            est_f_eisenklam_flow_1958.append(False)
      elif regime_turb == 'laminar':
         I_0_eisenklam_flow_1958[i] = 0
         Lambda_r_s_eisenklam_flow_1958[i] = np.nan
         est_f_eisenklam_flow_1958.append(False)
      else:
         raise ValueError('Invalid regime_turb for eisenklam_flow_1958:', regime_turb)
   
   i = i + 1

trip_eisenklam_flow_1958 = trip2_eisenklam_flow_1958

nubar = sumnu / N_nu
T     = temp_from_nu_water(nubar)
rho_l = rho_water(T)
rho_g = rho_ideal_gas(P_atm, T, MW_air)
nu_l  = nu_water(T)
nu_g  = rho_g * mu_ideal_gas(T, mu_0_air, T_0_air, C_air)
sigma = sigma_water(T)

We_l0_eisenklam_flow_1958 = rho_l * Ubar_0_eisenklam_flow_1958**2 * d_0_eisenklam_flow_1958 / sigma
Fr_0_eisenklam_flow_1958  = Ubar_0_eisenklam_flow_1958**2 / (g * d_0_eisenklam_flow_1958)
#K_c_eisenklam_flow_1958   = (P_atm - P_v_water(T)) / (0.5 * rho_l * Ubar_0_eisenklam_flow_1958**2)
K_c_eisenklam_flow_1958   = K_c(P_atm, P_v_water(T), rho_l, Ubar_0_eisenklam_flow_1958, K_L_wellrounded, f_eisenklam_flow_1958, L_0s_eisenklam_flow_1958)
Ma_g_eisenklam_flow_1958  = Ubar_0 / c_ideal_gas(gamma_air(T), T + T_zero, MW_air)

df_eisenklam_flow_1958 = pd.DataFrame({'trip': trip_eisenklam_flow_1958})
df_eisenklam_flow_1958['key']        = 'eisenklam_flow_1958'
df_eisenklam_flow_1958['alt key']    = 'hooper_experimental_1959'
df_eisenklam_flow_1958['liquid']     = 'water'
df_eisenklam_flow_1958['gas']        = 'air'
df_eisenklam_flow_1958['degas']      = np.nan
df_eisenklam_flow_1958['regulation'] = True # TODO Unclear if regulated or not. Presumably the "air storage bottles" mentioned on p. 3 are regulated.

df_eisenklam_flow_1958['bend']              = False # See fig. 2.
df_eisenklam_flow_1958['flow straightener'] = False # Screen shown in fig. 2, but no flow straightener.
df_eisenklam_flow_1958['screen']            = True
df_eisenklam_flow_1958['contraction shape'] = 'smooth' # See fig. 2.
df_eisenklam_flow_1958['c']                 = 25.4e-3 / d_0_eisenklam_flow_1958 # Fig. 2 says the OD of the pipe is 1". No ID is given, so I'm using the OD. This is not expected to be a bad choice because the contraction ratio is so large that the turbulence reduction effects should be insensitive to an error of this order.
#df_eisenklam_flow_1958['trip']              = trip_eisenklam_flow_1958
df_eisenklam_flow_1958['L_r/d_0']           = 0
df_eisenklam_flow_1958['eps_r/d_0']         = np.nan

df_eisenklam_flow_1958['L_0/d_0']       = L_0s_eisenklam_flow_1958 
df_eisenklam_flow_1958['roughness']     = 'smooth'
df_eisenklam_flow_1958['f page fig']    = page_fig_eisenklam_flow_1958
df_eisenklam_flow_1958['pipe material'] = 'glass' # p. 1: > precision bore glass tubes
df_eisenklam_flow_1958['est f']         = est_f_eisenklam_flow_1958
df_eisenklam_flow_1958['check FD']      = True # See fig. 11.

df_eisenklam_flow_1958['t/d_0']       = 4.27 # Estimated based on analysis of the photos.
df_eisenklam_flow_1958['L_tip/d_0']   = np.nan
df_eisenklam_flow_1958['end checked'] = np.nan

df_eisenklam_flow_1958['orientation']         = 'down' # Based on fig. 1
df_eisenklam_flow_1958['vibration isolation'] = np.nan
df_eisenklam_flow_1958['d_chamber/d_0']       = np.nan # Never stated. A vacuum chamber was used at some point, so it's likely the jets were enclosed in some sort of chamber.

df_eisenklam_flow_1958['We_l0']            = We_l0_eisenklam_flow_1958
df_eisenklam_flow_1958['We_l0 resolution'] = np.nan
df_eisenklam_flow_1958['Re_l0']            = Re_l0_eisenklam_flow_1958
df_eisenklam_flow_1958['Re_l0 resolution'] = np.nan
df_eisenklam_flow_1958['I_0']              = I_0_eisenklam_flow_1958
df_eisenklam_flow_1958['I_0 resolution']   = np.nan
df_eisenklam_flow_1958['Fr_0']             = Fr_0_eisenklam_flow_1958
df_eisenklam_flow_1958['Fr_0 resolution']  = np.nan
df_eisenklam_flow_1958['rho_s']            = rho_l / rho_g
df_eisenklam_flow_1958['rho_s resolution'] = np.nan
df_eisenklam_flow_1958['nu_s']             = nu_l / nu_g
df_eisenklam_flow_1958['nu_s resolution']  = np.nan
df_eisenklam_flow_1958['K_c']              = K_c_eisenklam_flow_1958
df_eisenklam_flow_1958['K_c resolution']   = np.nan
df_eisenklam_flow_1958['Ma_g']             = Ma_g_eisenklam_flow_1958
df_eisenklam_flow_1958['Ma_g resolution']  = np.nan
df_eisenklam_flow_1958['Lambda_r_s']       = Lambda_r_s_eisenklam_flow_1958
df_eisenklam_flow_1958['est rho_s']        = True
df_eisenklam_flow_1958['est nu_s']         = True

df_eisenklam_flow_1958['regime photo'] = regime_photo_eisenklam_flow_1958
df_eisenklam_flow_1958['regime L_b']   = np.nan
df_eisenklam_flow_1958['regime turb']  = regime_turb_eisenklam_flow_1958
df_eisenklam_flow_1958['est turb']     = est_turb_eisenklam_flow_1958

df_eisenklam_flow_1958['L_b/d_0']                      = np.nan
df_eisenklam_flow_1958['L_b/d_0 stat error']           = np.nan
df_eisenklam_flow_1958['L_b/d_0 resolution']           = np.nan
df_eisenklam_flow_1958['L_b method']                   = np.nan
df_eisenklam_flow_1958['L_b/d_0 page fig']             = np.nan
df_eisenklam_flow_1958['L_b/d_0 transcription method'] = np.nan

df_eisenklam_flow_1958['theta']                      = np.nan
df_eisenklam_flow_1958['theta stat error']           = np.nan
df_eisenklam_flow_1958['theta resolution']           = np.nan
df_eisenklam_flow_1958['theta page fig']             = np.nan
df_eisenklam_flow_1958['theta transcription method'] = np.nan

df_eisenklam_flow_1958['D_10/d_0']                   = np.nan
df_eisenklam_flow_1958['D_10/d_0 stat error']        = np.nan
df_eisenklam_flow_1958['D_10/d_0 resolution']        = np.nan
df_eisenklam_flow_1958['D_30/d_0']                   = np.nan
df_eisenklam_flow_1958['D_30/d_0 stat error']        = np.nan
df_eisenklam_flow_1958['D_30/d_0 resolution']        = np.nan
df_eisenklam_flow_1958['D_32/d_0']                   = np.nan
df_eisenklam_flow_1958['D_32/d_0 stat error']        = np.nan
df_eisenklam_flow_1958['D_32/d_0 resolution']        = np.nan
df_eisenklam_flow_1958['droplet x/d_0']              = np.nan
df_eisenklam_flow_1958['D/d_0 page fig']             = np.nan
df_eisenklam_flow_1958['D/d_0 transcription method'] = np.nan
df_eisenklam_flow_1958['D/d_0 measurement method']   = np.nan
df_eisenklam_flow_1958['v_d_bar/vp']                        = np.nan
df_eisenklam_flow_1958['v_d_bar/vp stat error']             = np.nan
df_eisenklam_flow_1958['v_d_bar/vp resolution']             = np.nan
df_eisenklam_flow_1958['v_d_bar/vp page fig']               = np.nan
df_eisenklam_flow_1958['v_d_bar/vp transcription method']   = np.nan
df_eisenklam_flow_1958['e_p']                        = np.nan
df_eisenklam_flow_1958['e_p page fig']               = np.nan
df_eisenklam_flow_1958['e_p transcription method']   = np.nan

df_eisenklam_flow_1958['x_trans/d_0']                      = x_transs_eisenklam_flow_1958
df_eisenklam_flow_1958['x_trans/d_0 page fig']             = page_fig_eisenklam_flow_1958
df_eisenklam_flow_1958['x_trans/d_0 transcription method'] = x_trans_transcription_method_eisenklam_flow_1958

df_eisenklam_flow_1958['x_i/d_0']                      = x_is_eisenklam_flow_1958
df_eisenklam_flow_1958['x_i/d_0 stat error']           = np.nan
df_eisenklam_flow_1958['x_i/d_0 resolution']           = np.nan # TODO
df_eisenklam_flow_1958['x_i/d_0 page fig']             = np.nan
df_eisenklam_flow_1958['x_i/d_0 transcription method'] = x_trans_transcription_method_eisenklam_flow_1958

df_eisenklam_flow_1958['x_e/d_0']                      = np.nan
df_eisenklam_flow_1958['x_e/d_0 stat error']           = np.nan
df_eisenklam_flow_1958['x_e/d_0 resolution']           = np.nan
df_eisenklam_flow_1958['x_e/d_0 page fig']             = np.nan
df_eisenklam_flow_1958['x_e/d_0 transcription method'] = np.nan

# p. 3: > Photographs of the issuing jets were taken using rear diffuse illumination from a short flash duration xenon discharge tube, but in one instance shadow photographs were taken using a masked discharge tube as the point source. Further, a high speed cine camera with synchronized stroboflash was used with rear diffuse illumination.
# p. 28: > average framing rate 3750 per second
# fig. 8 (pdf p. 41): > interval between frames 1/3750 sec.
df_eisenklam_flow_1958['photo filename']   = photo_filename_eisenklam_flow_1958
df_eisenklam_flow_1958['exposure time']    = 1/3750
df_eisenklam_flow_1958['flash']            = True
df_eisenklam_flow_1958['x_low']            = 0
df_eisenklam_flow_1958['x_mid']            = np.nan
df_eisenklam_flow_1958['x_high']           = np.nan
df_eisenklam_flow_1958['photo page fig']   = photo_page_fig_eisenklam_flow_1958
df_eisenklam_flow_1958['spectrum']         = 'visible'
df_eisenklam_flow_1958['lighting']         = 'diffuse background' # p. 3
df_eisenklam_flow_1958['background color'] = 'light'
df_eisenklam_flow_1958['grid']             = False
df_eisenklam_flow_1958['distance']         = 3*12*2.54e-2 # m, from fig. 3
df_eisenklam_flow_1958['f-stop']           = np.nan
df_eisenklam_flow_1958['focal length']     = np.nan
df_eisenklam_flow_1958['camera model']     = np.nan
df_eisenklam_flow_1958['sensor']           = 'film'

df_eisenklam_flow_1958 = df_eisenklam_flow_1958[cols]
summary_table(df_eisenklam_flow_1958)
df_jet_breakup = pd.concat([df_jet_breakup, df_eisenklam_flow_1958])

#######################
# chen_mechanics_1962 #
#######################

# TODO: Add transitional regimes (R2F, F2S, S2A).

T = None
d_0 = None
i = None
j = None

# DONE: Add photos for this.
# DONE: Determine regime from photos.
# TODO: Check if errata for Chen and Davis influenced anything in your script.

# Issues with this data:
# 1. chen_disintegration_1964 figure 3 and figure 6 are inconsistent. For example, the Re_l0 = 104,000, d_0 = 1/4 inch case does not appear to be labeled correctly in figure 6. Consequently, I obtained a copy of the dissertation, which had the data tabulated.
# 2. In the friction factor plot (fig. 23), 3/4" is listed twice. I assume the first was supposed to be 3/8" as the sizes appear to be in increasing order.
# 3. The computed Reynolds number is frequently inconsistent with the provided raw data. I assume the raw data is correct. The Weber number is consistent with the data. I assume the Reynolds number in the friction factor plot (fig. 23) is correct, though it might not be. The assumption that the friction factor plot Reynolds number is correct seems supported by the fact that many of data points in that plot appear to match the true Reynolds number of many of the tests. But many also do not. Perhaps some of the Reynolds numbers are correct and others are not. Either way, there is not much variation, so the correlation is expected to be reasonably accurate.
# 4. The photos in chen_mechanics_1962 p. 71 fig. 20 are completely different from the corresponding photos in chen_disintegration_1964 p. 182 fig. 3. Given that some of the photos in chen_mechanics_1962 appear inconsistent with the regime based on the breakup length and that chen_disintegration_1964 was published later, I assume that chen_disintegration_1964 was corrected.

# p. 31 \S c.: > maximum differences of 1:1.72 in the viscosity, and 1:1.52 in the surface tension were obtained

chen_mechanics_1962_f_csv = pd.read_csv('../data/chen_mechanics_1962/chen_mechanics_1962_fig_23.csv', sep=',', header=0)

# pp. 27-28: > In Figure 23 the comparative smoothness of the test pipe is shown indicating that the pipes were not perfectly smooth

d_0_chen_mechanics_1962   = chen_mechanics_1962_f_csv['d_0'] * 0.0254 # m
Re_l0_chen_mechanics_1962 = chen_mechanics_1962_f_csv['Re_l0']
f_chen_mechanics_1962     = chen_mechanics_1962_f_csv['f']

# Use a power law for the friction factors of the different pipes as a function of Re.

d_0_array = []
for d_0 in d_0_chen_mechanics_1962:
   d_0 = d_0 / 0.0254
   if not(d_0 in d_0_array):
      d_0_array.append(d_0)

preexponential = np.zeros(len(d_0_array))
exponent       = np.zeros(len(d_0_array))
i = 0
for d_0 in d_0_array:
   chen_mechanics_1962_f_df = chen_mechanics_1962_f_csv[chen_mechanics_1962_f_csv.d_0 == d_0]
   
   d_0 = d_0 * 0.254 # m
   Re_l0_chen_mechanics_1962 = chen_mechanics_1962_f_df['Re_l0']
   f_chen_mechanics_1962     = chen_mechanics_1962_f_df['f']
   
   coeffs = np.polyfit(log(Re_l0_chen_mechanics_1962), log(f_chen_mechanics_1962), 1)
   
   preexponential[i] = exp(coeffs[1])
   exponent[i]       = coeffs[0]
   
   #Re_l0_arr = np.logspace(4, 6, 100)
   #plt.loglog(Re_l0_chen_mechanics_1962, f_chen_mechanics_1962, 'x', Re_l0_arr, exp(coeffs[1]) * Re_l0_arr**coeffs[0], '-')
   #plt.xlabel('Re_{l0}')
   #plt.ylabel('f')
   #plt.grid()
   #plt.show()
   
   i = i + 1

d_0_array = np.array(d_0_array) * 0.3048 / 12

#df_chen_mechanics_1962_f = pd.DataFrame({'d_0': d_0_array,
                                    #'preexponential': preexponential,
                                    #'exponent': exponent})

chen_mechanics_1962_csv = pd.read_csv('../data/chen_mechanics_1962/chen_mechanics_1962.csv', sep=',', header=0)

d_0_chen_mechanics_1962            = chen_mechanics_1962_csv['d_0'] * 0.3048 / 12 # m
Ubar_0_chen_mechanics_1962         = chen_mechanics_1962_csv['Ubar_0'] * 0.3048 # m/s
L_bs_chen_mechanics_1962           = chen_mechanics_1962_csv['L_b/d_0']
n_chen_mechanics_1962              = chen_mechanics_1962_csv['n']
mu_l_chen_mechanics_1962           = chen_mechanics_1962_csv['mu_l'] * 47.880259 # Pa*s
rho_l_chen_mechanics_1962          = chen_mechanics_1962_csv['rho_l'] * 515.379 # kg/m^3
sigma_chen_mechanics_1962          = chen_mechanics_1962_csv['sigma'] * 14.5939 # N/m
fluid_chen_mechanics_1962          = chen_mechanics_1962_csv['liquid']
D_10_chen_mechanics_1962           = chen_mechanics_1962_csv['D_10'] # Seems to be initial drop size near the breakup point. The earlier spray formation is not a part of this. This could be useful for the range theory as the largest drops are important for range, and the drops in this case are those formed immediately after complete breakup. See chen_disintegration_1964 p. 181: >The photographs taken for [...]
regime_L_b_chen_mechanics_1962     = chen_mechanics_1962_csv['regime']
photo_filename_chen_mechanics_1962 = chen_mechanics_1962_csv['photo filename']
photo_page_fig_chen_mechanics_1962 = chen_mechanics_1962_csv['photo page fig']

rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g

i = 0
f_chen_mechanics_1962           = np.zeros(len(d_0_chen_mechanics_1962))
Re_l0_chen_mechanics_1962       = np.zeros(len(d_0_chen_mechanics_1962))
We_l0_chen_mechanics_1962       = np.zeros(len(d_0_chen_mechanics_1962))
Fr_0_chen_mechanics_1962        = np.zeros(len(d_0_chen_mechanics_1962))
rho_s_chen_mechanics_1962       = np.zeros(len(d_0_chen_mechanics_1962))
I_0_chen_mechanics_1962         = np.zeros(len(d_0_chen_mechanics_1962))
nu_s_chen_mechanics_1962        = np.zeros(len(d_0_chen_mechanics_1962))
K_c_chen_mechanics_1962         = np.zeros(len(d_0_chen_mechanics_1962))
Ma_g_chen_mechanics_1962        = np.zeros(len(d_0_chen_mechanics_1962))
Lambda_r_s_chen_mechanics_1962  = np.zeros(len(d_0_chen_mechanics_1962))
regime_turb_chen_mechanics_1962 = []
for d_0 in d_0_chen_mechanics_1962:
   rho_l  = rho_l_chen_mechanics_1962[i]
   Ubar_0 = Ubar_0_chen_mechanics_1962[i]
   mu_l   = mu_l_chen_mechanics_1962[i]
   
   Re_l0_chen_mechanics_1962[i] = rho_l * Ubar_0 * d_0 / mu_l
   
   j = np.argmin(abs(d_0_array - d_0))
   f_chen_mechanics_1962[i] = preexponential[j]*Re_l0_chen_mechanics_1962[i]**exponent[j]
   f = f_chen_mechanics_1962[i]
   
   f_laminar = friction_factor_laminar(Re_l0_chen_mechanics_1962[i])
   if (f / f_laminar > f_ratio_cutoff):
      regime_turb_chen_mechanics_1962.append('turbulent')
   else:
      regime_turb_chen_mechanics_1962.append('laminar')
   
   We_l0_chen_mechanics_1962[i]      = rho_l * Ubar_0**2 * d_0 / sigma
   Fr_0_chen_mechanics_1962[i]       = Ubar_0**2 / (g * d_0)
   I_0_chen_mechanics_1962[i]        = I_fully_developed(f)
   Lambda_r_s_chen_mechanics_1962[i] = Lambda_s('smooth', f, 'v')
   rho_s_chen_mechanics_1962[i]      = rho_l / rho_g
   nu_s_chen_mechanics_1962[i]       = nu_l / nu_g
   #K_c_chen_mechanics_1962[i]        = (P_atm - P_v_water(T_std)) / (0.5 * rho_l * Ubar_0**2)
   K_c_chen_mechanics_1962[i]        = K_c(P_atm, P_v_water(T_std), rho_l, Ubar_0, K_L_wellrounded, f, 100.)
   Ma_g_chen_mechanics_1962[i]       = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)
   
   #print j, Re_l0_chen_mechanics_1962[i], sqrt(We_l0_chen_mechanics_1962[i]), f_chen_mechanics_1962[i], I_0_chen_mechanics_1962[i]
   
   i = i + 1

# TODO: Error analysis
# fig. 9, p. 63 (pdf p. 73): standard deviation of L_b as a function of Re_l0
# tab. 1, p. 68 (pdf p. 78): N photos for L_b/d_0 determination
# tab. 5, p. 76 (pdf p. 86): droplet size standard deviation

df_chen_mechanics_1962 = pd.DataFrame({'key': 'chen_mechanics_1962',
                                       'liquid': fluid_chen_mechanics_1962})
df_chen_mechanics_1962['alt key']    = 'chen_disintegration_1964'
df_chen_mechanics_1962['gas']        = 'air'
df_chen_mechanics_1962['degas']      = np.nan
df_chen_mechanics_1962['regulation'] = 'pump with bypass' # See pp. 30-31, 33.

df_chen_mechanics_1962['bend']              = False # Guess based on fig. 11.
df_chen_mechanics_1962['flow straightener'] = False # Guess based on fig. 11.
df_chen_mechanics_1962['screen']            = False # Guess based on fig. 11.
df_chen_mechanics_1962['contraction shape'] = 'smooth' # p. 32: > A brass, bell-mouthed entrance block was used for the connection of the test pipe to the 8 inch I.D. approaching section as shown in Figure 11.
df_chen_mechanics_1962['c']                 = 8. / chen_mechanics_1962_csv['d_0'] # See bottom of p. 28 and also p. 32. 8 inch diameter inlet pipe.
df_chen_mechanics_1962['trip']              = False
df_chen_mechanics_1962['L_r/d_0']           = 0.
df_chen_mechanics_1962['eps_r/d_0']         = np.nan

df_chen_mechanics_1962['L_0/d_0']       = 100. # p. 16: "straight smooth glass tube of length 100 times the diameter"; p. 26: > straight, smooth pipes of different sizes and 100 diameters in length
df_chen_mechanics_1962['roughness']     = 'rough'
df_chen_mechanics_1962['f page fig']    = 'p. 75, fig. 23'
df_chen_mechanics_1962['pipe material'] = 'PVC pipe' # p. 31: > Four sizes of smooth Kraloy rigid unplasticized polyvinyl chloride pipe
df_chen_mechanics_1962['est f']         = False
df_chen_mechanics_1962['check FD']      = True # p. 26: > A part of this study thus involved measurement of the pressure distribution along the longitudinal axis of the pipe in order to insure that fulyl developed turbulent flow existed prior to the discharge from the pipe.
# p. 33: > This was confirmed through indirect observation of the pressure fluctuations, hydraulic gradient along the test pipe and also by periodic gravimetric measurement of the flow.

df_chen_mechanics_1962['t/d_0']       = np.nan # TODO Check photo 4.
df_chen_mechanics_1962['L_tip/d_0']   = np.nan
df_chen_mechanics_1962['end checked'] = True # p. 32: > To secure a sharp edge at the exit of the jet and to minimize flow irregularities, the test pipes were cut cleanly and sharply with great care. The cut was then cleaned with soft cloth and examined with a microscope for any irregularities or protuberances.

df_chen_mechanics_1962['orientation']         = 'horizontal' # chen_measuring_1961 seems to indicate the experiments were horizontal
df_chen_mechanics_1962['vibration isolation'] = True # p. 31: > The test pipe was supported by a metal frame which was separated from the rest of the system in order to keep it as free as possible from external vibrations.
df_chen_mechanics_1962['d_chamber/d_0']       = np.inf

df_chen_mechanics_1962['We_l0']      = We_l0_chen_mechanics_1962
df_chen_mechanics_1962['Re_l0']      = Re_l0_chen_mechanics_1962
df_chen_mechanics_1962['I_0']        = I_0_chen_mechanics_1962
df_chen_mechanics_1962['Fr_0']       = Fr_0_chen_mechanics_1962
df_chen_mechanics_1962['rho_s']      = rho_s_chen_mechanics_1962
df_chen_mechanics_1962['nu_s']       = nu_s_chen_mechanics_1962
df_chen_mechanics_1962['K_c']        = K_c_chen_mechanics_1962
df_chen_mechanics_1962['Ma_g']       = Ma_g_chen_mechanics_1962
df_chen_mechanics_1962['Lambda_r_s'] = Lambda_r_s_chen_mechanics_1962
df_chen_mechanics_1962['est rho_s']  = True
df_chen_mechanics_1962['est nu_s']   = True

df_chen_mechanics_1962['regime photo'] = np.nan
df_chen_mechanics_1962['regime L_b']   = regime_L_b_chen_mechanics_1962
df_chen_mechanics_1962['regime turb']  = regime_turb_chen_mechanics_1962
df_chen_mechanics_1962['est turb']     = False

df_chen_mechanics_1962['L_b/d_0']                      = L_bs_chen_mechanics_1962
df_chen_mechanics_1962['L_b/d_0 stat error']           = scipy.stats.t.ppf(1 - (1 - interval_probability_level) / 2, n_chen_mechanics_1962 - 1) * breakup_length_sigmas * L_bs_chen_mechanics_1962 / sqrt(n_chen_mechanics_1962) # Using phinney_breakup_1973's estimate for the standard deviation. TODO: Use fig. 9 instead.
df_chen_mechanics_1962['L_b/d_0 resolution']           = 0.5*2.54e-2 / d_0_chen_mechanics_1962 # See p. 30: > The reading was taken with an accuracy of \pm 0.5''
df_chen_mechanics_1962['L_b method']                   = 'photographic'
df_chen_mechanics_1962['L_b/d_0 page fig']             = 'p. 68, tab. 1; p. 75. fig. 23., p. 76. tab. 5' # TODO
df_chen_mechanics_1962['L_b/d_0 transcription method'] = 'table'

df_chen_mechanics_1962['theta']                      = np.nan
df_chen_mechanics_1962['theta stat error']           = np.nan
df_chen_mechanics_1962['theta resolution']           = np.nan
df_chen_mechanics_1962['theta page fig']             = np.nan
df_chen_mechanics_1962['theta transcription method'] = np.nan

df_chen_mechanics_1962['D_10/d_0']                   = D_10_chen_mechanics_1962 / d_0_chen_mechanics_1962
df_chen_mechanics_1962['D_10/d_0 stat error']        = np.nan # TODO
df_chen_mechanics_1962['D_10/d_0 resolution']        = np.nan
df_chen_mechanics_1962['D_32/d_0']                   = np.nan
df_chen_mechanics_1962['D_32/d_0 stat error']        = np.nan
df_chen_mechanics_1962['D_32/d_0 resolution']        = np.nan
df_chen_mechanics_1962['D_30/d_0']                   = np.nan
df_chen_mechanics_1962['D_30/d_0 stat error']        = np.nan
df_chen_mechanics_1962['D_30/d_0 resolution']        = np.nan
df_chen_mechanics_1962['droplet x/d_0']              = np.nan # TODO: Add location for D_{10} from chen_mechanics_1962/chen_disintegration_1964
df_chen_mechanics_1962['D/d_0 page fig']             = np.nan
df_chen_mechanics_1962['D/d_0 transcription method'] = 'table'
df_chen_mechanics_1962['D/d_0 measurement method']   = 'photographic'
df_chen_mechanics_1962['v_d_bar/vp']                        = np.nan
df_chen_mechanics_1962['v_d_bar/vp stat error']             = np.nan
df_chen_mechanics_1962['v_d_bar/vp resolution']             = np.nan
df_chen_mechanics_1962['v_d_bar/vp page fig']               = np.nan
df_chen_mechanics_1962['v_d_bar/vp transcription method']   = np.nan
df_chen_mechanics_1962['e_p']                        = np.nan
df_chen_mechanics_1962['e_p page fig']               = np.nan
df_chen_mechanics_1962['e_p transcription method']   = np.nan

df_chen_mechanics_1962['x_trans/d_0']                      = np.nan
df_chen_mechanics_1962['x_trans/d_0 page fig']             = np.nan
df_chen_mechanics_1962['x_trans/d_0 transcription method'] = np.nan

df_chen_mechanics_1962['x_i/d_0']                      = np.nan
df_chen_mechanics_1962['x_i/d_0 stat error']           = np.nan
df_chen_mechanics_1962['x_i/d_0 resolution']           = np.nan
df_chen_mechanics_1962['x_i/d_0 page fig']             = np.nan
df_chen_mechanics_1962['x_i/d_0 transcription method'] = np.nan

df_chen_mechanics_1962['x_e/d_0']                      = np.nan
df_chen_mechanics_1962['x_e/d_0 stat error']           = np.nan
df_chen_mechanics_1962['x_e/d_0 resolution']           = np.nan
df_chen_mechanics_1962['x_e/d_0 page fig']             = np.nan
df_chen_mechanics_1962['x_e/d_0 transcription method'] = np.nan

df_chen_mechanics_1962['photo filename']   = photo_filename_chen_mechanics_1962
df_chen_mechanics_1962['exposure time']    = 2e-6  # p. 29: > A high intensity photoflash light of 2 microsecond duration was used
df_chen_mechanics_1962['flash']            = True # p. 29: > A high intensity photoflash light of 2 microsecond duration was used
df_chen_mechanics_1962['x_low']            = np.nan
df_chen_mechanics_1962['x_mid']            = L_bs_chen_mechanics_1962
df_chen_mechanics_1962['x_high']           = np.nan
df_chen_mechanics_1962['photo page fig']   = photo_page_fig_chen_mechanics_1962
df_chen_mechanics_1962['spectrum']         = 'visible'
df_chen_mechanics_1962['lighting']         = 'diffuse background'
df_chen_mechanics_1962['background color'] = 'black' # p. 30
df_chen_mechanics_1962['grid']             = True # 2 inch square grid according to p. 29. Fig. 4 has a ruler.
df_chen_mechanics_1962['distance']         = np.nan
df_chen_mechanics_1962['f-stop']           = 2.8
df_chen_mechanics_1962['focal length']     = 50 # mm
df_chen_mechanics_1962['camera model']     = 'Exa with Tessar lens'
df_chen_mechanics_1962['sensor']           = '35 film, ISO 250 to 350'

df_chen_mechanics_1962 = df_chen_mechanics_1962[cols]
summary_table(df_chen_mechanics_1962)
df_jet_breakup = pd.concat([df_jet_breakup, df_chen_mechanics_1962])

#####################
# rupe_dynamic_1962 #
#####################

# TODO: Add transitional regimes (R2F, F2S, S2A).

T = None
d_0 = None
i = None
j = None

# Issues with this data:
# 1. Page 11R says the glycerin-water mix is 79% by weight glycerin, but figure 4 says 77% by weight glycerin. I'll use the properties listed in the text, and assume 79% by weight glycerin for the surface tension (which is not given).

rupe_dynamic_1962_csv = pd.read_csv('../data/rupe_dynamic_1962/rupe_dynamic_1962.csv', sep=',', header=0)

liquid_rupe_dynamic_1962         = rupe_dynamic_1962_csv['liquid']
page_fig_rupe_dynamic_1962       = rupe_dynamic_1962_csv['page fig']
d_0_rupe_dynamic_1962            = rupe_dynamic_1962_csv['d_0 (in)'] * 25.4e-3 # m
L_0s_rupe_dynamic_1962           = rupe_dynamic_1962_csv['L_0/d_0']
Ubar_0_rupe_dynamic_1962         = rupe_dynamic_1962_csv['Ubar_0 (ft/s)'] * 0.3048 # m/s
Re_l0_rupe_dynamic_1962          = rupe_dynamic_1962_csv['Re_l0']
trip_raw_rupe_dynamic_1962       = rupe_dynamic_1962_csv['trip']
L_rs_rupe_dynamic_1962           = rupe_dynamic_1962_csv['L_r/d_0']
eps_rs_rupe_dynamic_1962         = rupe_dynamic_1962_csv['eps_r/d_0']
regime_turb_rupe_dynamic_1962    = rupe_dynamic_1962_csv['regime turb']
regime_photo_rupe_dynamic_1962   = rupe_dynamic_1962_csv['regime photo']
xtranss_rupe_dynamic_1962        = rupe_dynamic_1962_csv['x_trans/d_0']
photo_filename_rupe_dynamic_1962 = rupe_dynamic_1962_csv['photo filename']

# Non-dimensionalize x_trans_eisenklam_flow_1958. Need to do loop due to nans in the raw data.
i = 0
T_rupe_dynamic_1962 = np.zeros(len(liquid_rupe_dynamic_1962))
trip_rupe_dynamic_1962 = []
sumnu = 0
N_nu  = 0
for fluid in liquid_rupe_dynamic_1962:
   if not(np.isnan(Ubar_0_rupe_dynamic_1962[i])):
      if fluid == 'water':
         nu = (Ubar_0_rupe_dynamic_1962[i] * d_0_rupe_dynamic_1962[i]) / Re_l0_rupe_dynamic_1962[i]
         sumnu = sumnu + nu
         N_nu = N_nu + 1
   
   if trip_raw_rupe_dynamic_1962[i]:
      trip_rupe_dynamic_1962.append(True)
   else:
      trip_rupe_dynamic_1962.append(False)
   
   i = i + 1

nubar   = sumnu / N_nu
T_water = temp_from_nu_water(nubar)

i = 0
xtranss_page_fig_rupe_dynamic_1962 = []
We_l0_rupe_dynamic_1962      = np.zeros(len(liquid_rupe_dynamic_1962))
Re_l0_rupe_dynamic_1962      = np.zeros(len(liquid_rupe_dynamic_1962))
Fr_0_rupe_dynamic_1962       = np.zeros(len(liquid_rupe_dynamic_1962))
rho_s_rupe_dynamic_1962      = np.zeros(len(liquid_rupe_dynamic_1962))
nu_s_rupe_dynamic_1962       = np.zeros(len(liquid_rupe_dynamic_1962))
K_c_rupe_dynamic_1962        = np.zeros(len(liquid_rupe_dynamic_1962))
Ma_g_rupe_dynamic_1962       = np.zeros(len(liquid_rupe_dynamic_1962))
for fluid in liquid_rupe_dynamic_1962:
   if fluid == 'water':
      rho_l = rho_water(T_water)
      rho_g = rho_ideal_gas(P_atm, T_water, MW_air)
      nu_l  = nu_water(T_water)
      nu_g  = rho_g * mu_ideal_gas(T_water, mu_0_air, T_0_air, C_air)
      sigma = sigma_water(T_water)
      P_v   = P_v_water(T_water)
   elif fluid == '77% by weight glycerin in water':
      # Surface tension, vapor pressure, and boiling point interpolated from glycerine_producers_association_physical_1975.
      # glycerine_producers_association_physical_1975 p. 11R has physical properties of glycerin-water mixture
      T             = 28.5 # C
      pctglycerine  = 79 # As for why this is not 77% as written above, see issue #1 with this data.
      
      rho_l = 1.2060e3 # kg/m^3
      mu_l  = 37.2e-3 # Pa*s, p. 11R, DONE: Where did this come from? Rupe p. 6R gives a source but not the number.
      nu_l  = mu_l / rho_l # m^2/s
      sigma = np.interp(pctglycerine, [61.44, 81.98],
                         [np.interp(T, np.array([20, 30, 40]), np.array([67.64, 66.68, 65.71]) * 0.001),
                          np.interp(T, np.array([17, 20, 30, 40]), np.array([65.41, 65.26, 64.66, 63.93])) * 0.001]) # N/m
      P_v   = np.interp(pctglycerine, np.array([11.69, 11.97, 21.63, 45.02, 61.15, 77.38, 88.03]), np.array([230.9, 229.3, 222.3, 199.9, 173.4, 128.1, 79.8])) * 133.322 # Pa, from physical_properties_of_glycerine_and_its_solutions p. 8, tab. 12
      
      #T_boil_glycerine = np.interp(pctglycerine, np.array([70, 80]), np.array([114.0, 121.5])) # C, p. 6, tab. 10
      # This book also has c_p for glycerine-water solution on p. 16 in table 26.
      # h_fg? Calculate via Clausius-Claperyon?
   else:
      raise ValueError('Invalid fluid for rupe_dynamic_1962', fluid)
   
   Ubar_0 = Ubar_0_rupe_dynamic_1962[i]
   d_0    = d_0_rupe_dynamic_1962[i]
   
   We_l0_rupe_dynamic_1962[i] = (rho_l * Ubar_0**2 * d_0) / sigma
   Re_l0_rupe_dynamic_1962[i] = Ubar_0 * d_0 / nu_l
   Fr_0_rupe_dynamic_1962[i]  = Ubar_0**2 / (g * d_0)
   rho_s_rupe_dynamic_1962[i] = rho_l / rho_g
   nu_s_rupe_dynamic_1962[i]  = nu_l / nu_g
   #K_c_rupe_dynamic_1962[i]   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_rupe_dynamic_1962[i]   = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_wellrounded, friction_factor_smooth(Re_l0_rupe_dynamic_1962[i]), L_0s_rupe_dynamic_1962[i])
   Ma_g_rupe_dynamic_1962[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_water), T_water + T_zero, MW_air)
   
   if not(np.isnan(xtranss_rupe_dynamic_1962[i])):
      xtranss_page_fig_rupe_dynamic_1962.append(page_fig_rupe_dynamic_1962[i])
   else:
      xtranss_page_fig_rupe_dynamic_1962.append(None)
   
   i = i + 1

I_0_rupe_dynamic_1962        = I_fully_developed_array(Re_l0_rupe_dynamic_1962, Re_turb)
Lambda_r_s_rupe_dynamic_1962 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_rupe_dynamic_1962, Re_turb), 'v')

df_rupe_dynamic_1962 = pd.DataFrame({'key': 'rupe_dynamic_1962',
                                     'liquid': liquid_rupe_dynamic_1962})
df_rupe_dynamic_1962['alt key']    = np.nan
df_rupe_dynamic_1962['gas']        = 'air'
df_rupe_dynamic_1962['degas']      = np.nan
df_rupe_dynamic_1962['regulation'] = True # p. 6: "conventional gas-pressurized system". Not sure what this means, but I assume it is regulated.

df_rupe_dynamic_1962['bend']              = False
df_rupe_dynamic_1962['flow straightener'] = False
df_rupe_dynamic_1962['contraction shape'] = 'smooth' # See p. 7R, fig. 2.
df_rupe_dynamic_1962['screen']            = True # pp. 6R-7L
df_rupe_dynamic_1962['c']                 = np.nan # TODO Estimate from fig. 2?
df_rupe_dynamic_1962['trip']              = trip_rupe_dynamic_1962
df_rupe_dynamic_1962['L_r/d_0']           = L_rs_rupe_dynamic_1962 # p. 8L
df_rupe_dynamic_1962['eps_r/d_0']         = eps_rs_rupe_dynamic_1962

df_rupe_dynamic_1962['L_0/d_0']       = L_0s_rupe_dynamic_1962
df_rupe_dynamic_1962['roughness']     = 'smooth'
df_rupe_dynamic_1962['f page fig']    = np.nan
df_rupe_dynamic_1962['pipe material'] = np.nan # Never mentioned?
df_rupe_dynamic_1962['est f']         = True
df_rupe_dynamic_1962['check FD']      = True

df_rupe_dynamic_1962['t/d_0']       = (3/16 - rupe_dynamic_1962_csv['d_0 (in)']) / (2 * rupe_dynamic_1962_csv['d_0 (in)']) # fig. 2, p. 7R has OD of tubes
df_rupe_dynamic_1962['L_tip/d_0']   = np.nan
df_rupe_dynamic_1962['end checked'] = True # While not specifically said, the tubes were "honed".

df_rupe_dynamic_1962['orientation']         = 'down' # See p. 6, fig. 1
df_rupe_dynamic_1962['vibration isolation'] = np.nan
df_rupe_dynamic_1962['d_chamber/d_0']       = np.nan # TODO Estimate from fig. 1?

df_rupe_dynamic_1962['We_l0']      = We_l0_rupe_dynamic_1962
df_rupe_dynamic_1962['Re_l0']      = Re_l0_rupe_dynamic_1962
df_rupe_dynamic_1962['I_0']        = I_0_rupe_dynamic_1962
df_rupe_dynamic_1962['Fr_0']       = Fr_0_rupe_dynamic_1962
df_rupe_dynamic_1962['rho_s']      = rho_s_rupe_dynamic_1962
df_rupe_dynamic_1962['nu_s']       = nu_s_rupe_dynamic_1962
df_rupe_dynamic_1962['K_c']        = K_c_rupe_dynamic_1962
df_rupe_dynamic_1962['Ma_g']       = Ma_g_rupe_dynamic_1962
df_rupe_dynamic_1962['Lambda_r_s'] = Lambda_r_s_rupe_dynamic_1962
df_rupe_dynamic_1962['est rho_s']  = True
df_rupe_dynamic_1962['est nu_s']   = True

df_rupe_dynamic_1962['regime photo']                 = regime_photo_rupe_dynamic_1962
df_rupe_dynamic_1962['regime L_b']                   = np.nan
df_rupe_dynamic_1962['regime turb']                  = regime_turb_rupe_dynamic_1962
df_rupe_dynamic_1962['est turb']                     = False

df_rupe_dynamic_1962['L_b/d_0']                      = np.nan
df_rupe_dynamic_1962['L_b/d_0 stat error']           = np.nan
df_rupe_dynamic_1962['L_b/d_0 resolution']           = np.nan
df_rupe_dynamic_1962['L_b method']                   = np.nan
df_rupe_dynamic_1962['L_b/d_0 page fig']             = np.nan
df_rupe_dynamic_1962['L_b/d_0 transcription method'] = np.nan

df_rupe_dynamic_1962['theta']                      = np.nan
df_rupe_dynamic_1962['theta stat error']           = np.nan
df_rupe_dynamic_1962['theta resolution']           = np.nan
df_rupe_dynamic_1962['theta page fig']             = np.nan
df_rupe_dynamic_1962['theta transcription method'] = np.nan

df_rupe_dynamic_1962['D_10/d_0']                   = np.nan
df_rupe_dynamic_1962['D_10/d_0 stat error']        = np.nan
df_rupe_dynamic_1962['D_10/d_0 resolution']        = np.nan
df_rupe_dynamic_1962['D_30/d_0']                   = np.nan
df_rupe_dynamic_1962['D_30/d_0 stat error']        = np.nan
df_rupe_dynamic_1962['D_30/d_0 resolution']        = np.nan
df_rupe_dynamic_1962['D_32/d_0']                   = np.nan
df_rupe_dynamic_1962['D_32/d_0 stat error']        = np.nan
df_rupe_dynamic_1962['D_32/d_0 resolution']        = np.nan
df_rupe_dynamic_1962['droplet x/d_0']              = np.nan
df_rupe_dynamic_1962['D/d_0 page fig']             = np.nan
df_rupe_dynamic_1962['D/d_0 transcription method'] = np.nan
df_rupe_dynamic_1962['D/d_0 measurement method']   = np.nan
df_rupe_dynamic_1962['v_d_bar/vp']                        = np.nan
df_rupe_dynamic_1962['v_d_bar/vp stat error']             = np.nan
df_rupe_dynamic_1962['v_d_bar/vp resolution']             = np.nan
df_rupe_dynamic_1962['v_d_bar/vp page fig']               = np.nan
df_rupe_dynamic_1962['v_d_bar/vp transcription method']   = np.nan
df_rupe_dynamic_1962['e_p']                        = np.nan
df_rupe_dynamic_1962['e_p page fig']               = np.nan
df_rupe_dynamic_1962['e_p transcription method']   = np.nan

df_rupe_dynamic_1962['x_trans/d_0']                      = xtranss_rupe_dynamic_1962
df_rupe_dynamic_1962['x_trans/d_0 page fig']             = xtranss_page_fig_rupe_dynamic_1962
df_rupe_dynamic_1962['x_trans/d_0 transcription method'] = 'photo'

df_rupe_dynamic_1962['x_i/d_0']                      = np.nan
df_rupe_dynamic_1962['x_i/d_0 stat error']           = np.nan
df_rupe_dynamic_1962['x_i/d_0 resolution']           = np.nan
df_rupe_dynamic_1962['x_i/d_0 page fig']             = np.nan
df_rupe_dynamic_1962['x_i/d_0 transcription method'] = np.nan

df_rupe_dynamic_1962['x_e/d_0']                      = np.nan
df_rupe_dynamic_1962['x_e/d_0 stat error']           = np.nan
df_rupe_dynamic_1962['x_e/d_0 resolution']           = np.nan
df_rupe_dynamic_1962['x_e/d_0 page fig']             = np.nan
df_rupe_dynamic_1962['x_e/d_0 transcription method'] = np.nan

df_rupe_dynamic_1962['photo filename'] = photo_filename_rupe_dynamic_1962
df_rupe_dynamic_1962['exposure time']  = 2e-6 # p. 7L: > effective flash time of approximately 2.0 $\mu$sec.
df_rupe_dynamic_1962['flash']            = True # See line above.
df_rupe_dynamic_1962['x_low']            = 0
df_rupe_dynamic_1962['x_mid']            = np.nan
df_rupe_dynamic_1962['x_high']           = np.nan
df_rupe_dynamic_1962['photo page fig']   = page_fig_rupe_dynamic_1962
df_rupe_dynamic_1962['spectrum']         = 'visible'
df_rupe_dynamic_1962['lighting']         = 'diffuse background and front' # p. 7L: > The jet photos were obtained with three of the tubes positioned behind a ground glass to the rear of the jet and one just above the camera lens so as to provide some front lighting
df_rupe_dynamic_1962['background color'] = 'light' # glass?
df_rupe_dynamic_1962['grid']             = False
df_rupe_dynamic_1962['distance']         = np.nan
df_rupe_dynamic_1962['f-stop']           = np.nan
df_rupe_dynamic_1962['focal length']     = np.nan
df_rupe_dynamic_1962['camera model']     = np.nan
df_rupe_dynamic_1962['sensor']           = np.nan

df_rupe_dynamic_1962 = df_rupe_dynamic_1962[cols]
summary_table(df_rupe_dynamic_1962)
df_jet_breakup = pd.concat([df_jet_breakup, df_rupe_dynamic_1962])

###############################
# skrebkov_turbulentnyye_1963 #
###############################

T = None
d_0 = None
i = None
j = None

skrebkov_turbulentnyye_1963_csv = pd.read_csv('../data/skrebkov_turbulentnyye_1963/skrebkov_turbulentnyye_1963_fig_1.csv', sep=',', header=0)

# Based on fig. 2, the spray angle given in this paper is the full angle.

vp_skrebkov_turbulentnyye_1963    = skrebkov_turbulentnyye_1963_csv['vp']
theta_skrebkov_turbulentnyye_1963 = skrebkov_turbulentnyye_1963_csv['theta'] * (2 * np.pi / 360)

Ubar_0 = 45 # m/s (p. 143 of translation)
d_0    = 7e-3 # m (p. 143 of translation)
T      = 25 # C (guess; no temperatures are given in the paper)
f_skrebkov_turbulentnyye_1963 = 8 * (vp_skrebkov_turbulentnyye_1963 / (1.2 * Ubar_0))**2

rho_l = rho_water(T)
mu_l  = mu_water(T)
sigma = sigma_water(T)
rho_g = rho_ideal_gas(P_atm, T, MW_air)
nu_g  = mu_ideal_gas(T, mu_0_air, T_0_air, C_air) / rho_g

i = 0
Re_l0_skrebkov_turbulentnyye_1963      = np.zeros(len(f_skrebkov_turbulentnyye_1963))
We_l0_skrebkov_turbulentnyye_1963      = np.zeros(len(f_skrebkov_turbulentnyye_1963))
Fr_0_skrebkov_turbulentnyye_1963       = np.zeros(len(f_skrebkov_turbulentnyye_1963))
rho_s_skrebkov_turbulentnyye_1963      = np.zeros(len(f_skrebkov_turbulentnyye_1963))
I_0_skrebkov_turbulentnyye_1963        = np.zeros(len(f_skrebkov_turbulentnyye_1963))
nu_s_skrebkov_turbulentnyye_1963       = np.zeros(len(f_skrebkov_turbulentnyye_1963))
K_c_skrebkov_turbulentnyye_1963        = np.zeros(len(f_skrebkov_turbulentnyye_1963))
Ma_g_skrebkov_turbulentnyye_1963       = np.zeros(len(f_skrebkov_turbulentnyye_1963))
Lambda_r_s_skrebkov_turbulentnyye_1963 = np.zeros(len(f_skrebkov_turbulentnyye_1963))
regime_turb_skrebkov_turbulentnyye_1963 = []
for f in f_skrebkov_turbulentnyye_1963:
   Re_l0_skrebkov_turbulentnyye_1963[i]      = rho_l * Ubar_0 * d_0 / mu_l
   We_l0_skrebkov_turbulentnyye_1963[i]      = rho_l * Ubar_0**2. * d_0 / sigma
   Fr_0_skrebkov_turbulentnyye_1963[i]       = Ubar_0**2. / (g * d_0)
   I_0_skrebkov_turbulentnyye_1963[i]        = I_fully_developed(f_skrebkov_turbulentnyye_1963[i])
   Lambda_r_s_skrebkov_turbulentnyye_1963[i] = Lambda_s('rough', f_skrebkov_turbulentnyye_1963[i], 'v')
   rho_s_skrebkov_turbulentnyye_1963[i] = rho_l / rho_g
   nu_s_skrebkov_turbulentnyye_1963[i]  = nu_l / nu_g
   #K_c_skrebkov_turbulentnyye_1963[i] = (P_atm - P_v_water(T)) / (0.5 * rho_l * Ubar_0**2.)
   K_c_skrebkov_turbulentnyye_1963[i] = K_c(P_atm, P_v_water(T), rho_l, Ubar_0, K_L_guess, f_skrebkov_turbulentnyye_1963[i], 150.)
   Ma_g_skrebkov_turbulentnyye_1963[i]  = Ubar_0 / c_ideal_gas(gamma_air(T), T + T_zero, MW_air)
   
   f_laminar = friction_factor_laminar(Re_l0_skrebkov_turbulentnyye_1963[i])
   if (f / f_laminar > f_ratio_cutoff):
      regime_turb_skrebkov_turbulentnyye_1963.append('turbulent')
   else:
      regime_turb_skrebkov_turbulentnyye_1963.append('laminar')
   
   i = i + 1

df_skrebkov_turbulentnyye_1963 = pd.DataFrame({'We_l0': We_l0_skrebkov_turbulentnyye_1963})
df_skrebkov_turbulentnyye_1963['key']     = 'skrebkov_turbulentnyye_1963'
df_skrebkov_turbulentnyye_1963['alt key'] = np.nan

df_skrebkov_turbulentnyye_1963['liquid']     = 'water'
df_skrebkov_turbulentnyye_1963['gas']        = 'air'
df_skrebkov_turbulentnyye_1963['degas']      = np.nan
df_skrebkov_turbulentnyye_1963['regulation'] = np.nan

df_skrebkov_turbulentnyye_1963['bend']              = np.nan
df_skrebkov_turbulentnyye_1963['flow straightener'] = np.nan
df_skrebkov_turbulentnyye_1963['screen']            = np.nan
df_skrebkov_turbulentnyye_1963['contraction shape'] = np.nan
df_skrebkov_turbulentnyye_1963['c']                 = np.nan
df_skrebkov_turbulentnyye_1963['trip']              = np.nan
df_skrebkov_turbulentnyye_1963['L_r/d_0']           = 0.
df_skrebkov_turbulentnyye_1963['eps_r/d_0']         = np.nan

df_skrebkov_turbulentnyye_1963['L_0/d_0']       = 150. # Estimate; p. 142 (English translation) says the nozzles had a length of 100 to 200 diameters
df_skrebkov_turbulentnyye_1963['roughness']     = 'rough'
df_skrebkov_turbulentnyye_1963['f page fig']    = 'p. 143, fig. 1'
df_skrebkov_turbulentnyye_1963['pipe material'] = np.nan
df_skrebkov_turbulentnyye_1963['est f']         = False
df_skrebkov_turbulentnyye_1963['check FD']      = False

df_skrebkov_turbulentnyye_1963['t/d_0']       = np.nan
df_skrebkov_turbulentnyye_1963['L_tip/d_0']   = np.nan
df_skrebkov_turbulentnyye_1963['end checked'] = np.nan

df_skrebkov_turbulentnyye_1963['orientation']         = np.nan
df_skrebkov_turbulentnyye_1963['vibration isolation'] = np.nan
df_skrebkov_turbulentnyye_1963['d_chamber/d_0']       = np.nan

#df_skrebkov_turbulentnyye_1963['We_l0']      = We_l0_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['Re_l0']      = Re_l0_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['Fr_0']       = Fr_0_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['I_0']        = I_0_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['rho_s']      = rho_s_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['nu_s']       = nu_s_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['K_c']        = K_c_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['Ma_g']       = Ma_g_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['Lambda_r_s'] = Lambda_r_s_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['est rho_s']  = True
df_skrebkov_turbulentnyye_1963['est nu_s']   = True

df_skrebkov_turbulentnyye_1963['regime photo'] = np.nan
df_skrebkov_turbulentnyye_1963['regime L_b']   = np.nan
df_skrebkov_turbulentnyye_1963['regime turb']  = regime_turb_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['est turb']     = False

df_skrebkov_turbulentnyye_1963['L_b/d_0']                      = np.nan
df_skrebkov_turbulentnyye_1963['L_b/d_0 stat error']           = np.nan
df_skrebkov_turbulentnyye_1963['L_b/d_0 resolution']           = np.nan
df_skrebkov_turbulentnyye_1963['L_b method']                   = np.nan
df_skrebkov_turbulentnyye_1963['L_b/d_0 page fig']             = np.nan
df_skrebkov_turbulentnyye_1963['L_b/d_0 transcription method'] = np.nan

df_skrebkov_turbulentnyye_1963['theta']                      = theta_skrebkov_turbulentnyye_1963
df_skrebkov_turbulentnyye_1963['theta stat error']           = np.nan
df_skrebkov_turbulentnyye_1963['theta resolution']           = np.nan
df_skrebkov_turbulentnyye_1963['theta page fig']             = 'p. 143, fig. 1'
df_skrebkov_turbulentnyye_1963['theta transcription method'] = 'figure'

df_skrebkov_turbulentnyye_1963['D_10/d_0']                   = np.nan
df_skrebkov_turbulentnyye_1963['D_10/d_0 stat error']        = np.nan
df_skrebkov_turbulentnyye_1963['D_10/d_0 resolution']        = np.nan
df_skrebkov_turbulentnyye_1963['D_30/d_0']                   = np.nan
df_skrebkov_turbulentnyye_1963['D_30/d_0 stat error']        = np.nan
df_skrebkov_turbulentnyye_1963['D_30/d_0 resolution']        = np.nan
df_skrebkov_turbulentnyye_1963['D_32/d_0']                   = np.nan
df_skrebkov_turbulentnyye_1963['D_32/d_0 stat error']        = np.nan
df_skrebkov_turbulentnyye_1963['D_32/d_0 resolution']        = np.nan
df_skrebkov_turbulentnyye_1963['droplet x/d_0']              = np.nan
df_skrebkov_turbulentnyye_1963['D/d_0 page fig']             = np.nan
df_skrebkov_turbulentnyye_1963['D/d_0 transcription method'] = np.nan
df_skrebkov_turbulentnyye_1963['D/d_0 measurement method']   = np.nan
df_skrebkov_turbulentnyye_1963['v_d_bar/vp']                        = np.nan
df_skrebkov_turbulentnyye_1963['v_d_bar/vp stat error']             = np.nan
df_skrebkov_turbulentnyye_1963['v_d_bar/vp resolution']             = np.nan
df_skrebkov_turbulentnyye_1963['v_d_bar/vp page fig']               = np.nan
df_skrebkov_turbulentnyye_1963['v_d_bar/vp transcription method']   = np.nan
df_skrebkov_turbulentnyye_1963['e_p']                        = np.nan
df_skrebkov_turbulentnyye_1963['e_p page fig']               = np.nan
df_skrebkov_turbulentnyye_1963['e_p transcription method']   = np.nan

df_skrebkov_turbulentnyye_1963['x_trans/d_0']                      = np.nan
df_skrebkov_turbulentnyye_1963['x_trans/d_0 page fig']             = np.nan
df_skrebkov_turbulentnyye_1963['x_trans/d_0 transcription method'] = np.nan

df_skrebkov_turbulentnyye_1963['x_i/d_0']                      = np.nan
df_skrebkov_turbulentnyye_1963['x_i/d_0 stat error']           = np.nan
df_skrebkov_turbulentnyye_1963['x_i/d_0 resolution']           = np.nan
df_skrebkov_turbulentnyye_1963['x_i/d_0 page fig']             = np.nan
df_skrebkov_turbulentnyye_1963['x_i/d_0 transcription method'] = np.nan

df_skrebkov_turbulentnyye_1963['x_e/d_0']                      = np.nan
df_skrebkov_turbulentnyye_1963['x_e/d_0 stat error']           = np.nan
df_skrebkov_turbulentnyye_1963['x_e/d_0 resolution']           = np.nan
df_skrebkov_turbulentnyye_1963['x_e/d_0 page fig']             = np.nan
df_skrebkov_turbulentnyye_1963['x_e/d_0 transcription method'] = np.nan

df_skrebkov_turbulentnyye_1963['photo filename']   = np.nan
df_skrebkov_turbulentnyye_1963['exposure time']    = np.nan
df_skrebkov_turbulentnyye_1963['flash']            = np.nan
df_skrebkov_turbulentnyye_1963['x_low']            = np.nan
df_skrebkov_turbulentnyye_1963['x_mid']            = np.nan
df_skrebkov_turbulentnyye_1963['x_high']           = np.nan
df_skrebkov_turbulentnyye_1963['photo page fig']   = np.nan
df_skrebkov_turbulentnyye_1963['spectrum']         = np.nan
df_skrebkov_turbulentnyye_1963['lighting']         = np.nan
df_skrebkov_turbulentnyye_1963['background color'] = np.nan
df_skrebkov_turbulentnyye_1963['grid']             = np.nan
df_skrebkov_turbulentnyye_1963['distance']         = np.nan
df_skrebkov_turbulentnyye_1963['f-stop']           = np.nan
df_skrebkov_turbulentnyye_1963['focal length']     = np.nan
df_skrebkov_turbulentnyye_1963['camera model']     = np.nan
df_skrebkov_turbulentnyye_1963['sensor']           = np.nan

df_skrebkov_turbulentnyye_1963 = df_skrebkov_turbulentnyye_1963[cols]
summary_table(df_skrebkov_turbulentnyye_1963)
df_jet_breakup = pd.concat([df_jet_breakup, df_skrebkov_turbulentnyye_1963])

########################
# grant_newtonian_1965 #
########################

T = None
d_0 = None
i = None
j = None

# TODO: Add photos to this. (Are these done?)
# DONE: Add laminar data to help improve the regime map.
# DONE: Missing: table A-14, p. 234 (pdf p. 252)

# grant_newtonian_1965 p. 64: > it was possible to maintain the experimental apparatus at 25$^\deg$C $\pm$ 0.1$^\deg$C. All experiments were conducted at this temperature.
# grant_newtonian_1965 pp. 64-65: > Operating pressures were measured with a series of Ashcroft Test Gages, accurate to within 0.25% of the scale range. Gages of the following ranges were available: 0-60 psi., 0-300 psi., 0-600 psi., 0-1000 psi., 0-1500 psi., and 0-2000 psi. [...] The maximum error in pressure measurement associated with the head of liquid in the reservoir was approximately 1.5 psi. [...] it was possible to obtain chamber pressures down to 0.20 in. Hg absolute.
# grant_newtonian_1965 p. 74: > Efflux times were measured with an electric timer which could be read to 0.1 second. A total reaction time of 0.3 second in starting and stopping the timer made it advisable that time measurements of durations less than 30 seconds be avoided.
# grant_newtonian_1965 pp. 77-78 / grant_newtonian_1966 p. 673L: > In most cases the jets were allowed to run for a minimum of 15 to 20 sec. before the first series of from four to eight photographs were taken. Thus, data points at each velocity represent the average of a minimum of four separate photographs

# p. 64: > By employing this rather elaborate thermostating arrangement, it was possible to maintain the experimental apparatus at 25$^\circ~\pm~0.1^\circ$C.
T_atm = 25 # Assuming the air is at or near the same temperature.

# The subatmospheric cases do not have a typo in the pressure. See p. 109 of the dissertation.

# Issues with this data:
# 1. The densities given for the subatmospheric cases seem to be incorrect. Those found inconsistent with the ideal gas law at the given pressures were recalculated to match the ideal gas law.

# nozzles 1, 4, 5, 8, 9
# photos: pp. 83-84 (nozzle 1), 94-95 (nozzle 8), 101 (nozzle 8), 114 (nozzle 4)

grant_newtonian_1965_csv = pd.read_csv('../data/grant_newtonian_1965/grant_newtonian_1965.csv', sep=',', header=0)

page_grant_newtonian_1965           = grant_newtonian_1965_csv['page']
table_grant_newtonian_1965          = grant_newtonian_1965_csv['table']
fluid_grant_newtonian_1965          = grant_newtonian_1965_csv['liquid']
rho_l_grant_newtonian_1965          = grant_newtonian_1965_csv['rho_l']
sigma_grant_newtonian_1965          = grant_newtonian_1965_csv['sigma']
mu_l_grant_newtonian_1965           = grant_newtonian_1965_csv['mu_l']
d_0_grant_newtonian_1965            = grant_newtonian_1965_csv['d_0']
Ubar_0_grant_newtonian_1965         = grant_newtonian_1965_csv['Ubar_0']
L_0s_grant_newtonian_1965           = grant_newtonian_1965_csv['L_0/d_0']
L_bs_grant_newtonian_1965           = grant_newtonian_1965_csv['L_b/d_0']
regime_L_b_grant_newtonian_1965     = grant_newtonian_1965_csv['regime L_b']
rho_g_grant_newtonian_1965          = grant_newtonian_1965_csv['rho_g (g/cm^3)'] * 1000
P_atm_grant_newtonian_1965          = grant_newtonian_1965_csv['P_atm (mm Hg)'] * 133.322
nozzle_grant_newtonian_1965         = grant_newtonian_1965_csv['nozzle']
photo_filename_grant_newtonian_1965 = grant_newtonian_1965_csv['photo filename']
photo_page_fig_grant_newtonian_1965 = grant_newtonian_1965_csv['photo page fig']
regime_photo_grant_newtonian_1965   = grant_newtonian_1965_csv['regime photo']
x_low_grant_newtonian_1965          = grant_newtonian_1965_csv['x low']
x_tr_s_grant_newtonian_1965         = grant_newtonian_1965_csv['x_tr/d_0']

i = 0
Re_l0_grant_newtonian_1965       = np.zeros(len(Ubar_0_grant_newtonian_1965))
We_l0_grant_newtonian_1965       = np.zeros(len(Ubar_0_grant_newtonian_1965))
Fr_0_grant_newtonian_1965        = np.zeros(len(Ubar_0_grant_newtonian_1965))
rho_s_grant_newtonian_1965       = np.zeros(len(Ubar_0_grant_newtonian_1965))
nu_s_grant_newtonian_1965        = np.zeros(len(Ubar_0_grant_newtonian_1965))
nu_s_grant_newtonian_1965        = np.zeros(len(Ubar_0_grant_newtonian_1965))
K_c_grant_newtonian_1965         = np.zeros(len(Ubar_0_grant_newtonian_1965))
Ma_g_grant_newtonian_1965        = np.zeros(len(Ubar_0_grant_newtonian_1965))
L_bs_res_grant_newtonian_1965    = np.zeros(len(Ubar_0_grant_newtonian_1965))
e_L_bs_grant_newtonian_1965      = np.array([])
Re_turb_grant_newtonian_1965_array = np.zeros(len(Ubar_0_grant_newtonian_1965))
regime_turb_grant_newtonian_1965 = []
page_fig_grant_newtonian_1965    = []
est_turb_grant_newtonian_1965    = []
for Ubar_0, regime_L_b in zip(Ubar_0_grant_newtonian_1965, regime_L_b_grant_newtonian_1965):
   rho_l  = rho_l_grant_newtonian_1965[i]
   sigma  = sigma_grant_newtonian_1965[i]
   mu_l   = mu_l_grant_newtonian_1965[i]
   nu_l   = mu_l / rho_l
   d_0    = d_0_grant_newtonian_1965[i]
   fluid  = fluid_grant_newtonian_1965[i]
   nozzle = nozzle_grant_newtonian_1965[i]
   
   # TODO: Update rupe_dynamics_1962 to use the same glycerin-water data. No duplication.
   if fluid == 'distilled water':
      P_v = P_v_water(T_atm)
   elif fluid == 'glycerin-water (88%-12%)':
      P_v = np.interp(88, np.array([11.69, 11.97, 21.63, 45.02, 61.15, 77.38, 88.03]), np.array([230.9, 229.3, 222.3, 199.9, 173.4, 128.1, 79.8])) * 133.322 # Pa, from glycerine_producers_association_physical_1975 p. 8, tab. 12
   elif fluid == 'glycerin-water (72%-28%)':
      P_v = np.interp(72, np.array([11.69, 11.97, 21.63, 45.02, 61.15, 77.38, 88.03]), np.array([230.9, 229.3, 222.3, 199.9, 173.4, 128.1, 79.8])) * 133.322 # Pa, from glycerine_producers_association_physical_1975 p. 8, tab. 12
   elif fluid == 'ethylene-glycol':
      # http://www.dow.com/ethyleneglycol/about/properties.htm
      # 68 F
      P_v = 0.06 * 133.322 # Pa
   elif fluid == 'ethanol-water (95%-5%)':
      # TODO: Get better data or a better model for ethanol-water vapor pressure.
      P_v = np.interp(95, np.array([0, 100]), np.array([P_v_water(T_atm), 5.95e3]))
   else:
      raise ValueError('Invalid fluid for grant_newtonian_1965:', fluid)
   
   if np.isnan(rho_g_grant_newtonian_1965[i]):
      rho_g       = rho_ideal_gas(P_atm, T_atm, MW_air)
      P_atm_grant = P_atm
   else:
      rho_g       = rho_g_grant_newtonian_1965[i]
      P_atm_grant = P_atm_grant_newtonian_1965[i]
      #print rho_g, rho_ideal_gas(P_atm_grant, T_atm, MW_air), P_atm_grant / 101325.
      if not(abs(rho_g - rho_ideal_gas(P_atm_grant, T_atm, MW_air)) / rho_g < 1.e-2):
         rho_g = rho_ideal_gas(P_atm_grant, T_atm, MW_air)
   
   nu_g  = mu_ideal_gas(T_atm, mu_0_air, T_0_air, C_air) / rho_g
   
   Re_l0_grant_newtonian_1965[i] = rho_l * Ubar_0 * d_0 / mu_l
   We_l0_grant_newtonian_1965[i] = rho_l * Ubar_0**2 * d_0 / sigma
   Fr_0_grant_newtonian_1965[i] = Ubar_0**2 / (g * d_0)
   if not(np.isnan(page_grant_newtonian_1965[i])):
      page_fig_grant_newtonian_1965.append('p. '+str(page_grant_newtonian_1965[i])+', tab. '+table_grant_newtonian_1965[i])
   else:
      page_fig_grant_newtonian_1965.append(np.nan)
   rho_s_grant_newtonian_1965[i] = rho_l / rho_g
   nu_s_grant_newtonian_1965[i]  = nu_l / nu_g
   #K_c_grant_newtonian_1965[i]   = (P_atm_grant - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_grant_newtonian_1965[i]   = K_c(P_atm_grant, P_v, rho_l, Ubar_0, K_L_reentrant, friction_factor_smooth(Re_l0_grant_newtonian_1965[i]), L_0s_grant_newtonian_1965[i])
   Ma_g_grant_newtonian_1965[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_atm), T_atm + T_zero, MW_air)
   
   # grant_newtonian_1965 p. 73  grant_newtonian_1966 p. 672L: > All breakup data represent measurements made directly from negatives. Depending upon the order of magnitude of the dimensions involved, a traveling microscope accurate to 0.001 cm., a calibrated eyepiece accurate to 0.01 cm., or a simple scale accurate to 0.1 cm. was employed.
   # Assuming worst case scenario of 0.1 cm for the resolution. I also assume that the resolution given is that in physical space, not distance on the negative itself.
   L_bs_res_grant_newtonian_1965[i] = 0.1e-2 / d_0
   
   # DONE: Add transition Re from p. 143 (pdf p. 160).
   # nozzles 1, 4, 5, 8, 9
   if nozzle == 1:
      Re_trans_grant_newtonian_1965 = 2320
      Re_turb_grant_newtonian_1965  = 2540
      est_turb_grant_newtonian_1965.append(False)
   elif nozzle == 4:
      Re_trans_grant_newtonian_1965 = 2240
      Re_turb_grant_newtonian_1965  = 2410
      est_turb_grant_newtonian_1965.append(False)
   elif (nozzle == 5) or (nozzle == 6): # No information given, use averages given in the table.
      Re_trans_grant_newtonian_1965 = 2360
      Re_turb_grant_newtonian_1965  = 2490
      est_turb_grant_newtonian_1965.append(True)
   elif nozzle == 8:
      Re_trans_grant_newtonian_1965 = 0.5 * (2636 + 2350)
      Re_turb_grant_newtonian_1965  = 0.5 * (2670 + 2440)
      est_turb_grant_newtonian_1965.append(False)
   elif nozzle == 9:
      Re_trans_grant_newtonian_1965 = 0.5 * (2320 + 2290)
      Re_turb_grant_newtonian_1965  = 0.5 * (2500 + 2350)
      est_turb_grant_newtonian_1965.append(False)
   else:
      raise ValueError('Invalid nozzle for grant_newtonian_1965:', nozzle)
   
   assert(Re_turb_grant_newtonian_1965 > 0)
   Re_turb_grant_newtonian_1965_array[i] = Re_turb_grant_newtonian_1965
   
   if Re_l0_grant_newtonian_1965[i] > Re_turb_grant_newtonian_1965:
      regime_turb_grant_newtonian_1965.append('turbulent')
      
      if regime_L_b != 'Rayleigh':
         e_L_bs_grant_newtonian_1965 = np.append(e_L_bs_grant_newtonian_1965, scipy.stats.t.ppf(1. - (1. - interval_probability_level) / 2., 4 - 1) * breakup_length_sigmas * L_bs_grant_newtonian_1965[i] / sqrt(4.))
      else:
         e_L_bs_grant_newtonian_1965 = np.append(e_L_bs_grant_newtonian_1965, np.nan)
   elif Re_l0_grant_newtonian_1965[i] > Re_trans_grant_newtonian_1965:
      regime_turb_grant_newtonian_1965.append('transitional')
      e_L_bs_grant_newtonian_1965 = np.append(e_L_bs_grant_newtonian_1965, np.nan)
   else:
      regime_turb_grant_newtonian_1965.append('laminar')
      e_L_bs_grant_newtonian_1965 = np.append(e_L_bs_grant_newtonian_1965, np.nan)
   
   i = i + 1

I_0_grant_newtonian_1965        = I_fully_developed_array(Re_l0_grant_newtonian_1965, Re_turb_grant_newtonian_1965_array)
Lambda_r_s_grant_newtonian_1965 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_grant_newtonian_1965, Re_turb_grant_newtonian_1965_array), 'v')

df_grant_newtonian_1965 = pd.DataFrame({'We_l0': We_l0_grant_newtonian_1965})
df_grant_newtonian_1965['key']        = 'grant_newtonian_1965'
df_grant_newtonian_1965['alt key']    = 'grant_newtonian_1966'
df_grant_newtonian_1965['liquid']     = fluid_grant_newtonian_1965
df_grant_newtonian_1965['gas']        = 'air'
df_grant_newtonian_1965['degas']      = np.nan # TODO
df_grant_newtonian_1965['regulation'] = np.nan # TODO

df_grant_newtonian_1965['bend']              = False # See fig. 4-2 on p. 71.
df_grant_newtonian_1965['flow straightener'] = False
df_grant_newtonian_1965['contraction shape'] = 're-entrant'
df_grant_newtonian_1965['screen']            = False
df_grant_newtonian_1965['c']                 = 0.375 * 2.54e-2 / d_0_grant_newtonian_1965 # p. 63
df_grant_newtonian_1965['trip']              = False
df_grant_newtonian_1965['L_r/d_0']           = 0
df_grant_newtonian_1965['eps_r/d_0']         = np.nan

df_grant_newtonian_1965['L_0/d_0']       = L_0s_grant_newtonian_1965
df_grant_newtonian_1965['roughness']     = 'smooth'
df_grant_newtonian_1965['f page fig']    = np.nan
df_grant_newtonian_1965['pipe material'] = 'stainless steel' # grant_newtonian_1966 p. 671R
df_grant_newtonian_1965['est f']         = True
df_grant_newtonian_1965['check FD']      = False # Best I can tell.

df_grant_newtonian_1965['t/d_0']       = np.nan # TODO Find from photos?
df_grant_newtonian_1965['L_tip/d_0']   = np.nan
df_grant_newtonian_1965['end checked'] = True # p. 69

# fig. 4-1 on p. 67 must be an overhead shot
df_grant_newtonian_1965['orientation']         = 'horizontal' # grant_newtonian_1965 p. 62 (pdf p. 79)
df_grant_newtonian_1965['vibration isolation'] = True # p. 68
df_grant_newtonian_1965['d_chamber/d_0']       = 3.5 * 2.54e-2 / d_0_grant_newtonian_1965 # p. 60

#df_grant_newtonian_1965['We_l0']      = We_l0_grant_newtonian_1965
df_grant_newtonian_1965['Re_l0']      = Re_l0_grant_newtonian_1965
df_grant_newtonian_1965['I_0']        = I_0_grant_newtonian_1965
df_grant_newtonian_1965['Fr_0']       = Fr_0_grant_newtonian_1965
df_grant_newtonian_1965['rho_s']      = rho_s_grant_newtonian_1965
df_grant_newtonian_1965['nu_s']       = nu_s_grant_newtonian_1965
df_grant_newtonian_1965['K_c']        = K_c_grant_newtonian_1965
df_grant_newtonian_1965['Ma_g']       = Ma_g_grant_newtonian_1965
df_grant_newtonian_1965['Lambda_r_s'] = Lambda_r_s_grant_newtonian_1965
df_grant_newtonian_1965['est rho_s']  = False # Assuming the air temperature was 25 C as was everything else.
df_grant_newtonian_1965['est nu_s']  = False # See one line above.

df_grant_newtonian_1965['regime photo'] = regime_photo_grant_newtonian_1965
df_grant_newtonian_1965['regime L_b']   = regime_L_b_grant_newtonian_1965
df_grant_newtonian_1965['regime turb']  = regime_turb_grant_newtonian_1965
df_grant_newtonian_1965['est turb']     = est_turb_grant_newtonian_1965

df_grant_newtonian_1965['L_b/d_0'] = L_bs_grant_newtonian_1965
df_grant_newtonian_1965['L_b/d_0 stat error'] = e_L_bs_grant_newtonian_1965 # grant_newtonian_1965 p. 78: "series of from 4 to 8 photographs was taken." No standard deviation was provided. The standard deviation of the breakup length is estimated from phinney_breakup_1973.
df_grant_newtonian_1965['L_b/d_0 resolution'] = L_bs_res_grant_newtonian_1965
df_grant_newtonian_1965['L_b method'] = 'photographic'
df_grant_newtonian_1965['L_b/d_0 page fig'] = page_fig_grant_newtonian_1965
df_grant_newtonian_1965['L_b/d_0 transcription method'] = 'table'

df_grant_newtonian_1965['theta']                      = np.nan 
df_grant_newtonian_1965['theta stat error']           = np.nan 
df_grant_newtonian_1965['theta resolution']           = np.nan 
df_grant_newtonian_1965['theta page fig']             = np.nan 
df_grant_newtonian_1965['theta transcription method'] = np.nan 

df_grant_newtonian_1965['D_10/d_0']                   = np.nan 
df_grant_newtonian_1965['D_10/d_0 stat error']        = np.nan 
df_grant_newtonian_1965['D_10/d_0 resolution']        = np.nan 
df_grant_newtonian_1965['D_30/d_0']                   = np.nan 
df_grant_newtonian_1965['D_30/d_0 stat error']        = np.nan 
df_grant_newtonian_1965['D_30/d_0 resolution']        = np.nan 
df_grant_newtonian_1965['D_32/d_0']                   = np.nan 
df_grant_newtonian_1965['D_32/d_0 stat error']        = np.nan 
df_grant_newtonian_1965['D_32/d_0 resolution']        = np.nan 
df_grant_newtonian_1965['droplet x/d_0']              = np.nan 
df_grant_newtonian_1965['D/d_0 page fig']             = np.nan 
df_grant_newtonian_1965['D/d_0 transcription method'] = np.nan 
df_grant_newtonian_1965['D/d_0 measurement method']   = np.nan
df_grant_newtonian_1965['v_d_bar/vp']                        = np.nan
df_grant_newtonian_1965['v_d_bar/vp stat error']             = np.nan
df_grant_newtonian_1965['v_d_bar/vp resolution']             = np.nan
df_grant_newtonian_1965['v_d_bar/vp page fig']               = np.nan
df_grant_newtonian_1965['v_d_bar/vp transcription method']   = np.nan
df_grant_newtonian_1965['e_p']                        = np.nan
df_grant_newtonian_1965['e_p page fig']               = np.nan
df_grant_newtonian_1965['e_p transcription method']   = np.nan

df_grant_newtonian_1965['x_trans/d_0']                      = x_tr_s_grant_newtonian_1965
df_grant_newtonian_1965['x_trans/d_0 page fig']             = np.nan
df_grant_newtonian_1965['x_trans/d_0 transcription method'] = np.nan

df_grant_newtonian_1965['x_i/d_0']                      = np.nan
df_grant_newtonian_1965['x_i/d_0 stat error']           = np.nan
df_grant_newtonian_1965['x_i/d_0 resolution']           = np.nan
df_grant_newtonian_1965['x_i/d_0 page fig']             = np.nan
df_grant_newtonian_1965['x_i/d_0 transcription method'] = np.nan

df_grant_newtonian_1965['x_e/d_0']                      = np.nan
df_grant_newtonian_1965['x_e/d_0 stat error']           = np.nan
df_grant_newtonian_1965['x_e/d_0 resolution']           = np.nan
df_grant_newtonian_1965['x_e/d_0 page fig']             = np.nan
df_grant_newtonian_1965['x_e/d_0 transcription method'] = np.nan

df_grant_newtonian_1965['photo filename']   = photo_filename_grant_newtonian_1965
df_grant_newtonian_1965['exposure time']    = 0.5e-6
df_grant_newtonian_1965['flash']            = True # grant_newtonian_1965 p. 72: > flash of 0.5 microsecond duration
df_grant_newtonian_1965['x_low']            = x_low_grant_newtonian_1965
df_grant_newtonian_1965['x_mid']            = np.nan # TODO
df_grant_newtonian_1965['x_high']           = np.nan # TODO
df_grant_newtonian_1965['photo page fig']   = photo_page_fig_grant_newtonian_1965
df_grant_newtonian_1965['spectrum']         = 'visible'
df_grant_newtonian_1965['lighting']         = 'background' # p. 71
df_grant_newtonian_1965['background color'] = np.nan
df_grant_newtonian_1965['grid']             = np.nan
df_grant_newtonian_1965['distance']         = '2 to 8 feet' # p. 73
df_grant_newtonian_1965['f-stop']           = 4.7 # minimum, f/11 or f/22 were used; see pp. 72-73
df_grant_newtonian_1965['focal length']     = 135 # mm
df_grant_newtonian_1965['camera model']     = '4x5 Graflex Crown Graphic' # p. 72
df_grant_newtonian_1965['sensor']           = 'ASA 400 4x5 Kodak Royal Pan film' # p. 72

df_grant_newtonian_1965 = df_grant_newtonian_1965[cols]
summary_table(df_grant_newtonian_1965)
df_jet_breakup = pd.concat([df_jet_breakup, df_grant_newtonian_1965])

#####################
# kusui_liquid_1969 #
#####################

T = None
d_0 = None
i = None
j = None

# Issues with this data:
# 1. Tu unknown aside from smooth pipes due to rough-to-smooth transition. Tu seems to be higher than Smits' data suggests it should be. The rough-to-smooth transition location may increase Tu if the transition is abrupt or poorly set up.

# DONE: It does not appear that they provide any information relevant for a quantitative error analysis. However, they did use the electrical conductivity method to measure breakup length, which is better than photographs (equivalent to very large N).

# TODO: Use interpolated friction factor data from paper rather than asymptotic value.
# TODO: Add data from figure 4.
# DONE: Adjust Tu estimate from kusui_liquid_1969 to take into account the 8.75 d_0 long smooth section at the end of the pipe. Combine standard turbulent modeling approaches and turbulent Bernoulli to get friction factor term? If how to solve the differential equations is not obvious, use numerics. A single space step may be sufficient.
# TODO: Add atomization breakup lengths, and categorize them as in the atomization regime.
# siuru_response_1975 p. 211: > Previous investigations have shown that the transition from rough to smooth wall flow is slower than for the smooth to rough transition.
# rough to smooth: antonia_response_1972 (BL), tani_response_1971, Smits

kusui_liquid_1969_csv = pd.read_csv('../data/kusui_liquid_1969/kusui_liquid_1969.csv', sep=',', header=0)

page_kusui_liquid_1969       = kusui_liquid_1969_csv['page']
fig_kusui_liquid_1969        = kusui_liquid_1969_csv['figure']
roughness_kusui_liquid_1969  = kusui_liquid_1969_csv['roughness']
T_l_kusui_liquid_1969        = kusui_liquid_1969_csv['T_l']
d_0_kusui_liquid_1969        = kusui_liquid_1969_csv['d_0']
f_kusui_liquid_1969          = kusui_liquid_1969_csv['f']
Re_l0_kusui_liquid_1969      = kusui_liquid_1969_csv['Re_l0']
L_tips_kusui_liquid_1969     = kusui_liquid_1969_csv['L_tip/d_0']
L_bs_kusui_liquid_1969       = kusui_liquid_1969_csv['L_b/d_0']
regime_L_b_kusui_liquid_1969 = kusui_liquid_1969_csv['regime L_b']

i = 0
We_l0_kusui_liquid_1969       = np.zeros(len(Re_l0_kusui_liquid_1969))
Fr_0_kusui_liquid_1969        = np.zeros(len(Re_l0_kusui_liquid_1969))
I_0_kusui_liquid_1969         = np.zeros(len(Re_l0_kusui_liquid_1969))
rho_s_kusui_liquid_1969       = np.zeros(len(Re_l0_kusui_liquid_1969))
nu_s_kusui_liquid_1969        = np.zeros(len(Re_l0_kusui_liquid_1969))
K_c_kusui_liquid_1969         = np.zeros(len(Re_l0_kusui_liquid_1969))
Ma_g_kusui_liquid_1969        = np.zeros(len(Re_l0_kusui_liquid_1969))
Lambda_r_s_kusui_liquid_1969  = np.zeros(len(Re_l0_kusui_liquid_1969))
page_fig_kusui_liquid_1969    = []
regime_turb_kusui_liquid_1969 = []
e_L_bs_kusui_liquid_1969      = []
for d_0 in d_0_kusui_liquid_1969:
   T_l       = T_l_kusui_liquid_1969[i]
   rho_l     = rho_water(T_l)
   sigma     = sigma_water(T_l)
   mu_l      = mu_water(T_l)
   nu_l      = nu_water(T_l)
   P_v       = P_v_water(T_l)
   f         = f_kusui_liquid_1969[i]
   Re_l0     = Re_l0_kusui_liquid_1969[i]
   d_0       = d_0_kusui_liquid_1969[i]
   roughness = roughness_kusui_liquid_1969[i]
   
   Ubar_0 = Re_l0 * nu_l / d_0
   
   Oh_l0 = mu_l / sqrt(rho_l * sigma * d_0)
   #We_l0_kusui_liquid_1969[i] = (Oh_l0 * Re_l0)**2.
   We_l0_kusui_liquid_1969[i] = rho_l * Ubar_0**2. * d_0 / sigma
   page_fig_kusui_liquid_1969.append('p. '+str(page_kusui_liquid_1969[i])+', fig. '+fig_kusui_liquid_1969[i])
   I_0_kusui_liquid_1969[i] = I_fully_developed(f)
   Lambda_r_s_kusui_liquid_1969[i] = Lambda_s(roughness, f, 'v')
   if roughness == 'rough':
      # interp1([7.5, 11.3], [15.24, 11.04], 8.75) / 22.07
      # Based on Smits data from bad photograph of APS talk. Assumes that Tu drops negligibly through 3.3 diameters from the nozzle.
      # Assumes that the fraction Tu drops is only a function of x_smooth/d_0. An actual turbulence model might do better.
      Tu_FDsmooth = I_fully_developed(friction_factor_smooth(Re_l0))
      Tu_FDrough = I_0_kusui_liquid_1969[i]
      r2s_L = 0.23568 # Tani
      Re_L  = 76000.
      r2s_H = 0.40001 # Smits
      Re_H  = 131000.
      a = (r2s_H -r2s_L) / (Re_H - Re_L)
      b = r2s_L - a * Re_L
      Re_max = (1 - b) / a
      Re_min = -b / a
      
      if Re_l0 > Re_max:
         print 'Re_max violated:', Re_l0, Re_max
         sys.exit()
      elif Re_l0 < Re_min:
         print 'Re_min violated:', Re_l0, Re_min
         #sys.exit()
         r2s = 0
      
      #r2s = a * Re_l0 + b
      r2s = 1
      #print Re_l0, Re_max, r2s
      I_0_kusui_liquid_1969[i] = sqrt(Tu_FDsmooth**2. + (Tu_FDrough**2. - Tu_FDsmooth**2.) * r2s) # 0.400 for Smits; 0.334 for average of Tani and Smits; ~0.7 seems to make the exponents of We and 1/Tu match, which is what I might expect based on an earlier theory.
      
      #1.0   0.53337
      #0.9   0.52979
      #0.7   0.42029
      #0.4   -1.51934
   
   rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
   nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
   
   Fr_0_kusui_liquid_1969[i]  = Ubar_0**2 / (g * d_0)
   rho_s_kusui_liquid_1969[i] = rho_l / rho_g
   nu_s_kusui_liquid_1969[i]  = nu_l / nu_g
   #K_c_kusui_liquid_1969[i]   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_kusui_liquid_1969[i]   = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_wellrounded + friction_factor_smooth(Re_l0) * L_tips_kusui_liquid_1969[i], f, 30.)
   Ma_g_kusui_liquid_1969[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)
   
   f_laminar = friction_factor_laminar(Re_l0)
   if (f / f_laminar > f_ratio_cutoff):
      regime_turb_kusui_liquid_1969.append('turbulent')
   else:
      regime_turb_kusui_liquid_1969.append('laminar')
   
   e_L_bs_kusui_liquid_1969.append(1.e-2 / d_0) # Guessed 1 cm precision.
   
   i = i + 1

df_kusui_liquid_1969 = pd.DataFrame({'We_l0': We_l0_kusui_liquid_1969})
df_kusui_liquid_1969['key']        = 'kusui_liquid_1969'
df_kusui_liquid_1969['alt key']    = 'kusui_liquid_1969-1'
df_kusui_liquid_1969['liquid']     = 'water'
df_kusui_liquid_1969['gas']        = 'air'
df_kusui_liquid_1969['degas']      = False
df_kusui_liquid_1969['regulation'] = 'none' # Page 1063L suggests that the air compressor is turned off during the test as to avoid flow pulsations.

df_kusui_liquid_1969['bend']              = False
df_kusui_liquid_1969['flow straightener'] = False
df_kusui_liquid_1969['contraction shape'] = 'smooth'
df_kusui_liquid_1969['screen']            = False
df_kusui_liquid_1969['c']                 = np.inf # Fig. 1 indicates the water came from a large reservoir.
df_kusui_liquid_1969['trip']              = False
df_kusui_liquid_1969['L_r/d_0']           = 0.
df_kusui_liquid_1969['eps_r/d_0']         = np.nan

df_kusui_liquid_1969['L_0/d_0']       = 30.
df_kusui_liquid_1969['roughness']     = roughness_kusui_liquid_1969
df_kusui_liquid_1969['f page fig']    = 'p. 1065, fig. 5'
df_kusui_liquid_1969['pipe material'] = np.nan # The material does not appear to be mentioned.
df_kusui_liquid_1969['est f']         = False
df_kusui_liquid_1969['check FD']      = True

df_kusui_liquid_1969['t/d_0']       = np.nan # Appears large based on fig. 1.
df_kusui_liquid_1969['L_tip/d_0']   = L_tips_kusui_liquid_1969
df_kusui_liquid_1969['end checked'] = True # Based on the discussion it seems likely the author made sure to use a good tip.

df_kusui_liquid_1969['orientation']         = 'horizontal' # Presumably the same as kusui_liquid_1968
df_kusui_liquid_1969['vibration isolation'] = np.nan 
df_kusui_liquid_1969['d_chamber/d_0']       = np.inf

df_kusui_liquid_1969['L_b/d_0']                      = L_bs_kusui_liquid_1969
df_kusui_liquid_1969['L_b/d_0 page fig']             = page_fig_kusui_liquid_1969
df_kusui_liquid_1969['L_b/d_0 stat error']           = np.nan
df_kusui_liquid_1969['L_b/d_0 resolution']           = e_L_bs_kusui_liquid_1969
df_kusui_liquid_1969['L_b method']                   = 'electrical'
df_kusui_liquid_1969['L_b/d_0 transcription method'] = 'figure'

#df_kusui_liquid_1969['We_l0']      = We_l0_kusui_liquid_1969
df_kusui_liquid_1969['Re_l0']      = Re_l0_kusui_liquid_1969
df_kusui_liquid_1969['I_0']        = I_0_kusui_liquid_1969
df_kusui_liquid_1969['Fr_0']       = Fr_0_kusui_liquid_1969
df_kusui_liquid_1969['rho_s']      = rho_s_kusui_liquid_1969
df_kusui_liquid_1969['nu_s']       = nu_s_kusui_liquid_1969
df_kusui_liquid_1969['K_c']        = K_c_kusui_liquid_1969
df_kusui_liquid_1969['Ma_g']       = Ma_g_kusui_liquid_1969
df_kusui_liquid_1969['Lambda_r_s'] = Lambda_r_s_kusui_liquid_1969
df_kusui_liquid_1969['est rho_s']  = True
df_kusui_liquid_1969['est nu_s']   = True

df_kusui_liquid_1969['regime photo'] = np.nan
df_kusui_liquid_1969['regime L_b']   = regime_L_b_kusui_liquid_1969
df_kusui_liquid_1969['regime turb']  = regime_turb_kusui_liquid_1969
df_kusui_liquid_1969['est turb']     = False

df_kusui_liquid_1969['theta']                      = np.nan
df_kusui_liquid_1969['theta stat error']           = np.nan
df_kusui_liquid_1969['theta resolution']           = np.nan
df_kusui_liquid_1969['theta page fig']             = np.nan
df_kusui_liquid_1969['theta transcription method'] = np.nan

df_kusui_liquid_1969['D_10/d_0']                   = np.nan
df_kusui_liquid_1969['D_10/d_0 stat error']        = np.nan
df_kusui_liquid_1969['D_10/d_0 resolution']        = np.nan
df_kusui_liquid_1969['D_30/d_0']                   = np.nan
df_kusui_liquid_1969['D_30/d_0 stat error']        = np.nan
df_kusui_liquid_1969['D_30/d_0 resolution']        = np.nan
df_kusui_liquid_1969['D_32/d_0']                   = np.nan
df_kusui_liquid_1969['D_32/d_0 stat error']        = np.nan
df_kusui_liquid_1969['D_32/d_0 resolution']        = np.nan
df_kusui_liquid_1969['droplet x/d_0']              = np.nan
df_kusui_liquid_1969['D/d_0 page fig']             = np.nan
df_kusui_liquid_1969['D/d_0 transcription method'] = np.nan
df_kusui_liquid_1969['D/d_0 measurement method']   = np.nan
df_kusui_liquid_1969['v_d_bar/vp']                        = np.nan
df_kusui_liquid_1969['v_d_bar/vp stat error']             = np.nan
df_kusui_liquid_1969['v_d_bar/vp resolution']             = np.nan
df_kusui_liquid_1969['v_d_bar/vp page fig']               = np.nan
df_kusui_liquid_1969['v_d_bar/vp transcription method']   = np.nan
df_kusui_liquid_1969['e_p']                        = np.nan
df_kusui_liquid_1969['e_p page fig']               = np.nan
df_kusui_liquid_1969['e_p transcription method']   = np.nan

df_kusui_liquid_1969['x_trans/d_0']                      = np.nan
df_kusui_liquid_1969['x_trans/d_0 page fig']             = np.nan
df_kusui_liquid_1969['x_trans/d_0 transcription method'] = np.nan

df_kusui_liquid_1969['x_i/d_0']                      = np.nan
df_kusui_liquid_1969['x_i/d_0 stat error']           = np.nan
df_kusui_liquid_1969['x_i/d_0 resolution']           = np.nan
df_kusui_liquid_1969['x_i/d_0 page fig']             = np.nan
df_kusui_liquid_1969['x_i/d_0 transcription method'] = np.nan

df_kusui_liquid_1969['x_e/d_0']                      = np.nan
df_kusui_liquid_1969['x_e/d_0 stat error']           = np.nan
df_kusui_liquid_1969['x_e/d_0 resolution']           = np.nan
df_kusui_liquid_1969['x_e/d_0 page fig']             = np.nan
df_kusui_liquid_1969['x_e/d_0 transcription method'] = np.nan

df_kusui_liquid_1969['photo filename']   = np.nan
df_kusui_liquid_1969['exposure time']    = np.nan
df_kusui_liquid_1969['flash']            = np.nan
df_kusui_liquid_1969['x_low']            = np.nan
df_kusui_liquid_1969['x_mid']            = np.nan
df_kusui_liquid_1969['x_high']           = np.nan
df_kusui_liquid_1969['photo page fig']   = np.nan
df_kusui_liquid_1969['spectrum']         = np.nan
df_kusui_liquid_1969['lighting']         = np.nan
df_kusui_liquid_1969['background color'] = np.nan
df_kusui_liquid_1969['grid']             = np.nan
df_kusui_liquid_1969['distance']         = np.nan
df_kusui_liquid_1969['f-stop']           = np.nan
df_kusui_liquid_1969['focal length']     = np.nan
df_kusui_liquid_1969['camera model']     = np.nan
df_kusui_liquid_1969['sensor']           = np.nan

df_kusui_liquid_1969 = df_kusui_liquid_1969[cols]
summary_table(df_kusui_liquid_1969)
df_jet_breakup = pd.concat([df_jet_breakup, df_kusui_liquid_1969])

##########################
# phinney_stability_1970 #
##########################

# Issues with this data:
# 1. The Ohnesorge number implied by the fluid properties for solution 8, nozzle 2 is inconsistent with that printed in the tables. The data also is fairly far off that for other authors for this case. Solution 8 does not appear in phinney_stability_1972, so I can't use that to get insight into this. I believe that this should have said nozzle 1 as the Ohnesorge number matches in that case.

# phinney_stability_1972 has the same data. See p. 433R.

# p. 6: > it was felt that pressure measurements were accurate to within two percent.
# p. 9: > [...] The combined errors in elapsed time and weight amount to an error in exit velocity on the order of 2 percent.
# p. 9: > The accuracy of position measurement was approximately .2 cm.

phinney_stability_1970_csv = pd.read_csv('../data/phinney_stability_1970/phinney_stability_1970.csv', sep=',', header=0)

page_fig_phinney_stability_1970   = phinney_stability_1970_csv['page fig']
nozzle_phinney_stability_1970     = phinney_stability_1970_csv['nozzle']
solution_phinney_stability_1970   = phinney_stability_1970_csv['solution']
Ubar_0_phinney_stability_1970     = phinney_stability_1970_csv['Ubar_0 (cm/s)'] * 1e-2 # m/s
L_b_phinney_stability_1970        = phinney_stability_1970_csv['L_b (cm)'] * 1e-2 # m
regime_L_b_phinney_stability_1970 = phinney_stability_1970_csv['regime L_b']

i = 0
d_0_phinney_stability_1970    = np.zeros(len(Ubar_0_phinney_stability_1970))
L_0s_phinney_stability_1970   = np.zeros(len(Ubar_0_phinney_stability_1970))
e_L_bs_phinney_stability_1970 = np.zeros(len(Ubar_0_phinney_stability_1970))
for nozzle in nozzle_phinney_stability_1970:
   if nozzle == 1:
      d_0_phinney_stability_1970[i]    = 0.125e-2 # m
      L_0s_phinney_stability_1970[i]   = 17.845e-2 / d_0_phinney_stability_1970[i]
   elif nozzle == 2:
      d_0_phinney_stability_1970[i]    = 0.05041e-2 # m
      L_0s_phinney_stability_1970[i]   = 7.43e-2 / d_0_phinney_stability_1970[i]
   elif nozzle == 3:
      d_0_phinney_stability_1970[i]    = 0.05571e-2 # m
      L_0s_phinney_stability_1970[i]   = 123.0e-2 / d_0_phinney_stability_1970[i]
   else:
      raise ValueError('Invalid nozzle for phinney_stability_1970.')
   
   e_L_bs_phinney_stability_1970[i] = 0.2e-2 / d_0_phinney_stability_1970[i] # p. 9: > The accuracy of position measurement was approximately .2 cm.
   i = i + 1

i = 0

rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g

We_l0_phinney_stability_1970 = np.zeros(len(Ubar_0_phinney_stability_1970))
Re_l0_phinney_stability_1970 = np.zeros(len(Ubar_0_phinney_stability_1970))
Fr_0_phinney_stability_1970  = np.zeros(len(Ubar_0_phinney_stability_1970))
rho_s_phinney_stability_1970 = np.zeros(len(Ubar_0_phinney_stability_1970))
nu_s_phinney_stability_1970  = np.zeros(len(Ubar_0_phinney_stability_1970))
K_c_phinney_stability_1970   = np.zeros(len(Ubar_0_phinney_stability_1970))
Ma_g_phinney_stability_1970  = np.zeros(len(Ubar_0_phinney_stability_1970))
liquid_phinney_stability_1970 = []
for solution, nozzle in zip(solution_phinney_stability_1970, nozzle_phinney_stability_1970): # p. 27
   d_0    = d_0_phinney_stability_1970[i]
   Ubar_0 = Ubar_0_phinney_stability_1970[i]
   liquid_phinney_stability_1970.append('Karo syrup, water, and salt mixture (solution '+str(solution)+')')
   
   if solution == 1:
      rho_l = 1.294e3  # kg/m^3
      sigma = 76.1e-3  # N/m
      mu_l  = 0.470e-1 # Pa*s
   elif solution == 2:
      rho_l = 1.2575e3 # kg/m^3
      sigma = 72.1e-3  # N/m
      mu_l  = 0.238e-1 # Pa*s
   elif solution == 3:
      rho_l = 1.152e3  # kg/m^3
      sigma = 64.1e-3  # N/m
      mu_l  = 0.038e-1 # Pa*s
   elif solution == 4:
      rho_l = 1.131e3   # kg/m^3
      sigma = 64.4e-3   # N/m
      mu_l  = 0.0272e-1 # Pa*s
   elif solution == 5:
      rho_l = 1.236e3  # kg/m^3
      sigma = 75.7e-3  # N/m
      mu_l  = 0.14e-1  # Pa*s
   elif solution == 6:
      rho_l = 1.192e3  # kg/m^3
      sigma = 70.3e-3  # N/m
      mu_l  = 0.068e-1 # Pa*s
   elif solution == 7:
      rho_l = 1.386e3 # kg/m^3
      sigma = 83.0e-3 # N/m
      mu_l  = 4.10e-1 # Pa*s
   elif solution == 8:
      rho_l = 1.359e3 # kg/m^3
      sigma = 86.7e-3 # N/m
      mu_l  = 2.51e-1 # Pa*s
   elif solution == 9:
      rho_l = 1.330e3 # kg/m^3
      sigma = 84.5e-3 # N/m
      mu_l  = 0.85e-1 # Pa*s
   elif solution == 10:
      rho_l = 1.051e3    # kg/m^3
      sigma = 72.5e-3    # N/m
      mu_l  = 0.01023e-1 # Pa*s
   elif solution == 11:
      rho_l = 1.076e3   # kg/m^3
      sigma = 66.0e-3   # N/m
      mu_l  = 0.0164e-1 # Pa*s
   elif solution == 12:
      rho_l = 1.0242e3   # kg/m^3
      sigma = 69.5e-3    # N/m
      mu_l  = 0.01069e-1 # Pa*s
   else:
      raise ValueError('Invalid solution for phinney_stability_1970.')
   
   # Check Ohnesorge number consistent with given fluid properties against that printed in the tables.
   Oh_l0 = mu_l / np.sqrt(rho_l * sigma * d_0)
   
   if solution == 1:
      if nozzle == 1:
         assert(abs(Oh_l0 - 0.133963) < 1.e-4)
   elif solution == 2:
      if nozzle == 1:
         assert(abs(Oh_l0 - 7.06969e-2) < 1.e-4)
      elif nozzle == 2:
         assert(abs(Oh_l0 - 0.111326) < 1.e-4)
   elif solution == 3:
      if nozzle == 1:
         assert(abs(Oh_l0 - 1.25076e-2) < 1.e-4)
      elif nozzle == 2:
         assert(abs(Oh_l0 - 1.96956e-2) < 1.e-4)
   elif solution == 5:
      if nozzle == 1:
         assert(abs(Oh_l0 - 0.040937) < 1.e-4)
      elif nozzle == 2:
         assert(abs(Oh_l0 - 6.44633e-2) < 1.e-4)
   elif solution == 6:
      if nozzle == 1:
         assert(abs(Oh_l0 - 0.020993) < 1.e-4)
      elif nozzle == 1:
         assert(abs(Oh_l0 - 3.30576e-2) < 1.e-4)
   elif solution == 8:
      if nozzle == 1: # printed as nozzle 2 in the data
         assert(abs(Oh_l0 - 0.654033) < 1.e-4)
   elif solution == 9:
      if nozzle == 1:
         assert(abs(Oh_l0 - 0.226783) < 1.e-4)
      elif nozzle == 2:
         assert(abs(Oh_l0 - 0.357114) < 1.e-4)
   elif solution == 10:
      if nozzle == 3:
         assert(abs(Oh_l0 - 4.96522e-3) < 1.e-4)
   elif solution == 11:
      if nozzle == 3:
         assert(abs(Oh_l0 - 8.24516e-3) < 1.e-4)
   elif solution == 12:
      if nozzle == 3:
         assert(abs(Oh_l0 - 5.36817e-3) < 1.e-4)
   else:
      raise ValueError('Invalid solution for phinney_stability_1970.')
   
   assert(rho_l < 1500.)
   assert(rho_l > 1000.)
   assert(sigma < 1.e-1)
   assert(sigma > 1.e-3)
   assert(mu_l < 5.e-1)
   assert(mu_l > 0.01e-1)
   
   nu_l = mu_l / rho_l
   P_v  = P_v_water(T_std) # TODO: Later change this to be for the correct fluid. This is for water when it actually needs a mixture of water and Karo syrup.
   
   We_l0_phinney_stability_1970[i] = rho_l * Ubar_0**2 * d_0 / sigma
   Re_l0_phinney_stability_1970[i] = rho_l * Ubar_0 * d_0 / mu_l
   Fr_0_phinney_stability_1970[i]  = Ubar_0**2 / (g * d_0)
   rho_s_phinney_stability_1970[i] = rho_l / rho_g
   nu_s_phinney_stability_1970[i]  = nu_l / nu_g
   #K_c_phinney_stability_1970[i]   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_phinney_stability_1970[i]   = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_reentrant, friction_factor_smooth(Re_l0_phinney_stability_1970[i]), L_0s_phinney_stability_1970[i])
   Ma_g_phinney_stability_1970[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)
   
   # p. 8: > The density is measured [...] The accuracy is certainly within one percent.
   # p. 8: > > A series of readings on fresh samples were taken until it was felt that the value of $\sigma$ could be defined within $\pm.5$ percent.
   # p. 9: > For any given test (on one pipe) all values of $\mu$ would be within a band $\pm 1$ percent. For different pipes the values of \$mu$ would be within 3 percent.
   i = i + 1

I_0_phinney_stability_1970        = I_fully_developed_array(Re_l0_phinney_stability_1970, Re_turb)
Lambda_r_s_phinney_stability_1970 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_phinney_stability_1970, Re_turb), 'v')
L_bs_phinney_stability_1970       = L_b_phinney_stability_1970 / d_0_phinney_stability_1970

df_phinney_stability_1970 = pd.DataFrame({'We_l0': We_l0_phinney_stability_1970})
df_phinney_stability_1970['key']     = 'phinney_stability_1970'
df_phinney_stability_1970['alt key'] = 'phinney_stability_1972'

df_phinney_stability_1970['liquid']     = liquid_phinney_stability_1970 # p. 8
df_phinney_stability_1970['gas']        = 'air'
df_phinney_stability_1970['degas']      = np.nan
df_phinney_stability_1970['regulation'] = True # p. 6

df_phinney_stability_1970['bend']              = False
df_phinney_stability_1970['flow straightener'] = False
df_phinney_stability_1970['screen']            = False
df_phinney_stability_1970['contraction shape'] = 're-entrant' # fig. 3, pdf p. 53
df_phinney_stability_1970['c']                 = 0.625*2.54e-2 / d_0_phinney_stability_1970 # 0.625 inch diameter inlet pipe
df_phinney_stability_1970['trip']              = False
df_phinney_stability_1970['L_r/d_0']           = 0
df_phinney_stability_1970['eps_r/d_0']         = np.nan

df_phinney_stability_1970['L_0/d_0']       = L_0s_phinney_stability_1970
df_phinney_stability_1970['roughness']     = 'smooth'
df_phinney_stability_1970['f page fig']    = np.nan
df_phinney_stability_1970['pipe material'] = 'glass' # p. 6
df_phinney_stability_1970['est f']         = True
df_phinney_stability_1970['check FD']      = np.nan

df_phinney_stability_1970['t/d_0']       = (0.25 * 2.54e-2 - d_0_phinney_stability_1970) / (2 * d_0_phinney_stability_1970) # fig. 3 says the nominal outer diameter of the glass tubes is 1/4"
df_phinney_stability_1970['L_tip/d_0']   = np.nan
df_phinney_stability_1970['end checked'] = True # p. 6: > Special care was taken to make sure the neds of the nozzles were square.

df_phinney_stability_1970['orientation']         = 'horizontal'
df_phinney_stability_1970['vibration isolation'] = True # p. 6: > To minimize disturbances...
df_phinney_stability_1970['d_chamber/d_0']       = np.inf

#df_phinney_stability_1970['We_l0']      = We_l0_phinney_stability_1970
df_phinney_stability_1970['Re_l0']      = Re_l0_phinney_stability_1970
df_phinney_stability_1970['Fr_0']       = Fr_0_phinney_stability_1970
df_phinney_stability_1970['I_0']        = I_0_phinney_stability_1970
df_phinney_stability_1970['rho_s']      = rho_s_phinney_stability_1970
df_phinney_stability_1970['nu_s']       = nu_s_phinney_stability_1970
df_phinney_stability_1970['K_c']        = K_c_phinney_stability_1970
df_phinney_stability_1970['Ma_g']       = Ma_g_phinney_stability_1970
df_phinney_stability_1970['Lambda_r_s'] = Lambda_r_s_phinney_stability_1970
df_phinney_stability_1970['est rho_s']  = True
df_phinney_stability_1970['est nu_s']   = True

df_phinney_stability_1970['regime photo'] = np.nan
df_phinney_stability_1970['regime L_b']   = regime_L_b_phinney_stability_1970
df_phinney_stability_1970['regime turb']  = 'laminar'
df_phinney_stability_1970['est turb']     = False

df_phinney_stability_1970['L_b/d_0']                      = L_bs_phinney_stability_1970
df_phinney_stability_1970['L_b/d_0 stat error']           = np.nan
df_phinney_stability_1970['L_b/d_0 resolution']           = e_L_bs_phinney_stability_1970
df_phinney_stability_1970['L_b method']                   = 'electrical'
df_phinney_stability_1970['L_b/d_0 page fig']             = page_fig_phinney_stability_1970
df_phinney_stability_1970['L_b/d_0 transcription method'] = 'table'

df_phinney_stability_1970['theta']                      = np.nan
df_phinney_stability_1970['theta stat error']           = np.nan
df_phinney_stability_1970['theta resolution']           = np.nan
df_phinney_stability_1970['theta page fig']             = np.nan
df_phinney_stability_1970['theta transcription method'] = np.nan

df_phinney_stability_1970['D_10/d_0']                   = np.nan
df_phinney_stability_1970['D_10/d_0 stat error']        = np.nan
df_phinney_stability_1970['D_10/d_0 resolution']        = np.nan
df_phinney_stability_1970['D_30/d_0']                   = np.nan
df_phinney_stability_1970['D_30/d_0 stat error']        = np.nan
df_phinney_stability_1970['D_30/d_0 resolution']        = np.nan
df_phinney_stability_1970['D_32/d_0']                   = np.nan
df_phinney_stability_1970['D_32/d_0 stat error']        = np.nan
df_phinney_stability_1970['D_32/d_0 resolution']        = np.nan
df_phinney_stability_1970['droplet x/d_0']              = np.nan
df_phinney_stability_1970['D/d_0 page fig']             = np.nan
df_phinney_stability_1970['D/d_0 transcription method'] = np.nan
df_phinney_stability_1970['D/d_0 measurement method']   = np.nan
df_phinney_stability_1970['v_d_bar/vp']                        = np.nan
df_phinney_stability_1970['v_d_bar/vp stat error']             = np.nan
df_phinney_stability_1970['v_d_bar/vp resolution']             = np.nan
df_phinney_stability_1970['v_d_bar/vp page fig']               = np.nan
df_phinney_stability_1970['v_d_bar/vp transcription method']   = np.nan
df_phinney_stability_1970['e_p']                        = np.nan
df_phinney_stability_1970['e_p page fig']               = np.nan
df_phinney_stability_1970['e_p transcription method']   = np.nan

df_phinney_stability_1970['x_trans/d_0']                      = np.nan
df_phinney_stability_1970['x_trans/d_0 page fig']             = np.nan
df_phinney_stability_1970['x_trans/d_0 transcription method'] = np.nan

df_phinney_stability_1970['x_i/d_0']                      = np.nan
df_phinney_stability_1970['x_i/d_0 stat error']           = np.nan
df_phinney_stability_1970['x_i/d_0 resolution']           = np.nan
df_phinney_stability_1970['x_i/d_0 page fig']             = np.nan
df_phinney_stability_1970['x_i/d_0 transcription method'] = np.nan

df_phinney_stability_1970['x_e/d_0']                      = np.nan
df_phinney_stability_1970['x_e/d_0 stat error']           = np.nan
df_phinney_stability_1970['x_e/d_0 resolution']           = np.nan
df_phinney_stability_1970['x_e/d_0 page fig']             = np.nan
df_phinney_stability_1970['x_e/d_0 transcription method'] = np.nan

df_phinney_stability_1970['photo filename']   = np.nan
df_phinney_stability_1970['exposure time']    = np.nan
df_phinney_stability_1970['flash']            = np.nan
df_phinney_stability_1970['x_low']            = np.nan
df_phinney_stability_1970['x_mid']            = np.nan
df_phinney_stability_1970['x_high']           = np.nan
df_phinney_stability_1970['photo page fig']   = np.nan
df_phinney_stability_1970['spectrum']         = np.nan
df_phinney_stability_1970['lighting']         = np.nan
df_phinney_stability_1970['background color'] = np.nan
df_phinney_stability_1970['grid']             = np.nan
df_phinney_stability_1970['distance']         = np.nan
df_phinney_stability_1970['f-stop']           = np.nan
df_phinney_stability_1970['focal length']     = np.nan
df_phinney_stability_1970['camera model']     = np.nan
df_phinney_stability_1970['sensor']           = np.nan

df_phinney_stability_1970 = df_phinney_stability_1970[cols]
summary_table(df_phinney_stability_1970)
df_jet_breakup = pd.concat([df_jet_breakup, df_phinney_stability_1970])

########################
# phinney_breakup_1973 #
########################

# TODO: Add a column for data to remove from the dataframe before processing. Mark all fluid II data as to be removed, due to the error in the reporting of the fluid properties. Mark all inconsistent fluid I data to be removed. This should allow at least some data to improve the regime map and breakup length correlation.

# It appears that no data points are duplicated with phinney_stability_1970 or phinney_breakup_1975.

# p. 694: > The accuracies of the fluid properties and nozzle diameters are easily within 1%.

# DONE: Ask Dr. Phinney about the problems with this data.
# Issues with this data:
# 1. It is unclear how they reduced the surface tension so much for fluid II if it is merely salt water. Temperature? I am removing fluid II for now.
# 2. The stability parameter as defined in figure 2 on p. 692 appears to be wrong. However, I have verified that the correct definition is what was used when producing figure 3 by comparing against data from chen_mechanics_1962.

phinney_breakup_1973_csv = pd.read_csv('../data/phinney_breakup_1973/phinney_breakup_1973.csv', sep=',', header=0)

nozzle_phinney_breakup_1973     = phinney_breakup_1973_csv['nozzle']
fluid_phinney_breakup_1973      = phinney_breakup_1973_csv['fluid']
Re_l0_phinney_breakup_1973      = phinney_breakup_1973_csv['Re_l0']
lambda_phinney_breakup_1973     = phinney_breakup_1973_csv['lambda']
regime_L_b_phinney_breakup_1973 = phinney_breakup_1973_csv['regime L_b']

i = 0
d_0_phinney_breakup_1973    = np.zeros(len(Re_l0_phinney_breakup_1973))
L_0s_phinney_breakup_1973   = np.zeros(len(Re_l0_phinney_breakup_1973))
e_L_bs_phinney_breakup_1973 = np.zeros(len(Re_l0_phinney_breakup_1973))
for nozzle in nozzle_phinney_breakup_1973: # p. 694, table 2
   if nozzle == 'A':
      d_0_phinney_breakup_1973[i]    = 0.493e-2 # m
      L_0s_phinney_breakup_1973[i]   = 103.5
   elif nozzle == 'B':
      d_0_phinney_breakup_1973[i]    = 0.206e-2 # m
      L_0s_phinney_breakup_1973[i]   = 103.1
   elif nozzle == 'C':
      d_0_phinney_breakup_1973[i]    = 0.125e-2 # m
      L_0s_phinney_breakup_1973[i]   = 102.7
   elif nozzle == 'D':
      d_0_phinney_breakup_1973[i]    = 0.0504e-2 # m
      L_0s_phinney_breakup_1973[i]   = 103.5
   else:
      raise ValueError('Invalid nozzle for phinney_breakup_1973: ' + nozzle)
   
   e_L_bs_phinney_breakup_1973[i] = 0.2e-3 / d_0_phinney_breakup_1973[i] # phinney_stability_1970 p. 9: > The accuracy of position measurement was approximately .2 cm. ==> Or could use smaller value from phinney_breakup_1975 p. 997L
   i = i + 1

i = 0

# P = 1 atm according to table 1
rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g

We_l0_phinney_breakup_1973 = np.zeros(len(Re_l0_phinney_breakup_1973))
Fr_0_phinney_breakup_1973  = np.zeros(len(Re_l0_phinney_breakup_1973))
rho_s_phinney_breakup_1973 = np.zeros(len(Re_l0_phinney_breakup_1973))
nu_s_phinney_breakup_1973  = np.zeros(len(Re_l0_phinney_breakup_1973))
K_c_phinney_breakup_1973   = np.zeros(len(Re_l0_phinney_breakup_1973))
Ma_g_phinney_breakup_1973  = np.zeros(len(Re_l0_phinney_breakup_1973))
L_bs_phinney_breakup_1973  = np.zeros(len(Re_l0_phinney_breakup_1973))
#regime_L_b_phinney_breakup_1973  = []
regime_turb_phinney_breakup_1973 = []
liquid_phinney_breakup_1973      = []
for fluid in fluid_phinney_breakup_1973: # p. 693, table 1
   Re_l0 = Re_l0_phinney_breakup_1973[i]
   d_0   = d_0_phinney_breakup_1973[i]
   
   if fluid == 'I':
      rho_l = 1.106e3 # kg/m^3
      sigma = 77.5e-3 # N/m
      mu_l  = 1.34e-3 # Pa*s
      liquid_phinney_breakup_1973.append('salt water (I)') # p. 693
   elif fluid == 'II':
      rho_l = 1.2575e3 # kg/m^3
      sigma = 37.4e-3  # N/m
      mu_l  = 1.41e-3  # Pa*s
      liquid_phinney_breakup_1973.append('salt water (II)') # p. 693
   else:
      raise ValueError('Invalid fluid for phinney_breakup_1973: ' + fluid)
   
   We_l0  = (Re_l0 * mu_l)**2 / (rho_l * sigma * d_0)
   Ubar_0 = sqrt(We_l0 * sigma / (rho_l * d_0))
   Oh_l0  = sqrt(We_l0) / Re_l0
   
   nu_l = mu_l / rho_l
   P_v  = P_v_water(T_std) # TODO: Later change this to be for the correct fluid.
   
   We_l0_phinney_breakup_1973[i] = We_l0
   Fr_0_phinney_breakup_1973[i] = Ubar_0**2 / (g * d_0)
   rho_s_phinney_breakup_1973[i] = rho_l / rho_g
   nu_s_phinney_breakup_1973[i]  = nu_l / nu_g
   #K_c_phinney_breakup_1973[i]   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_phinney_breakup_1973[i]   = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_reentrant, friction_factor_smooth(Re_l0), L_0s_phinney_breakup_1973[i])
   Ma_g_phinney_breakup_1973[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)
   
   L_bs_phinney_breakup_1973[i] = lambda_phinney_breakup_1973[i] * sqrt(We_l0) * (1 + 3 * Oh_l0)
   
   # Wrong definition (but as in fig. 2):
   #L_bs_phinney_breakup_1973[i] = lambda_phinney_breakup_1973[i] * sqrt(We_l0 * (1 + 3 * Oh_l0))
   
   # if Re_l0 < 5.e3:
      # regime_L_b_phinney_breakup_1973.append('Rayleigh')
   # else:
      # regime_L_b_phinney_breakup_1973.append('R2S')
   # #else:
      # #regime_L_b_phinney_breakup_1973.append('second wind-induced')
   
   if Re_l0 < Re_trans:
      regime_turb_phinney_breakup_1973.append('laminar')
   elif Re_l0 < Re_turb:
      regime_turb_phinney_breakup_1973.append('transitional')
   else:
      regime_turb_phinney_breakup_1973.append('turbulent')
   
   # phinney_stability_1970 p. 8: > The density is measured [...] The accuracy is certainly within one percent.
   # phinney_stability_1970 p. 8: > > A series of readings on fresh samples were taken until it was felt that the value of $\sigma$ could be defined within $\pm.5$ percent.
   # phinney_stability_1970 p. 9: > For any given test (on one pipe) all values of $\mu$ would be within a band $\pm 1$ percent. For different pipes the values of \$mu$ would be within 3 percent.
   i = i + 1

I_0_phinney_breakup_1973        = I_fully_developed_array(Re_l0_phinney_breakup_1973, Re_turb)
Lambda_r_s_phinney_breakup_1973 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_phinney_breakup_1973, Re_turb), 'v')

# TODO: Figure out why one case is so much lower than the others, and why this data seems inconsistent with the others.
# A and B are consistent with each other, but not C and D.
# C and D are consistent with each other.
# The two fluids are inconsistent with each other.
#plt.loglog(We_l0_phinney_breakup_1973, L_bs_phinney_breakup_1973, 'x')
#plt.xlabel('We_{l0}')
#plt.ylabel('L_b/d_0')
#plt.grid()
#plt.show()

df_phinney_breakup_1973 = pd.DataFrame({'We_l0': We_l0_phinney_breakup_1973})
df_phinney_breakup_1973['key']     = 'phinney_breakup_1973'
df_phinney_breakup_1973['alt key'] = np.nan

#df_phinney_breakup_1973['liquid']     = 'salt water' # p. 693
df_phinney_breakup_1973['liquid']     = liquid_phinney_breakup_1973
df_phinney_breakup_1973['gas']        = 'air'
df_phinney_breakup_1973['degas']      = np.nan
df_phinney_breakup_1973['regulation'] = True # phinney_stability_1970 p. 6

df_phinney_breakup_1973['bend']              = False
df_phinney_breakup_1973['flow straightener'] = False
df_phinney_breakup_1973['screen']            = False
df_phinney_breakup_1973['contraction shape'] = 're-entrant' # phinney_stability_1970 fig. 3, pdf p. 53
df_phinney_breakup_1973['c']                 = 0.625*2.54e-2 / d_0_phinney_breakup_1973 # 0.625 inch diameter inlet pipe
df_phinney_breakup_1973['trip']              = False
df_phinney_breakup_1973['L_r/d_0']           = 0
df_phinney_breakup_1973['eps_r/d_0']         = np.nan

df_phinney_breakup_1973['L_0/d_0']       = L_0s_phinney_breakup_1973
df_phinney_breakup_1973['roughness']     = 'smooth'
df_phinney_breakup_1973['f page fig']    = np.nan
df_phinney_breakup_1973['pipe material'] = 'glass' # p. 694, table 2
df_phinney_breakup_1973['est f']         = True
df_phinney_breakup_1973['check FD']      = np.nan

df_phinney_breakup_1973['t/d_0']       = (0.25 * 2.54e-2 - d_0_phinney_breakup_1973) / (2 * d_0_phinney_breakup_1973) # fig. 3 says the nominal outer diameter of the glass tubes is 1/4"
df_phinney_breakup_1973['L_tip/d_0']   = np.nan
df_phinney_breakup_1973['end checked'] = True # phinney_stability_1970 p. 6: > Special care was taken to make sure the neds of the nozzles were square.

df_phinney_breakup_1973['orientation']         = 'horizontal'
df_phinney_breakup_1973['vibration isolation'] = True # phinney_stability_1970 p. 6: > To minimize disturbances...
df_phinney_breakup_1973['d_chamber/d_0']       = np.inf

#df_phinney_breakup_1973['We_l0']      = We_l0_phinney_breakup_1973
df_phinney_breakup_1973['Re_l0']      = Re_l0_phinney_breakup_1973
df_phinney_breakup_1973['Fr_0']       = Fr_0_phinney_breakup_1973
df_phinney_breakup_1973['I_0']        = I_0_phinney_breakup_1973
df_phinney_breakup_1973['rho_s']      = rho_s_phinney_breakup_1973
df_phinney_breakup_1973['nu_s']       = nu_s_phinney_breakup_1973
df_phinney_breakup_1973['K_c']        = K_c_phinney_breakup_1973
df_phinney_breakup_1973['Ma_g']       = Ma_g_phinney_breakup_1973
df_phinney_breakup_1973['Lambda_r_s'] = Lambda_r_s_phinney_breakup_1973
df_phinney_breakup_1973['est rho_s']  = True
df_phinney_breakup_1973['est nu_s']   = True

df_phinney_breakup_1973['regime photo'] = np.nan
df_phinney_breakup_1973['regime L_b']   = regime_L_b_phinney_breakup_1973
df_phinney_breakup_1973['regime turb']  = regime_turb_phinney_breakup_1973
df_phinney_breakup_1973['est turb']     = True

df_phinney_breakup_1973['L_b/d_0']                      = L_bs_phinney_breakup_1973
df_phinney_breakup_1973['L_b/d_0 stat error']           = np.nan
df_phinney_breakup_1973['L_b/d_0 resolution']           = e_L_bs_phinney_breakup_1973
df_phinney_breakup_1973['L_b method']                   = 'electrical'
df_phinney_breakup_1973['L_b/d_0 page fig']             = 'p. 695, fig. 3'
df_phinney_breakup_1973['L_b/d_0 transcription method'] = 'figure'

df_phinney_breakup_1973['theta']                      = np.nan
df_phinney_breakup_1973['theta stat error']           = np.nan
df_phinney_breakup_1973['theta resolution']           = np.nan
df_phinney_breakup_1973['theta page fig']             = np.nan
df_phinney_breakup_1973['theta transcription method'] = np.nan

df_phinney_breakup_1973['D_10/d_0']                   = np.nan
df_phinney_breakup_1973['D_10/d_0 stat error']        = np.nan
df_phinney_breakup_1973['D_10/d_0 resolution']        = np.nan
df_phinney_breakup_1973['D_30/d_0']                   = np.nan
df_phinney_breakup_1973['D_30/d_0 stat error']        = np.nan
df_phinney_breakup_1973['D_30/d_0 resolution']        = np.nan
df_phinney_breakup_1973['D_32/d_0']                   = np.nan
df_phinney_breakup_1973['D_32/d_0 stat error']        = np.nan
df_phinney_breakup_1973['D_32/d_0 resolution']        = np.nan
df_phinney_breakup_1973['droplet x/d_0']              = np.nan
df_phinney_breakup_1973['D/d_0 page fig']             = np.nan
df_phinney_breakup_1973['D/d_0 transcription method'] = np.nan
df_phinney_breakup_1973['D/d_0 measurement method']   = np.nan
df_phinney_breakup_1973['v_d_bar/vp']                        = np.nan
df_phinney_breakup_1973['v_d_bar/vp stat error']             = np.nan
df_phinney_breakup_1973['v_d_bar/vp resolution']             = np.nan
df_phinney_breakup_1973['v_d_bar/vp page fig']               = np.nan
df_phinney_breakup_1973['v_d_bar/vp transcription method']   = np.nan
df_phinney_breakup_1973['e_p']                        = np.nan
df_phinney_breakup_1973['e_p page fig']               = np.nan
df_phinney_breakup_1973['e_p transcription method']   = np.nan

df_phinney_breakup_1973['x_trans/d_0']                      = np.nan
df_phinney_breakup_1973['x_trans/d_0 page fig']             = np.nan
df_phinney_breakup_1973['x_trans/d_0 transcription method'] = np.nan

df_phinney_breakup_1973['x_i/d_0']                      = np.nan
df_phinney_breakup_1973['x_i/d_0 stat error']           = np.nan
df_phinney_breakup_1973['x_i/d_0 resolution']           = np.nan
df_phinney_breakup_1973['x_i/d_0 page fig']             = np.nan
df_phinney_breakup_1973['x_i/d_0 transcription method'] = np.nan

df_phinney_breakup_1973['x_e/d_0']                      = np.nan
df_phinney_breakup_1973['x_e/d_0 stat error']           = np.nan
df_phinney_breakup_1973['x_e/d_0 resolution']           = np.nan
df_phinney_breakup_1973['x_e/d_0 page fig']             = np.nan
df_phinney_breakup_1973['x_e/d_0 transcription method'] = np.nan

df_phinney_breakup_1973['photo filename']   = np.nan
df_phinney_breakup_1973['exposure time']    = np.nan
df_phinney_breakup_1973['flash']            = np.nan
df_phinney_breakup_1973['x_low']            = np.nan
df_phinney_breakup_1973['x_mid']            = np.nan
df_phinney_breakup_1973['x_high']           = np.nan
df_phinney_breakup_1973['photo page fig']   = np.nan
df_phinney_breakup_1973['spectrum']         = np.nan
df_phinney_breakup_1973['lighting']         = np.nan
df_phinney_breakup_1973['background color'] = np.nan
df_phinney_breakup_1973['grid']             = np.nan
df_phinney_breakup_1973['distance']         = np.nan
df_phinney_breakup_1973['f-stop']           = np.nan
df_phinney_breakup_1973['focal length']     = np.nan
df_phinney_breakup_1973['camera model']     = np.nan
df_phinney_breakup_1973['sensor']           = np.nan

df_phinney_breakup_1973 = df_phinney_breakup_1973[cols]
df_phinney_breakup_1973 = df_phinney_breakup_1973[df_phinney_breakup_1973['liquid'] == 'salt water (I)'] # Filtering out fluid II as the fluid properties are likely wrong, among other things.
summary_table(df_phinney_breakup_1973)
df_jet_breakup = pd.concat([df_jet_breakup, df_phinney_breakup_1973])

########################
# phinney_breakup_1975 #
########################

# Has turbulent breakup data for low pressure atmospheres.

# Issues with this data:
# 1. Very hard to figure out which point is which in figure 2. I called the author on the phone and they could not remember if there's a NOL report with tabulated data from this. I am going to try filing a FOIA for that.
# 2. 

# TODO: Ask Dr. Phinney about the air and saline temperature. Assuming T_std for now.

# p. 997L:
# > The inside diameter of the cylinder is 28 cm; the length is 120 cm.
# > The accuracy of the length measurements is about 1 mm.

# p. 997R
# glass tubing
# table salt and water
# rho_l = 1.062 g/cm^3    = 1062 kg/cm^3
# mu_l  = 1.15 centipoise = 1.15e-3 Pa*s
# sigma = 73.1 dyne/cm    = 73.1e-3 N/m

phinney_breakup_1975_csv = pd.read_csv('../data/phinney_breakup_1975/phinney_breakup_1975_fig_2.csv', sep=',', header=0)

nozzle_phinney_breakup_1975 = phinney_breakup_1975_csv['nozzle']
P_atm_phinney_breakup_1975  = phinney_breakup_1975_csv['P_atm (mm)'] * 133.322 # Pa
We_l0_phinney_breakup_1975  = phinney_breakup_1975_csv['We_l0']
L_bs_phinney_breakup_1975   = phinney_breakup_1975_csv['L_b/d_0']

i = 0
d_0_phinney_breakup_1975    = np.zeros(len(We_l0_phinney_breakup_1975))
L_0s_phinney_breakup_1975   = np.zeros(len(We_l0_phinney_breakup_1975))
e_L_bs_phinney_breakup_1975 = np.zeros(len(We_l0_phinney_breakup_1975))
for nozzle in nozzle_phinney_breakup_1975:
   # p. 997R, table 1
   if nozzle == 'B':
      d_0_phinney_breakup_1975[i]    = 0.206e-2
      L_0s_phinney_breakup_1975[i]   = 20.3e-2 / d_0_phinney_breakup_1975[i]
   elif nozzle == 'G':
      d_0_phinney_breakup_1975[i]  = 0.1024e-2
      L_0s_phinney_breakup_1975[i] = 10.5e-2 / d_0_phinney_breakup_1975[i]
   else:
      raise ValueError('Invalid nozzle for phinney_breakup_1975.')
   
   e_L_bs_phinney_breakup_1975[i] = 2.e-3 / d_0_phinney_breakup_1975[i] # Not using 1 mm from phinney_breakup_1975 p. 997L. Using 2 mm from other Phinney studies.
   i = i + 1

i = 0

# p. 997R: table salt and water
rho_l = 1.062e3  # kg/m^3
sigma = 73.1e-3  # N/m
mu_l  = 1.15e-3  # Pa*s
nu_l = mu_l / rho_l
P_v  = P_v_water(T_std) # TODO: Later change this to be for the correct fluid. This is for water when it actually needs salt water.

Re_l0_phinney_breakup_1975      = np.zeros(len(We_l0_phinney_breakup_1975))
Fr_0_phinney_breakup_1975       = np.zeros(len(We_l0_phinney_breakup_1975))
rho_s_phinney_breakup_1975      = np.zeros(len(We_l0_phinney_breakup_1975))
nu_s_phinney_breakup_1975       = np.zeros(len(We_l0_phinney_breakup_1975))
K_c_phinney_breakup_1975        = np.zeros(len(We_l0_phinney_breakup_1975))
Ma_g_phinney_breakup_1975       = np.zeros(len(We_l0_phinney_breakup_1975))
regime_L_b_phinney_breakup_1975  = []
regime_turb_phinney_breakup_1975 = []
for We_l0 in We_l0_phinney_breakup_1975:
   d_0    = d_0_phinney_breakup_1975[i]
   Ubar_0 = sqrt(We_l0 * sigma / (rho_l * d_0))
   
   rho_g = rho_ideal_gas(P_atm_phinney_breakup_1975[i], T_std, MW_air)
   nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
   
   Re_l0_phinney_breakup_1975[i] = rho_l * Ubar_0 * d_0 / mu_l
   Fr_0_phinney_breakup_1975[i] = Ubar_0**2 / (g * d_0)
   rho_s_phinney_breakup_1975[i] = rho_l / rho_g
   nu_s_phinney_breakup_1975[i]  = nu_l / nu_g
   #K_c_phinney_breakup_1975[i]   = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_phinney_breakup_1975[i]   = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_reentrant, friction_factor_smooth(Re_l0_phinney_breakup_1975[i]), L_0s_phinney_breakup_1975[i])
   Ma_g_phinney_breakup_1975[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air)
   
   if We_l0 > 3000.: # where does this number come from?
      regime_L_b_phinney_breakup_1975.append('second wind-induced')
   elif We_l0 > 1000.: # Set 2019-06-27.
      regime_L_b_phinney_breakup_1975.append('R2S')
   else:
      regime_L_b_phinney_breakup_1975.append('Rayleigh')
   
   if Re_l0 > Re_turb:
      regime_turb_phinney_breakup_1975.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_phinney_breakup_1975.append('transitional')
   else:
      regime_turb_phinney_breakup_1975.append('laminar')
   
   # phinney_stability_1970 p. 8: > The density is measured [...] The accuracy is certainly within one percent.
   # phinney_stability_1970 p. 8: > > A series of readings on fresh samples were taken until it was felt that the value of $\sigma$ could be defined within $\pm.5$ percent.
   # phinney_stability_1970 p. 9: > For any given test (on one pipe) all values of $\mu$ would be within a band $\pm 1$ percent. For different pipes the values of \$mu$ would be within 3 percent.
   i = i + 1

I_0_phinney_breakup_1975 = I_fully_developed_array(Re_l0_phinney_breakup_1975, Re_turb)
Lambda_r_s_phinney_breakup_1975 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_phinney_breakup_1975, Re_turb), 'v')

df_phinney_breakup_1975 = pd.DataFrame({'We_l0': We_l0_phinney_breakup_1975})
df_phinney_breakup_1975['key']        = 'phinney_breakup_1975'
df_phinney_breakup_1975['alt key']    = np.nan

df_phinney_breakup_1975['liquid']     = 'salt water'
df_phinney_breakup_1975['gas']        = 'air'
df_phinney_breakup_1975['degas']      = np.nan # TODO
df_phinney_breakup_1975['regulation'] = True # phinney_stability_1970 p. 6

df_phinney_breakup_1975['bend']              = False
df_phinney_breakup_1975['flow straightener'] = False
df_phinney_breakup_1975['screen']            = False
df_phinney_breakup_1975['contraction shape'] = 're-entrant' # phinney_stability_1970 fig. 3, pdf p. 53
df_phinney_breakup_1975['c']                 = 0.625*2.54e-2 / d_0_phinney_breakup_1975 # 0.625 inch diameter inlet pipe
df_phinney_breakup_1975['trip']              = False
df_phinney_breakup_1975['L_r/d_0']           = 0
df_phinney_breakup_1975['eps_r/d_0']         = np.nan

df_phinney_breakup_1975['L_0/d_0']       = L_0s_phinney_breakup_1975
df_phinney_breakup_1975['roughness']     = 'smooth'
df_phinney_breakup_1975['f page fig']    = np.nan
df_phinney_breakup_1975['pipe material'] = 'glass'
df_phinney_breakup_1975['est f']         = True
df_phinney_breakup_1975['check FD']      = np.nan

df_phinney_breakup_1975['t/d_0']       = (0.25 * 2.54e-2 - d_0_phinney_breakup_1975) / (2 * d_0_phinney_breakup_1975) # phinney_stability_1970 fig. 3 says the nominal outer diameter of the glass tubes is 1/4"
df_phinney_breakup_1975['L_tip/d_0']   = 0
df_phinney_breakup_1975['end checked'] = True

df_phinney_breakup_1975['orientation']         = 'horizontal'
df_phinney_breakup_1975['vibration isolation'] = True # phinney_stability_1970 p. 6: > To minimize disturbances...
df_phinney_breakup_1975['d_chamber/d_0']       = 28e-2 / d_0_phinney_breakup_1975 # p. 997L

df_phinney_breakup_1975['L_b/d_0']                      = L_bs_phinney_breakup_1975
df_phinney_breakup_1975['L_b/d_0 page fig']             = 'p. 998, fig. 2'
df_phinney_breakup_1975['L_b/d_0 stat error']           = np.nan
df_phinney_breakup_1975['L_b/d_0 resolution']           = e_L_bs_phinney_breakup_1975
df_phinney_breakup_1975['L_b method']                   = 'electrical'
df_phinney_breakup_1975['L_b/d_0 transcription method'] = 'figure'

#df_phinney_breakup_1975['We_l0']      = We_l0_phinney_breakup_1975
df_phinney_breakup_1975['Re_l0']      = Re_l0_phinney_breakup_1975
df_phinney_breakup_1975['I_0']        = I_0_phinney_breakup_1975
df_phinney_breakup_1975['Fr_0']       = Fr_0_phinney_breakup_1975
df_phinney_breakup_1975['rho_s']      = rho_s_phinney_breakup_1975
df_phinney_breakup_1975['nu_s']       = nu_s_phinney_breakup_1975
df_phinney_breakup_1975['K_c']        = K_c_phinney_breakup_1975
df_phinney_breakup_1975['Ma_g']       = Ma_g_phinney_breakup_1975
df_phinney_breakup_1975['Lambda_r_s'] = Lambda_r_s_phinney_breakup_1975
df_phinney_breakup_1975['est rho_s']  = True
df_phinney_breakup_1975['est nu_s']   = True

df_phinney_breakup_1975['regime photo'] = np.nan
df_phinney_breakup_1975['regime L_b']   = regime_L_b_phinney_breakup_1975
df_phinney_breakup_1975['regime turb']  = regime_turb_phinney_breakup_1975
df_phinney_breakup_1975['est turb']     = True

df_phinney_breakup_1975['theta']                      = np.nan
df_phinney_breakup_1975['theta stat error']           = np.nan
df_phinney_breakup_1975['theta resolution']           = np.nan
df_phinney_breakup_1975['theta page fig']             = np.nan
df_phinney_breakup_1975['theta transcription method'] = np.nan

df_phinney_breakup_1975['D_10/d_0']                   = np.nan
df_phinney_breakup_1975['D_10/d_0 stat error']        = np.nan
df_phinney_breakup_1975['D_10/d_0 resolution']        = np.nan
df_phinney_breakup_1975['D_30/d_0']                   = np.nan
df_phinney_breakup_1975['D_30/d_0 stat error']        = np.nan
df_phinney_breakup_1975['D_30/d_0 resolution']        = np.nan
df_phinney_breakup_1975['D_32/d_0']                   = np.nan
df_phinney_breakup_1975['D_32/d_0 stat error']        = np.nan
df_phinney_breakup_1975['D_32/d_0 resolution']        = np.nan
df_phinney_breakup_1975['droplet x/d_0']              = np.nan
df_phinney_breakup_1975['D/d_0 page fig']             = np.nan
df_phinney_breakup_1975['D/d_0 transcription method'] = np.nan
df_phinney_breakup_1975['D/d_0 measurement method']   = np.nan
df_phinney_breakup_1975['v_d_bar/vp']                        = np.nan
df_phinney_breakup_1975['v_d_bar/vp stat error']             = np.nan
df_phinney_breakup_1975['v_d_bar/vp resolution']             = np.nan
df_phinney_breakup_1975['v_d_bar/vp page fig']               = np.nan
df_phinney_breakup_1975['v_d_bar/vp transcription method']   = np.nan
df_phinney_breakup_1975['e_p']                        = np.nan
df_phinney_breakup_1975['e_p page fig']               = np.nan
df_phinney_breakup_1975['e_p transcription method']   = np.nan

df_phinney_breakup_1975['x_trans/d_0']                      = np.nan
df_phinney_breakup_1975['x_trans/d_0 page fig']             = np.nan
df_phinney_breakup_1975['x_trans/d_0 transcription method'] = np.nan

df_phinney_breakup_1975['x_i/d_0']                      = np.nan
df_phinney_breakup_1975['x_i/d_0 stat error']           = np.nan
df_phinney_breakup_1975['x_i/d_0 resolution']           = np.nan
df_phinney_breakup_1975['x_i/d_0 page fig']             = np.nan
df_phinney_breakup_1975['x_i/d_0 transcription method'] = np.nan

df_phinney_breakup_1975['x_e/d_0']                      = np.nan
df_phinney_breakup_1975['x_e/d_0 stat error']           = np.nan
df_phinney_breakup_1975['x_e/d_0 resolution']           = np.nan
df_phinney_breakup_1975['x_e/d_0 page fig']             = np.nan
df_phinney_breakup_1975['x_e/d_0 transcription method'] = np.nan

df_phinney_breakup_1975['photo filename']   = np.nan
df_phinney_breakup_1975['exposure time']    = np.nan
df_phinney_breakup_1975['flash']            = np.nan
df_phinney_breakup_1975['x_low']            = np.nan
df_phinney_breakup_1975['x_mid']            = np.nan
df_phinney_breakup_1975['x_high']           = np.nan
df_phinney_breakup_1975['photo page fig']   = np.nan
df_phinney_breakup_1975['spectrum']         = np.nan
df_phinney_breakup_1975['lighting']         = np.nan
df_phinney_breakup_1975['background color'] = np.nan
df_phinney_breakup_1975['grid']             = np.nan
df_phinney_breakup_1975['distance']         = np.nan
df_phinney_breakup_1975['f-stop']           = np.nan
df_phinney_breakup_1975['focal length']     = np.nan
df_phinney_breakup_1975['camera model']     = np.nan
df_phinney_breakup_1975['sensor']           = np.nan

df_phinney_breakup_1975 = df_phinney_breakup_1975[cols]
summary_table(df_phinney_breakup_1975)
df_jet_breakup = pd.concat([df_jet_breakup, df_phinney_breakup_1975])

#############################
# sterling_instability_1975 #
#############################

# TODO: Add tabulated data from sterling_instability_1969.
# TODO: Add wavelength data (sterling_instability_1969 table VII). This has a different camera setup, so you should change the camera section for these. See sterling_instability_1975 p. 488.

sterling_instability_1975_csv = pd.read_csv('../data/sterling_instability_1975/sterling_instability_1969.csv', sep=',', header=0)

temp_sterling_instability_1975  = 23. # C, from sterling_instability_1975 table 1 (assumed gas temperature is the same as the fluid temperature)
P_atm_sterling_instability_1975 = P_atm # sterling_instability_1969 p. 81 just says "1 atm".

liquid_sterling_instability_1975     = sterling_instability_1975_csv['liquid']
sigma_sterling_instability_1975      = sterling_instability_1975_csv['surface tension (dynes/cm)'] * 1.e-3 # N/m
rho_l_sterling_instability_1975      = sterling_instability_1975_csv['density (g/cc)'] * 1.e3 # kg/m^3
mu_l_sterling_instability_1975       = sterling_instability_1975_csv['viscosity (cp)'] * 1.e-3 # N*s/m^2
d_0_sterling_instability_1975        = sterling_instability_1975_csv['d_0 (cm)'] * 1.e-2 # m
L_0s_sterling_instability_1975       = sterling_instability_1975_csv['L_0/d_0']
Ubar_0_sterling_instability_1975     = sterling_instability_1975_csv['Ubar_0 (cm/s)'] * 1.e-2 # m/s
xbavg_sterling_instability_1975      = sterling_instability_1975_csv['xbavg (cm)'] * 1.e-2 # m
table_sterling_instability_1975      = sterling_instability_1975_csv['table']
page_sterling_instability_1975       = sterling_instability_1975_csv['page']
regime_L_b_sterling_instability_1975 = sterling_instability_1975_csv['regime L_b']

xbavg_s_sterling_instability_1975 = xbavg_sterling_instability_1975 / d_0_sterling_instability_1975

rho_g_sterling_instability_1975 = rho_ideal_gas(P_atm_sterling_instability_1975, temp_sterling_instability_1975, MW_air)
nu_g_sterling_instability_1975  = mu_ideal_gas(temp_sterling_instability_1975, mu_0_air, T_0_air, C_air) / rho_g_sterling_instability_1975

nu_l_sterling_instability_1975 = mu_l_sterling_instability_1975 / rho_l_sterling_instability_1975
P_v_sterling_instability_1975  = P_v_water(temp_sterling_instability_1975)

We_l0_sterling_instability_1975      = rho_l_sterling_instability_1975 * Ubar_0_sterling_instability_1975**2. * d_0_sterling_instability_1975 / sigma_sterling_instability_1975
Re_l0_sterling_instability_1975      = Ubar_0_sterling_instability_1975 * d_0_sterling_instability_1975 / nu_l_sterling_instability_1975
f_sterling_instability_1975 = friction_factor_smooth_array(Re_l0_sterling_instability_1975)
Fr_0_sterling_instability_1975       = Ubar_0_sterling_instability_1975**2. / (d_0_sterling_instability_1975 * g)
rho_s_sterling_instability_1975      = rho_l_sterling_instability_1975 / rho_g_sterling_instability_1975
nu_s_sterling_instability_1975       = nu_l_sterling_instability_1975 / nu_g_sterling_instability_1975
K_c_sterling_instability_1975        = K_c(P_atm_sterling_instability_1975, P_v_sterling_instability_1975, rho_l_sterling_instability_1975, Ubar_0_sterling_instability_1975, K_L_sharpedged, friction_factor_smooth_array(Re_l0_sterling_instability_1975), L_0s_sterling_instability_1975)

Ma_g_sterling_instability_1975        = []
Lambda_r_s_sterling_instability_1975  = []
regime_turb_sterling_instability_1975 = []
I_0_sterling_instability_1975         = []
for Ubar_0, Re_l0, L_0s in zip(Ubar_0_sterling_instability_1975, Re_l0_sterling_instability_1975, L_0s_sterling_instability_1975):
   #print Re_l0
   assert(Re_l0 > 10.)
   if L_0s > 30.:
      f = friction_factor_smooth(Re_l0)
      Ma_g_sterling_instability_1975.append(Ubar_0 / c_ideal_gas(gamma_air(temp_sterling_instability_1975), temp_sterling_instability_1975 + T_zero, MW_air))
      
      if Re_l0 > Re_turb:
         regime_turb_sterling_instability_1975.append('turbulent')
         Lambda_r_s_sterling_instability_1975.append(Lambda_s('smooth', f, 'v'))
         I_0_sterling_instability_1975.append(I_fully_developed(f))
         
         assert(I_fully_developed(f) < 0.15)
      elif Re_l0 > Re_trans:
         regime_turb_sterling_instability_1975.append('transitional')
         Lambda_r_s_sterling_instability_1975.append(np.nan)
         I_0_sterling_instability_1975.append(np.nan)
      else:
         regime_turb_sterling_instability_1975.append('laminar')
         Lambda_r_s_sterling_instability_1975.append(np.nan)
         I_0_sterling_instability_1975.append(np.nan)
   else:
      Ma_g_sterling_instability_1975.append(np.nan)
      Lambda_r_s_sterling_instability_1975.append(np.nan)
      I_0_sterling_instability_1975.append(np.nan)
      regime_turb_sterling_instability_1975.append(np.nan)

df_sterling_instability_1975 = pd.DataFrame({'key': 'sterling_instability_1975',
                                       'We_l0': We_l0_sterling_instability_1975})
df_sterling_instability_1975['liquid']     = liquid_sterling_instability_1975
df_sterling_instability_1975['alt key']    = 'sterling_instability_1969'
df_sterling_instability_1975['gas']        = 'air'
df_sterling_instability_1975['degas']      = False
df_sterling_instability_1975['regulation'] = True # sterling_instability_1969 p. 90

df_sterling_instability_1975['bend']              = False
df_sterling_instability_1975['flow straightener'] = True
df_sterling_instability_1975['screen']            = True
df_sterling_instability_1975['contraction shape'] = 'smooth'
df_sterling_instability_1975['c']                 = 1.5 * 2.54e-2 / d_0_sterling_instability_1975 # fig. 3
df_sterling_instability_1975['trip']              = False
df_sterling_instability_1975['L_r/d_0']           = 0
df_sterling_instability_1975['eps_r/d_0']         = np.nan

df_sterling_instability_1975['L_0/d_0']       = L_0s_sterling_instability_1975
df_sterling_instability_1975['roughness']     = 'smooth'
df_sterling_instability_1975['f page fig']    = np.nan
df_sterling_instability_1975['pipe material'] = 'steel'
df_sterling_instability_1975['est f']         = True
df_sterling_instability_1975['check FD']      = False

df_sterling_instability_1975['t/d_0']       = np.nan
df_sterling_instability_1975['L_tip/d_0']   = np.nan
df_sterling_instability_1975['end checked'] = np.nan

df_sterling_instability_1975['orientation']         = 'horizontal' # fig. 1
df_sterling_instability_1975['vibration isolation'] = True # fig. 1
df_sterling_instability_1975['d_chamber/d_0']       = np.nan

#df_sterling_instability_1975['We_l0']      = We_l0_sterling_instability_1975
df_sterling_instability_1975['Re_l0']      = Re_l0_sterling_instability_1975
df_sterling_instability_1975['I_0']        = I_0_sterling_instability_1975
df_sterling_instability_1975['Fr_0']       = Fr_0_sterling_instability_1975
df_sterling_instability_1975['rho_s']      = rho_s_sterling_instability_1975
df_sterling_instability_1975['nu_s']       = nu_s_sterling_instability_1975
df_sterling_instability_1975['K_c']        = K_c_sterling_instability_1975
df_sterling_instability_1975['Ma_g']       = Ma_g_sterling_instability_1975
df_sterling_instability_1975['Lambda_r_s'] = Lambda_r_s_sterling_instability_1975
df_sterling_instability_1975['est rho_s']  = True
df_sterling_instability_1975['est nu_s']   = True

df_sterling_instability_1975['regime photo'] = np.nan
df_sterling_instability_1975['regime L_b']   = regime_L_b_sterling_instability_1975
df_sterling_instability_1975['regime turb']  = regime_turb_sterling_instability_1975
df_sterling_instability_1975['est turb']     = True

df_sterling_instability_1975['L_b/d_0']                      = xbavg_s_sterling_instability_1975
df_sterling_instability_1975['L_b/d_0 stat error']           = np.nan
df_sterling_instability_1975['L_b/d_0 resolution']           = np.nan
df_sterling_instability_1975['L_b method']                   = 'photographic'
df_sterling_instability_1975['L_b/d_0 page fig']             = np.nan
df_sterling_instability_1975['L_b/d_0 transcription method'] = np.nan

df_sterling_instability_1975['theta']                      = np.nan
df_sterling_instability_1975['theta stat error']           = np.nan
df_sterling_instability_1975['theta resolution']           = np.nan
df_sterling_instability_1975['theta page fig']             = np.nan
df_sterling_instability_1975['theta transcription method'] = np.nan

df_sterling_instability_1975['D_10/d_0']                        = np.nan
df_sterling_instability_1975['D_10/d_0 stat error']             = np.nan
df_sterling_instability_1975['D_10/d_0 resolution']             = np.nan
df_sterling_instability_1975['D_32/d_0']                        = np.nan
df_sterling_instability_1975['D_32/d_0 stat error']             = np.nan
df_sterling_instability_1975['D_32/d_0 resolution']             = np.nan
df_sterling_instability_1975['D_30/d_0']                        = np.nan
df_sterling_instability_1975['D_30/d_0 stat error']             = np.nan
df_sterling_instability_1975['D_30/d_0 resolution']             = np.nan
df_sterling_instability_1975['droplet x/d_0']                   = np.nan
df_sterling_instability_1975['D/d_0 page fig']                  = np.nan
df_sterling_instability_1975['D/d_0 transcription method']      = np.nan
df_sterling_instability_1975['D/d_0 measurement method']        = np.nan
df_sterling_instability_1975['v_d_bar/vp']                      = np.nan
df_sterling_instability_1975['v_d_bar/vp stat error']           = np.nan
df_sterling_instability_1975['v_d_bar/vp resolution']           = np.nan
df_sterling_instability_1975['v_d_bar/vp page fig']             = np.nan
df_sterling_instability_1975['v_d_bar/vp transcription method'] = np.nan
df_sterling_instability_1975['e_p']                             = np.nan
df_sterling_instability_1975['e_p page fig']                    = np.nan
df_sterling_instability_1975['e_p transcription method']        = np.nan

df_sterling_instability_1975['x_trans/d_0']                      = np.nan
df_sterling_instability_1975['x_trans/d_0 page fig']             = np.nan
df_sterling_instability_1975['x_trans/d_0 transcription method'] = np.nan

df_sterling_instability_1975['x_i/d_0']                      = np.nan
df_sterling_instability_1975['x_i/d_0 stat error']           = np.nan
df_sterling_instability_1975['x_i/d_0 resolution']           = np.nan
df_sterling_instability_1975['x_i/d_0 page fig']             = np.nan
df_sterling_instability_1975['x_i/d_0 transcription method'] = np.nan

df_sterling_instability_1975['x_e/d_0']                      = np.nan
df_sterling_instability_1975['x_e/d_0 stat error']           = np.nan
df_sterling_instability_1975['x_e/d_0 resolution']           = np.nan
df_sterling_instability_1975['x_e/d_0 page fig']             = np.nan
df_sterling_instability_1975['x_e/d_0 transcription method'] = np.nan

df_sterling_instability_1975['photo filename']   = np.nan
df_sterling_instability_1975['exposure time']    = 0.8e-6 # s, sterling_instability_1969 p. 92
df_sterling_instability_1975['flash']            = True
df_sterling_instability_1975['x_low']            = np.nan
df_sterling_instability_1975['x_mid']            = np.nan
df_sterling_instability_1975['x_high']           = np.nan
df_sterling_instability_1975['photo page fig']   = np.nan
df_sterling_instability_1975['spectrum']         = 'visible'
df_sterling_instability_1975['lighting']         = 'stroboscope'
df_sterling_instability_1975['background color'] = np.nan
df_sterling_instability_1975['grid']             = np.nan
df_sterling_instability_1975['distance']         = np.nan
df_sterling_instability_1975['f-stop']           = np.nan
df_sterling_instability_1975['focal length']     = np.nan
df_sterling_instability_1975['camera model']     = '16 mm Arriflex motion-picture camera' # p. 487
df_sterling_instability_1975['sensor']           = '16 mm Plus-X film' # p. 487

df_sterling_instability_1975 = df_sterling_instability_1975[cols]
summary_table(df_sterling_instability_1975)

# remove all short nozzle data
df_sterling_instability_1975 = df_sterling_instability_1975[df_sterling_instability_1975['L_0/d_0'] > 30.]
#print df_sterling_instability_1975['L_0/d_0']
#print df_sterling_instability_1975['I_0']
df_jet_breakup = pd.concat([df_jet_breakup, df_sterling_instability_1975])

##########################
# reitz_atomization_1978 #
##########################

# S: steady
# T: transient

# Nozzles I and II. See tab. 3.1 p. 159.
# photo system details: pp. 65-66 (pdf pp. 89-90), 79 (pdf p. 103)

# pdf p. 121

# p. 78 (pdf p. 102): > intact jet length $x_I$. This length is defined as the distance from the nozzle exit to the point where jet divergence is observed to begin.

# TODO: Figure 4.10, p. 202
# TODO: Figure 5.3, p. 220 (pdf p. 244)
# TODO: Figure 5.4, p. 221 (pdf p. 245)
# TODO: Figure 5.6, p. 223 (pdf p. 247)
# TODO: Figure 5.8, p. 225 (pdf p. 249)
# TODO: Figure 5.10, p. 227 (pdf p. 251)
# TODO: Figure 5.12, p. 229 (pdf p. 253)
# TODO: Figure 5.13, p. 230 (pdf p. 254): cavitation number vs. theta
# TODO: Figure 5.14, p. 231 (pdf p. 255): gas jet theta
# TODO: Figure C-35, p. C38 (pdf p. 320): photos
# TODO: Figure C-41, p. C45 (pdf p. 327): photo (1bS?); p. 68 suggests this might have a polymer added
# TODO: reitz_dependence_1979 p. 12 fig. 11a

reitz_atomization_1978_csv = pd.read_csv('../data/reitz_atomization_1978/reitz_atomization_1978.csv', sep=',', header=0)

#series_reitz_atomization_1978         = reitz_atomization_1978_csv['series']
nozzle_reitz_atomization_1978         = reitz_atomization_1978_csv['nozzle']
rho_g_reitz_atomization_1978          = reitz_atomization_1978_csv['rho_g (kg/m^3)']
P_atm_reitz_atomization_1978          = reitz_atomization_1978_csv['P_atm (psia)'] * 6894.76 # Pa
gas_reitz_atomization_1978            = reitz_atomization_1978_csv['gas']
liquid_reitz_atomization_1978         = reitz_atomization_1978_csv['liquid']
K_reitz_atomization_1978              = reitz_atomization_1978_csv['K']
Re_l0_reitz_atomization_1978          = reitz_atomization_1978_csv['Re_l0']
We_g0_reitz_atomization_1978          = reitz_atomization_1978_csv['We_g0']
We_l0_reitz_atomization_1978          = reitz_atomization_1978_csv['We_l0']
Oh_l0_reitz_atomization_1978          = reitz_atomization_1978_csv['Oh_l0']
theta_reitz_atomization_1978          = (2 * np.pi / 360) * reitz_atomization_1978_csv['theta (deg.)'] # full cone angle as seen in reitz_dependence_1979 p. 8, fig. 7
x_is_reitz_atomization_1978           = reitz_atomization_1978_csv['x_i/d_0']
photo_filename_reitz_atomization_1978 = reitz_atomization_1978_csv['photo filename']

# p. 56 (pdf p. 80): > The details of the construction of nozzle I are as follows: it comprised a 2.92 cm length of stainless steel hypodermic tubing which had an inside diameter of 0.034 cm and an outside diameter of 0.066 cm and was silver soldered flush with the nozzle exit plane into a 0.066 cm centralled drilled hole in the brass plug. Its upper end (the nozzle inlet) communicated with a 0.79 cm diameter centrally drilled hole 1.2 cm deep whose bottom had been radiused with a 0.38 cm radius cutter.

i = 0
d_0_reitz_atomization_1978  = np.zeros(len(theta_reitz_atomization_1978))
L_0s_reitz_atomization_1978 = np.zeros(len(theta_reitz_atomization_1978))
ts_reitz_atomization_1978   = np.zeros(len(theta_reitz_atomization_1978))
contraction_shape_reitz_atomization_1978 = []
for nozzle in nozzle_reitz_atomization_1978: # see table 3.1 on p. 159 (pdf p. 183)
   if nozzle == 'I':
      d_0_reitz_atomization_1978[i]  = 0.34e-3 # m, p. 186
      L_0s_reitz_atomization_1978[i] = 85
      ts_reitz_atomization_1978[i]   = np.inf # p. 186 (fig. 3.14); assuming the other part is much wider
      contraction_shape_reitz_atomization_1978.append('sudden') # p. 186, fig. 3.14
   elif nozzle == 'II':
      d_0_reitz_atomization_1978[i]  = 0.34e-3 # m, p. 186
      L_0s_reitz_atomization_1978[i] = 49.3 # more accurate number from p. 227
      ts_reitz_atomization_1978[i]   = 0.66e-3 / (2 * d_0_reitz_atomization_1978[i]) - 1 # p. 186 (fig. 3.14)
      contraction_shape_reitz_atomization_1978.append('converging cone') # p. 186, fig. 3.14
   else:
      raise ValueError('Invalid nozzle for reitz_atomization_1978:', nozzle)
   
   i = i + 1

i = 0
Fr_0_reitz_atomization_1978       = np.zeros(len(theta_reitz_atomization_1978))
rho_s_reitz_atomization_1978      = np.zeros(len(theta_reitz_atomization_1978))
nu_s_reitz_atomization_1978       = np.zeros(len(theta_reitz_atomization_1978))
K_c_reitz_atomization_1978        = np.zeros(len(theta_reitz_atomization_1978))
Ma_g_reitz_atomization_1978       = np.zeros(len(theta_reitz_atomization_1978))
x_low_reitz_atomization_1978      = np.zeros(len(theta_reitz_atomization_1978))
x_is_res_reitz_atomization_1978   = np.zeros(len(theta_reitz_atomization_1978))
theta_res_reitz_atomization_1978  = np.zeros(len(theta_reitz_atomization_1978))
photo_page_fig_reitz_atomization_1978   = []
spectrum_reitz_atomization_1978         = []
background_color_reitz_atomization_1978 = []
grid_reitz_atomization_1978             = []
regime_turb_reitz_atomization_1978      = []
for rho_g in rho_g_reitz_atomization_1978:
   d_0   = d_0_reitz_atomization_1978[i]
   #rho_g = rho_g_reitz_atomization_1978[i]
   Re_l0 = Re_l0_reitz_atomization_1978[i]
   We_l0 = We_l0_reitz_atomization_1978[i]
   
   # reitz_atomization_1978 p. 162 (pdf p. 186) has liquid and gas property data. Presumably these were the values used.
   # reitz_dependence_1979 p. 5R: > The chamber gases were [...] and the tests were made at about 300$^\circ$K.
   
   gas    = gas_reitz_atomization_1978[i]
   liquid = liquid_reitz_atomization_1978[i]
   
   # TODO: Fix problem with gas viscosity?
   if gas == 'air':
      #rho_g = 1.205 # kg/m^3
      mu_g  = 1.83e-5 # Pa*s
      nu_g  = mu_g / rho_g
      #nu_g = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
      c_g   = 343 # m/s
   elif gas == 'N2':
      #rho_g = 1.165 # kg/m^3
      mu_g  = 1.70e-5 # Pa*s
      nu_g  = mu_g / rho_g
      #nu_g = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
      c_g   = 349 # m/s
   else:
      raise ValueError('Invalid gas for reitz_atomization_1978:', gas)
   
   if liquid == 'water':
      rho_l = 998 # kg/m^3
      mu_l  = 0.001 # Pa*s
      nu_l  = mu_l / rho_l
      sigma = 72.8e-3 # N/m
      P_v   = 0.36 * 6894.76 # Pa
   elif liquid == '50% glycerol-water':
      rho_l = 1125 # kg/m^3
      mu_l  = 0.006 # Pa*s
      nu_l  = mu_l / rho_l
      sigma = 70.0e-3 # N/m
      P_v   = 0.5 * (0.36 + 2.0e-6) * 6894.76 # Pa
      # TODO: Get a better vapor pressure estimate later.
   else:
      raise ValueError('Invalid liquid for reitz_atomization_1978:', liquid)
   
   #Ubar_0 = (Re_l0 * mu_l) / (rho_l * d_0)
   Ubar_0 = sqrt((We_l0 * sigma) / (rho_l * d_0))
   #print Ubar_0
   
   Fr_0_reitz_atomization_1978[i]  = Ubar_0**2 / (g * d_0)
   rho_s_reitz_atomization_1978[i] = rho_l / rho_g
   nu_s_reitz_atomization_1978[i]  = nu_l / nu_g
   #K_c_reitz_atomization_1978[i]   = (P_atm_reitz_atomization_1978[i] - P_v) / (0.5 * rho_l * Ubar_0**2)
   if contraction_shape_reitz_atomization_1978[i] == 'sudden':
      K_L_reitz = K_L_sharpedged
   elif contraction_shape_reitz_atomization_1978[i] == 'converging cone':
      K_L_reitz = K_L_wellrounded
   else:
      sys.exit('Invalid contraction shape for Reitz.')
   K_c_reitz_atomization_1978[i]   = K_c(P_atm_reitz_atomization_1978[i], P_v, rho_l, Ubar_0, K_L_reitz, friction_factor_smooth(Re_l0), L_0s_reitz_atomization_1978[i])
   Ma_g_reitz_atomization_1978[i]  = Ubar_0 / c_g
   
   # reitz_dependence_1979 p. 8L: > Dimensions taken from the steady state photographs have accuracies of 0.5$^\circ$ for jet divergence angles and 25% d_0 for jet intact lengths.
   theta_res_reitz_atomization_1978[i] = (2 * np.pi / 360) * 0.5
   x_is_res_reitz_atomization_1978[i]  = 0.25
   
   # TODO: Check that this is working correctly.
   if photo_filename_reitz_atomization_1978[i] is np.nan:
      x_low_reitz_atomization_1978[i] = np.nan
      photo_page_fig_reitz_atomization_1978.append(np.nan)
      spectrum_reitz_atomization_1978.append(np.nan)
      background_color_reitz_atomization_1978.append(np.nan)
      grid_reitz_atomization_1978.append(np.nan)
   else:
      x_low_reitz_atomization_1978[i] = 0
      photo_page_fig_reitz_atomization_1978.append('p. 71, fig. 6')
      spectrum_reitz_atomization_1978.append('visible')
      background_color_reitz_atomization_1978.append('light')
      grid_reitz_atomization_1978.append(False)
   
   if Re_l0 > Re_turb:
      regime_turb_reitz_atomization_1978.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_reitz_atomization_1978.append('transitional')
   else:
      regime_turb_reitz_atomization_1978.append('laminar')
   
   i = i + 1

I_0_reitz_atomization_1978 = I_fully_developed_array(Re_l0_reitz_atomization_1978, Re_turb)
Lambda_r_s_reitz_atomization_1978 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_reitz_atomization_1978, Re_turb), 'v')

df_reitz_atomization_1978 = pd.DataFrame({'We_l0': We_l0_reitz_atomization_1978})
df_reitz_atomization_1978['key']        = 'reitz_atomization_1978'
df_reitz_atomization_1978['alt key']    = 'reitz_dependence_1979'
df_reitz_atomization_1978['liquid']     = liquid_reitz_atomization_1978
df_reitz_atomization_1978['gas']        = gas_reitz_atomization_1978
df_reitz_atomization_1978['degas']      = np.nan # TODO: I vaguely recall seeing in his dissertation that he did something to degas the liquid.
df_reitz_atomization_1978['regulation'] = 'pressure regulator'

df_reitz_atomization_1978['bend']              = False
df_reitz_atomization_1978['flow straightener'] = np.nan # TODO
df_reitz_atomization_1978['contraction shape'] = contraction_shape_reitz_atomization_1978
df_reitz_atomization_1978['screen']            = np.nan # TODO
df_reitz_atomization_1978['c']                 = np.nan # TODO
df_reitz_atomization_1978['trip']              = np.nan # TODO
df_reitz_atomization_1978['L_r/d_0']           = 0
df_reitz_atomization_1978['eps_r/d_0']         = np.nan

df_reitz_atomization_1978['L_0/d_0']       = L_0s_reitz_atomization_1978
df_reitz_atomization_1978['roughness']     = 'smooth'
df_reitz_atomization_1978['f page fig']    = np.nan
df_reitz_atomization_1978['pipe material'] = 'stainless steel'
df_reitz_atomization_1978['est f']         = True
df_reitz_atomization_1978['check FD']      = False

df_reitz_atomization_1978['t/d_0']       = ts_reitz_atomization_1978
df_reitz_atomization_1978['L_tip/d_0']   = np.nan
df_reitz_atomization_1978['end checked'] = np.nan

df_reitz_atomization_1978['orientation']         = np.nan # Check p. 175-177 (pdf p. 199-201)
df_reitz_atomization_1978['vibration isolation'] = np.nan # TODO
df_reitz_atomization_1978['d_chamber/d_0']       = 500 # reitz_dependence_1979 p. 2L

df_reitz_atomization_1978['regime photo'] = np.nan
df_reitz_atomization_1978['regime L_b']   = np.nan
df_reitz_atomization_1978['regime turb']  = regime_turb_reitz_atomization_1978
df_reitz_atomization_1978['est turb']     = True

#df_reitz_atomization_1978['We_l0']      = We_l0_reitz_atomization_1978
df_reitz_atomization_1978['Re_l0']      = Re_l0_reitz_atomization_1978
df_reitz_atomization_1978['I_0']        = I_0_reitz_atomization_1978
df_reitz_atomization_1978['Fr_0']       = Fr_0_reitz_atomization_1978
df_reitz_atomization_1978['rho_s']      = rho_s_reitz_atomization_1978
df_reitz_atomization_1978['nu_s']       = nu_s_reitz_atomization_1978
df_reitz_atomization_1978['K_c']        = K_c_reitz_atomization_1978
df_reitz_atomization_1978['Ma_g']       = Ma_g_reitz_atomization_1978
df_reitz_atomization_1978['Lambda_r_s'] = Lambda_r_s_reitz_atomization_1978
df_reitz_atomization_1978['est rho_s']  = True
df_reitz_atomization_1978['est nu_s']   = True

df_reitz_atomization_1978['L_b/d_0']                      = np.nan
df_reitz_atomization_1978['L_b/d_0 stat error']           = np.nan
df_reitz_atomization_1978['L_b/d_0 resolution']           = np.nan
df_reitz_atomization_1978['L_b method']                   = np.nan
df_reitz_atomization_1978['L_b/d_0 page fig']             = np.nan
df_reitz_atomization_1978['L_b/d_0 transcription method'] = np.nan

df_reitz_atomization_1978['theta']                      = theta_reitz_atomization_1978
df_reitz_atomization_1978['theta stat error']           = np.nan # TODO
df_reitz_atomization_1978['theta resolution']           = theta_res_reitz_atomization_1978
df_reitz_atomization_1978['theta page fig']             = 'p. 161, tab. 4.2'
df_reitz_atomization_1978['theta transcription method'] = 'table'

df_reitz_atomization_1978['D_10/d_0']                   = np.nan
df_reitz_atomization_1978['D_10/d_0 stat error']        = np.nan
df_reitz_atomization_1978['D_10/d_0 resolution']        = np.nan
df_reitz_atomization_1978['D_30/d_0']                   = np.nan
df_reitz_atomization_1978['D_30/d_0 stat error']        = np.nan
df_reitz_atomization_1978['D_30/d_0 resolution']        = np.nan
df_reitz_atomization_1978['D_32/d_0']                   = np.nan
df_reitz_atomization_1978['D_32/d_0 stat error']        = np.nan
df_reitz_atomization_1978['D_32/d_0 resolution']        = np.nan
df_reitz_atomization_1978['droplet x/d_0']              = np.nan
df_reitz_atomization_1978['D/d_0 page fig']             = np.nan
df_reitz_atomization_1978['D/d_0 transcription method'] = np.nan
df_reitz_atomization_1978['D/d_0 measurement method']   = np.nan
df_reitz_atomization_1978['v_d_bar/vp']                        = np.nan
df_reitz_atomization_1978['v_d_bar/vp stat error']             = np.nan
df_reitz_atomization_1978['v_d_bar/vp resolution']             = np.nan
df_reitz_atomization_1978['v_d_bar/vp page fig']               = np.nan
df_reitz_atomization_1978['v_d_bar/vp transcription method']   = np.nan
df_reitz_atomization_1978['e_p']                        = np.nan
df_reitz_atomization_1978['e_p page fig']               = np.nan
df_reitz_atomization_1978['e_p transcription method']   = np.nan

df_reitz_atomization_1978['x_trans/d_0']                      = np.nan
df_reitz_atomization_1978['x_trans/d_0 page fig']             = np.nan
df_reitz_atomization_1978['x_trans/d_0 transcription method'] = np.nan

df_reitz_atomization_1978['x_i/d_0']                      = x_is_reitz_atomization_1978
df_reitz_atomization_1978['x_i/d_0 stat error']           = np.nan
df_reitz_atomization_1978['x_i/d_0 resolution']           = x_is_res_reitz_atomization_1978
df_reitz_atomization_1978['x_i/d_0 page fig']             = np.nan
df_reitz_atomization_1978['x_i/d_0 transcription method'] = np.nan

df_reitz_atomization_1978['x_e/d_0']                      = np.nan
df_reitz_atomization_1978['x_e/d_0 stat error']           = np.nan
df_reitz_atomization_1978['x_e/d_0 resolution']           = np.nan
df_reitz_atomization_1978['x_e/d_0 page fig']             = np.nan
df_reitz_atomization_1978['x_e/d_0 transcription method'] = np.nan

df_reitz_atomization_1978['photo filename']   = np.nan # TODO
df_reitz_atomization_1978['exposure time']    = np.nan
df_reitz_atomization_1978['flash']            = np.nan
df_reitz_atomization_1978['x_low']            = np.nan
df_reitz_atomization_1978['x_mid']            = np.nan
df_reitz_atomization_1978['x_high']           = np.nan
df_reitz_atomization_1978['photo page fig']   = np.nan
df_reitz_atomization_1978['spectrum']         = np.nan
df_reitz_atomization_1978['lighting']         = np.nan
df_reitz_atomization_1978['background color'] = np.nan
df_reitz_atomization_1978['grid']             = np.nan
df_reitz_atomization_1978['distance']         = np.nan
df_reitz_atomization_1978['f-stop']           = np.nan
df_reitz_atomization_1978['focal length']     = np.nan
df_reitz_atomization_1978['camera model']     = np.nan
df_reitz_atomization_1978['sensor']           = np.nan

df_reitz_atomization_1978 = df_reitz_atomization_1978[cols]
summary_table(df_reitz_atomization_1978)
df_jet_breakup = pd.concat([df_jet_breakup, df_reitz_atomization_1978])

#######################
# hoyt_pipe-exit_1980 #
#######################

# WON'T: Measure spray angle and add it. Use the average of the half angle on both sides to estimate.

# TODO: Estimate pipe thickness based on photo.
# 0.43 in = 217 px
# t = 17.5 px (average of both sides) = 0.035 in

# liquid = water
# stainless steel pipe
# d_0 = 0.43 in
# L_d/d_0 = 70
# Ubar_0 = 86 ft/sec
# Re_l0 = 3e5 ==> TODO: Find viscosity from this, then temperature, then other properties
# flash duration = 15e-6 s
# TODO: See other papers for information about photographic setup.
# alt key = hoyt_structure_1974 or taylor_water_1983

hoyt_pipeexit_1980_csv = pd.read_csv('../data/hoyt_pipe-exit_1980/hoyt_pipe.csv', sep=',', header=0)

citation_key_hoyt_pipeexit_1980   = hoyt_pipeexit_1980_csv['citation key']
Ubar_0_hoyt_pipeexit_1980         = hoyt_pipeexit_1980_csv['U_0 (ft/s)'] * 0.3048 # m/s
Re_l0_approx_hoyt_pipeexit_1980   = hoyt_pipeexit_1980_csv['Re_l0']
x_s_low_hoyt_pipeexit_1980        = hoyt_pipeexit_1980_csv['x/d_0 1']
photo_filename_hoyt_pipeexit_1980 = hoyt_pipeexit_1980_csv['photo filename 1']
photo_page_fig_hoyt_pipeexit_1980 = hoyt_pipeexit_1980_csv['photo page fig 1']
regime_photo_hoyt_pipeexit_1980 = hoyt_pipeexit_1980_csv['regime photo']

d_0 = 0.43 * 2.54e-2 # m

nu_l = 86. * 0.3048 * d_0 / 3.e5
T = temp_from_nu_water(nu_l)
assert(T > 0.)
assert(T < 30.)

sigma = sigma_water(T)
rho_l = rho_water(T)
rho_g = rho_ideal_gas(P_atm, T, MW_air)
nu_g  = mu_ideal_gas(T, mu_0_air, T_0_air, C_air) / rho_g

We_l0_hoyt_pipeexit_1980       = []
Re_l0_hoyt_pipeexit_1980       = []
Fr_0_hoyt_pipeexit_1980        = []
K_c_hoyt_pipeexit_1980         = []
Ma_g_hoyt_pipeexit_1980        = []
for Ubar_0, Re_l0 in zip(Ubar_0_hoyt_pipeexit_1980, Re_l0_approx_hoyt_pipeexit_1980):
   if np.isnan(Ubar_0):
      assert(not(np.isnan(Re_l0)))
      Ubar_0 = Re_l0 * nu_l / d_0
   elif np.isnan(Re_l0):
      assert(not(np.isnan(Ubar_0)))
      Re_l0 = Ubar_0 * d_0 / nu_l
   else:
      print('Invalid hoyt_pipe-exit_1980 data.')
      sys.exit()
   
   We_l0 = rho_l * Ubar_0**2. * d_0 / sigma
   
   We_l0_hoyt_pipeexit_1980.append(We_l0)
   Re_l0_hoyt_pipeexit_1980.append(Re_l0)
   Fr_0_hoyt_pipeexit_1980.append(Ubar_0**2 / (g * d_0))
   #K_c_hoyt_pipeexit_1980.append((P_atm - P_v) / (0.5 * rho_l * Ubar_0**2))
   K_c_hoyt_pipeexit_1980.append(K_c(P_atm, P_v, rho_l, Ubar_0, K_L_guess, friction_factor_smooth(Re_l0), 70.))
   Ma_g_hoyt_pipeexit_1980.append(Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_air))

I_0_hoyt_pipeexit_1980        = I_fully_developed_array(Re_l0_hoyt_pipeexit_1980, Re_turb)
Lambda_r_s_hoyt_pipeexit_1980 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_hoyt_pipeexit_1980, Re_turb), 'v')

df_hoyt_pipeexit_1980 = pd.DataFrame({'We_l0': We_l0_hoyt_pipeexit_1980})
df_hoyt_pipeexit_1980['key']        = citation_key_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['alt key']    = np.nan
df_hoyt_pipeexit_1980['liquid']     = 'water'
df_hoyt_pipeexit_1980['gas']        = 'air'
df_hoyt_pipeexit_1980['degas']      = np.nan
df_hoyt_pipeexit_1980['regulation'] = np.nan

df_hoyt_pipeexit_1980['bend']              = np.nan
df_hoyt_pipeexit_1980['flow straightener'] = np.nan
df_hoyt_pipeexit_1980['contraction shape'] = np.nan
df_hoyt_pipeexit_1980['screen']            = np.nan
df_hoyt_pipeexit_1980['c']                 = np.nan
df_hoyt_pipeexit_1980['trip']              = np.nan
df_hoyt_pipeexit_1980['L_r/d_0']           = 0.
df_hoyt_pipeexit_1980['eps_r/d_0']         = np.nan

df_hoyt_pipeexit_1980['L_0/d_0']       = 70.
df_hoyt_pipeexit_1980['roughness']     = 'smooth'
df_hoyt_pipeexit_1980['f page fig']    = np.nan
df_hoyt_pipeexit_1980['pipe material'] = 'stainless steel'
df_hoyt_pipeexit_1980['est f']         = True
df_hoyt_pipeexit_1980['check FD']      = False

df_hoyt_pipeexit_1980['t/d_0']       = 0.035 / 0.43
df_hoyt_pipeexit_1980['L_tip/d_0']   = np.nan
df_hoyt_pipeexit_1980['end checked'] = np.nan

df_hoyt_pipeexit_1980['orientation']         = np.nan
df_hoyt_pipeexit_1980['vibration isolation'] = np.nan
df_hoyt_pipeexit_1980['d_chamber/d_0']       = np.nan

df_hoyt_pipeexit_1980['regime photo'] = regime_photo_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['regime L_b']   = np.nan
df_hoyt_pipeexit_1980['regime turb']  = 'turbulent'
df_hoyt_pipeexit_1980['est turb']     = False

#df_hoyt_pipeexit_1980['We_l0']      = We_l0_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['Re_l0']      = Re_l0_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['I_0']        = I_0_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['Fr_0']       = Fr_0_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['rho_s']      = rho_l / rho_g
df_hoyt_pipeexit_1980['nu_s']       = nu_l / nu_g
df_hoyt_pipeexit_1980['K_c']        = K_c_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['Ma_g']       = Ma_g_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['Lambda_r_s'] = Lambda_r_s_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['est rho_s']  = True
df_hoyt_pipeexit_1980['est nu_s']   = True

df_hoyt_pipeexit_1980['L_b/d_0']                      = np.nan
df_hoyt_pipeexit_1980['L_b/d_0 stat error']           = np.nan
df_hoyt_pipeexit_1980['L_b/d_0 resolution']           = np.nan
df_hoyt_pipeexit_1980['L_b method']                   = np.nan
df_hoyt_pipeexit_1980['L_b/d_0 page fig']             = np.nan
df_hoyt_pipeexit_1980['L_b/d_0 transcription method'] = np.nan

df_hoyt_pipeexit_1980['theta']                      = np.nan
df_hoyt_pipeexit_1980['theta stat error']           = np.nan
df_hoyt_pipeexit_1980['theta resolution']           = np.nan
df_hoyt_pipeexit_1980['theta page fig']             = np.nan
df_hoyt_pipeexit_1980['theta transcription method'] = np.nan

df_hoyt_pipeexit_1980['D_10/d_0']                   = np.nan
df_hoyt_pipeexit_1980['D_10/d_0 stat error']        = np.nan
df_hoyt_pipeexit_1980['D_10/d_0 resolution']        = np.nan
df_hoyt_pipeexit_1980['D_30/d_0']                   = np.nan
df_hoyt_pipeexit_1980['D_30/d_0 stat error']        = np.nan
df_hoyt_pipeexit_1980['D_30/d_0 resolution']        = np.nan
df_hoyt_pipeexit_1980['D_32/d_0']                   = np.nan
df_hoyt_pipeexit_1980['D_32/d_0 stat error']        = np.nan
df_hoyt_pipeexit_1980['D_32/d_0 resolution']        = np.nan
df_hoyt_pipeexit_1980['droplet x/d_0']              = np.nan
df_hoyt_pipeexit_1980['D/d_0 page fig']             = np.nan
df_hoyt_pipeexit_1980['D/d_0 transcription method'] = np.nan
df_hoyt_pipeexit_1980['D/d_0 measurement method']   = np.nan
df_hoyt_pipeexit_1980['v_d_bar/vp']                        = np.nan
df_hoyt_pipeexit_1980['v_d_bar/vp stat error']             = np.nan
df_hoyt_pipeexit_1980['v_d_bar/vp resolution']             = np.nan
df_hoyt_pipeexit_1980['v_d_bar/vp page fig']               = np.nan
df_hoyt_pipeexit_1980['v_d_bar/vp transcription method']   = np.nan
df_hoyt_pipeexit_1980['e_p']                        = np.nan
df_hoyt_pipeexit_1980['e_p page fig']               = np.nan
df_hoyt_pipeexit_1980['e_p transcription method']   = np.nan

df_hoyt_pipeexit_1980['x_trans/d_0']                      = np.nan
df_hoyt_pipeexit_1980['x_trans/d_0 page fig']             = np.nan
df_hoyt_pipeexit_1980['x_trans/d_0 transcription method'] = np.nan

df_hoyt_pipeexit_1980['x_i/d_0']                      = np.nan
df_hoyt_pipeexit_1980['x_i/d_0 stat error']           = np.nan
df_hoyt_pipeexit_1980['x_i/d_0 resolution']           = np.nan
df_hoyt_pipeexit_1980['x_i/d_0 page fig']             = np.nan
df_hoyt_pipeexit_1980['x_i/d_0 transcription method'] = np.nan

df_hoyt_pipeexit_1980['x_e/d_0']                      = np.nan
df_hoyt_pipeexit_1980['x_e/d_0 stat error']           = np.nan
df_hoyt_pipeexit_1980['x_e/d_0 resolution']           = np.nan
df_hoyt_pipeexit_1980['x_e/d_0 page fig']             = np.nan
df_hoyt_pipeexit_1980['x_e/d_0 transcription method'] = np.nan

df_hoyt_pipeexit_1980['photo filename']   = photo_filename_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['exposure time']    = 15e-6 # s
df_hoyt_pipeexit_1980['flash']            = True
df_hoyt_pipeexit_1980['x_low']            = x_s_low_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['x_mid']            = np.nan
df_hoyt_pipeexit_1980['x_high']           = np.nan
df_hoyt_pipeexit_1980['photo page fig']   = photo_page_fig_hoyt_pipeexit_1980
df_hoyt_pipeexit_1980['spectrum']         = np.nan
df_hoyt_pipeexit_1980['lighting']         = np.nan
df_hoyt_pipeexit_1980['background color'] = np.nan
df_hoyt_pipeexit_1980['grid']             = np.nan
df_hoyt_pipeexit_1980['distance']         = np.nan
df_hoyt_pipeexit_1980['f-stop']           = np.nan
df_hoyt_pipeexit_1980['focal length']     = np.nan
df_hoyt_pipeexit_1980['camera model']     = np.nan
df_hoyt_pipeexit_1980['sensor']           = np.nan

df_hoyt_pipeexit_1980 = df_hoyt_pipeexit_1980[cols]
summary_table(df_hoyt_pipeexit_1980)
df_jet_breakup = pd.concat([df_jet_breakup, df_hoyt_pipeexit_1980])

############################
# iciek_hydrodynamics_1982 #
############################

# TODO: Add.

##################################################
# kim_investigation_1983 / kim_condensation_1989 #
##################################################

# Issues with this data:
# 1. kim_condensation_1989 p. 1069L suggests the diameter of the rough tubes is 5.5 mm, but kim_investigation_1983 p. 68 fig. 5.8's caption says 4 mm. The caption in the 1989 paper does not give the rough tube diameter. The 1989 paper appears to have some corrections, e.g., the roughness is missing (i.e., not typeset and just left blank) in the figure in the 1983 dissertation. Also, the caption says "$L = 10 \pm 5$ cm", which is not true for 5.8b, so that's another reason to trust the 1989 paper over the 1983 dissertation. Table 5.5 on p. 98 of kim_condensation_1989 gives the diameters as 5.4 mm and 5.5 mm for the two different rough tubes, along with the L_0/d_0. This table gives a different value of the roughness than the caption, though, so it's unclear which tube is correct or even if the photo corresponds to one of these tubes. I'll assume the friction factor of the tube with the closer roughness is correct, so the friction factor is 0.054. I also assume that the photo corresponds to the fully rough regime.
# TODO: Check if the Reynolds number for the rough pipe photo is high enough to be in the fully rough regime for both roughnesses given.
# 2. kim_condensation_1989 p. 1068R suggests that "nozzle lengths of greater than 30 diameters are used", but kim_investigation_1983 p. 98 tab. 5 suggests the rough tubes have the aspect ratio L_0/d_0 = 19. I'll assume the other tubes have lengths of at least 30 diameters.
# 3. It's unclear whether the scales given for figures b and c start at the start of the line or the nozzle. I assume they start at the nozzle.

# Is L the length from the nozzle or the length of the pipe? The 1989 paper's nomenclature suggests length from the nozzle.

# nozzle pipe length:
# kim_condensation_1989 p. 1068R
# kim_investigation_1983 p. 98 tab. 5

# liquid property data:
# kim_investigation_1983 pp. 57-58 (pdf p. 75-76).
# Just reproduces data from elsewhere.

# atmospheric pressure:
# kim_investigation_1983 p. 55 mentions "barometric pressure and ambient temperature" as variables which were recorded.
# kim_condensation_1989 p. 1068R: > The pressure in the test chamber is measured by a mercury manometer. [...] The reference pressure is maintained at less than 50 $\mu$m Hg absolute"
P_atm_kim_investigation_1983 = 6772.777 # Pa

# atmospheric temperature:
# kim_investigation_1983 p. xvi: > saturation temperatures 10 - 30$^\deg$C
# kim_investigation_1983 p. 8: > vapor saturation temperature of 20-40$^\circ$C
# Get atmospheric temperature from saturation temperature at given pressure? Or assume the saturation temperature equals the gas temperature? Let's go with 30 C.
T_kim_investigation_1983 = 30 # C

# precision:
# kim_condensation_1989 p. 1069L: > A Fluke Data Logger with a built-in cold junction compensator and least count of 0.1 $^\circ$ F (0.06 $^\circ$ C) is used to measure the thermocouple outlets.

# kim_investigation_1983 p. 8: > The initial turbulence level in the jet is controlled by using long tube nozzles so as to be able to use established relations for fully developed turbulent flow for the estimation of initial turbulence level. An attempt was made to vary the initial turbulence level by using rough nozzles.

# kim_investigation_1983 pp. 174-200 (pdf pp. 192-)

kim_investigation_1983_csv = pd.read_csv('../data/kim_investigation_1983/kim_investigation_1983_photos.csv', sep=',', header=0)

liquid_kim_investigation_1983         = kim_investigation_1983_csv['fluid']
pipe_material_kim_investigation_1983  = kim_investigation_1983_csv['pipe material']
roughness_kim_investigation_1983      = kim_investigation_1983_csv['roughness']
page_fig_kim_investigation_1983       = kim_investigation_1983_csv['page fig']
Ubar_0_kim_investigation_1983         = kim_investigation_1983_csv['U_0']
x_low_kim_investigation_1983          = kim_investigation_1983_csv['x_low']
x_high_kim_investigation_1983         = kim_investigation_1983_csv['x_high']
f_kim_investigation_1983              = kim_investigation_1983_csv['f']
f_page_fig_kim_investigation_1983     = kim_investigation_1983_csv['f page fig']
T_l_kim_investigation_1983            = kim_investigation_1983_csv['T_l (C)']
d_0_kim_investigation_1983            = kim_investigation_1983_csv['d_0 (m)']
L_0_kim_investigation_1983            = kim_investigation_1983_csv['L_0 (m)']
x_i_kim_investigation_1983            = kim_investigation_1983_csv['x_i (m)']
photo_filename_kim_investigation_1983 = kim_investigation_1983_csv['photo filename']

x_is_kim_investigation_1983 = x_i_kim_investigation_1983 / d_0_kim_investigation_1983

rho_g = rho_ideal_gas(P_atm_kim_investigation_1983, T_kim_investigation_1983, MW_air)
nu_g  = mu_ideal_gas(T_kim_investigation_1983, mu_0_air, T_0_air, C_air) / rho_g

i = 0

We_l0_kim_investigation_1983      = np.zeros(len(Ubar_0_kim_investigation_1983))
Re_l0_kim_investigation_1983      = np.zeros(len(Ubar_0_kim_investigation_1983))
I_0_kim_investigation_1983        = np.zeros(len(Ubar_0_kim_investigation_1983))
Fr_0_kim_investigation_1983       = np.zeros(len(Ubar_0_kim_investigation_1983))
rho_s_kim_investigation_1983      = np.zeros(len(Ubar_0_kim_investigation_1983))
nu_s_kim_investigation_1983       = np.zeros(len(Ubar_0_kim_investigation_1983))
K_c_kim_investigation_1983        = np.zeros(len(Ubar_0_kim_investigation_1983))
Ma_g_kim_investigation_1983       = np.zeros(len(Ubar_0_kim_investigation_1983))
Lambda_r_s_kim_investigation_1983 = np.zeros(len(Ubar_0_kim_investigation_1983))
L_0s_kim_investigation_1983       = np.zeros(len(Ubar_0_kim_investigation_1983))
f_page_fig_kim_investigation_1983  = []
regime_turb_kim_investigation_1983 = []
est_turb_kim_investigation_1983    = []
est_f_kim_investigation_1983       = []
for Ubar_0 in Ubar_0_kim_investigation_1983:
   d_0       = d_0_kim_investigation_1983[i]
   liquid    = liquid_kim_investigation_1983[i]
   T_l       = T_l_kim_investigation_1983[i]
   roughness = roughness_kim_investigation_1983[i]
   
   if liquid == 'water':
      rho_l = rho_water(T_l)
      sigma = sigma_water(T_l)
      nu_l  = nu_water(T_l)
      P_v   = P_v_water(T_l)
   elif liquid == 'ethanol':
      # https://en.wikipedia.org/wiki/Ethanol_(data_page)
      rho_l = 0.78860e3 # kg/m^3
      sigma = 22.39e-3 # N/m
      nu_l  = 1.2e-3 / rho_l
      P_v   = 5.95e3 # Pa
   else:
      raise ValueError('Invalid liquid for kim_investigation_1983:', liquid)
   
   We_l0_kim_investigation_1983[i] = rho_l * Ubar_0**2 * d_0 / sigma
   Re_l0_kim_investigation_1983[i] = rho_l * Ubar_0 * d_0 / mu_l
   Fr_0_kim_investigation_1983[i]  = Ubar_0**2 / (g * d_0)
   rho_s_kim_investigation_1983[i] = rho_l / rho_g
   nu_s_kim_investigation_1983[i]  = nu_l / nu_g
   Ma_g_kim_investigation_1983[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_kim_investigation_1983), T_kim_investigation_1983 + T_zero, MW_N2)
   
   if L_0_kim_investigation_1983[i] > 0:
      L_0s_kim_investigation_1983[i] = L_0_kim_investigation_1983[i] / d_0
   else:
      L_0s_kim_investigation_1983[i] = 30
   
   Re_l0 = Re_l0_kim_investigation_1983[i]
   f = f_kim_investigation_1983[i]
   
   if f > 0:
      if (f / f_laminar > f_ratio_cutoff):
         regime_turb_kim_investigation_1983.append('turbulent')
      else:
         regime_turb_kim_investigation_1983.append('laminar')
      
      est_turb_kim_investigation_1983.append(False)
      est_f_kim_investigation_1983.append(False)
   else:
      if Re_l0 > Re_turb:
         regime_turb_kim_investigation_1983.append('turbulent')
      elif Re_l0 > Re_trans:
         regime_turb_kim_investigation_1983.append('transitional')
      else:
         regime_turb_kim_investigation_1983.append('laminar')
      
      est_turb_kim_investigation_1983.append(True)
      est_f_kim_investigation_1983.append(True)
   
   # friction factor:
   # kim_investigation_1983 p. 98
   # kim_condensation_1989 p. 1069L: > two roughenened nozzles were prepared from 5.5 mm I.D. tube to give values of $k_s/d$ of 0.026 and 0.035 (fully rough friction factors of 0.054 and 0.060, respectively)
   # TODO: Fix this so that it assigns the right I_0. Was this fixed?
   if f_kim_investigation_1983[i] > 0:
      I_0_kim_investigation_1983[i] = I_fully_developed(f_kim_investigation_1983[i])
      f_page_fig_kim_investigation_1983.append('p. 98') # TODO
      Lambda_r_s_kim_investigation_1983[i] = Lambda_s(roughness, f_kim_investigation_1983[i], 'v')
   else:
      I_0_kim_investigation_1983[i] = I_fully_developed(friction_factor_smooth(Re_l0_kim_investigation_1983[i]))
      f_page_fig_kim_investigation_1983.append(np.nan)
      Lambda_r_s_kim_investigation_1983[i] = Lambda_s(roughness, friction_factor_smooth(Re_l0_kim_investigation_1983[i]), 'v')
   
   #K_c_kim_investigation_1983[i] = (P_atm_kim_investigation_1983 - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_kim_investigation_1983[i] = K_c(P_atm_kim_investigation_1983, P_v, rho_l, Ubar_0, K_L_guess, f, L_0s_kim_investigation_1983[i])
   
   i = i + 1

df_kim_investigation_1983 = pd.DataFrame({'We_l0': We_l0_kim_investigation_1983})
df_kim_investigation_1983['key']        = 'kim_investigation_1983'
df_kim_investigation_1983['alt key']    = 'kim_condensation_1989'
df_kim_investigation_1983['liquid']     = liquid_kim_investigation_1983
df_kim_investigation_1983['gas']        = 'air'
df_kim_investigation_1983['degas']      = True # kim_condensation_1989 p. 1070L: > 1 Deaeration and system cooldown [...]
df_kim_investigation_1983['regulation'] = True # kim_condensation_1989 p. 1070L

df_kim_investigation_1983['bend']              = np.nan
df_kim_investigation_1983['flow straightener'] = np.nan # p. 48 mentions a "calming section"
df_kim_investigation_1983['contraction shape'] = np.nan
df_kim_investigation_1983['screen']            = np.nan
df_kim_investigation_1983['c']                 = np.nan
df_kim_investigation_1983['trip']              = True # kim_condensation_1989 p. 1068R: > we chose to use a long nozzle with an inlet turbulence trip
df_kim_investigation_1983['L_r/d_0']           = 0
df_kim_investigation_1983['eps_r/d_0']         = np.nan

df_kim_investigation_1983['L_0/d_0']       = L_0s_kim_investigation_1983
df_kim_investigation_1983['roughness']     = roughness_kim_investigation_1983
df_kim_investigation_1983['f page fig']    = np.nan # TODO
df_kim_investigation_1983['pipe material'] = 'glass'
df_kim_investigation_1983['est f']         = est_f_kim_investigation_1983
df_kim_investigation_1983['check FD']      = False

df_kim_investigation_1983['t/d_0']       = np.nan
df_kim_investigation_1983['L_tip/d_0']   = np.nan
df_kim_investigation_1983['end checked'] = np.nan

df_kim_investigation_1983['orientation']         = 'down' # presumably given that the vapor and jet enter at the "top"
df_kim_investigation_1983['vibration isolation'] = np.nan
df_kim_investigation_1983['d_chamber/d_0']       = 12e-3 / d_0_kim_investigation_1983 # kim_condensation_1989 p. 1068R

df_kim_investigation_1983['regime photo'] = 'second wind-induced'
df_kim_investigation_1983['regime L_b']   = np.nan
df_kim_investigation_1983['regime turb']  = regime_turb_kim_investigation_1983
df_kim_investigation_1983['est turb']     = est_turb_kim_investigation_1983

#df_kim_investigation_1983['We_l0']      = We_l0_kim_investigation_1983
df_kim_investigation_1983['Re_l0']      = Re_l0_kim_investigation_1983
df_kim_investigation_1983['I_0']        = I_0_kim_investigation_1983
df_kim_investigation_1983['Fr_0']       = Fr_0_kim_investigation_1983
df_kim_investigation_1983['rho_s']      = rho_s_kim_investigation_1983
df_kim_investigation_1983['nu_s']       = nu_s_kim_investigation_1983
df_kim_investigation_1983['K_c']        = K_c_kim_investigation_1983
df_kim_investigation_1983['Ma_g']       = Ma_g_kim_investigation_1983
df_kim_investigation_1983['Lambda_r_s'] = Lambda_r_s_kim_investigation_1983
df_kim_investigation_1983['est rho_s']  = True
df_kim_investigation_1983['est nu_s']   = True

df_kim_investigation_1983['L_b/d_0']                      = np.nan
df_kim_investigation_1983['L_b/d_0 stat error']           = np.nan
df_kim_investigation_1983['L_b/d_0 resolution']           = np.nan
df_kim_investigation_1983['L_b method']                   = np.nan
df_kim_investigation_1983['L_b/d_0 page fig']             = np.nan
df_kim_investigation_1983['L_b/d_0 transcription method'] = np.nan

df_kim_investigation_1983['theta']                      = np.nan
df_kim_investigation_1983['theta stat error']           = np.nan
df_kim_investigation_1983['theta resolution']           = np.nan
df_kim_investigation_1983['theta page fig']             = np.nan
df_kim_investigation_1983['theta transcription method'] = np.nan

df_kim_investigation_1983['D_10/d_0']                   = np.nan
df_kim_investigation_1983['D_10/d_0 stat error']        = np.nan
df_kim_investigation_1983['D_10/d_0 resolution']        = np.nan
df_kim_investigation_1983['D_30/d_0']                   = np.nan
df_kim_investigation_1983['D_30/d_0 stat error']        = np.nan
df_kim_investigation_1983['D_30/d_0 resolution']        = np.nan
df_kim_investigation_1983['D_32/d_0']                   = np.nan
df_kim_investigation_1983['D_32/d_0 stat error']        = np.nan
df_kim_investigation_1983['D_32/d_0 resolution']        = np.nan
df_kim_investigation_1983['droplet x/d_0']              = np.nan
df_kim_investigation_1983['D/d_0 page fig']             = np.nan
df_kim_investigation_1983['D/d_0 transcription method'] = np.nan
df_kim_investigation_1983['D/d_0 measurement method']   = np.nan
df_kim_investigation_1983['v_d_bar/vp']                        = np.nan
df_kim_investigation_1983['v_d_bar/vp stat error']             = np.nan
df_kim_investigation_1983['v_d_bar/vp resolution']             = np.nan
df_kim_investigation_1983['v_d_bar/vp page fig']               = np.nan
df_kim_investigation_1983['v_d_bar/vp transcription method']   = np.nan
df_kim_investigation_1983['e_p']                        = np.nan
df_kim_investigation_1983['e_p page fig']               = np.nan
df_kim_investigation_1983['e_p transcription method']   = np.nan

df_kim_investigation_1983['x_trans/d_0']                      = np.nan
df_kim_investigation_1983['x_trans/d_0 page fig']             = np.nan
df_kim_investigation_1983['x_trans/d_0 transcription method'] = np.nan

df_kim_investigation_1983['x_i/d_0']                      = np.nan
df_kim_investigation_1983['x_i/d_0 stat error']           = np.nan
df_kim_investigation_1983['x_i/d_0 resolution']           = np.nan
df_kim_investigation_1983['x_i/d_0 page fig']             = np.nan
df_kim_investigation_1983['x_i/d_0 transcription method'] = np.nan

df_kim_investigation_1983['x_e/d_0']                      = np.nan
df_kim_investigation_1983['x_e/d_0 stat error']           = np.nan
df_kim_investigation_1983['x_e/d_0 resolution']           = np.nan
df_kim_investigation_1983['x_e/d_0 page fig']             = np.nan
df_kim_investigation_1983['x_e/d_0 transcription method'] = np.nan

# kim_investigation_1983 p. 65 (pdf p. 83)
df_kim_investigation_1983['photo filename']   = photo_filename_kim_investigation_1983
df_kim_investigation_1983['exposure time']    = 1/3000
df_kim_investigation_1983['flash']            = True
df_kim_investigation_1983['x_low']            = x_low_kim_investigation_1983
df_kim_investigation_1983['x_mid']            = np.nan
df_kim_investigation_1983['x_high']           = x_high_kim_investigation_1983
df_kim_investigation_1983['photo page fig']   = page_fig_kim_investigation_1983
df_kim_investigation_1983['spectrum']         = 'visible'
df_kim_investigation_1983['lighting']         = np.nan # Unknown
df_kim_investigation_1983['background color'] = 'light'
df_kim_investigation_1983['grid']             = True
df_kim_investigation_1983['distance']         = np.nan
df_kim_investigation_1983['f-stop']           = np.nan
df_kim_investigation_1983['focal length']     = np.nan
df_kim_investigation_1983['camera model']     = 'Redlake Hycam 41-004'
df_kim_investigation_1983['sensor']           = '16 mm Kodak 4-X reversal film'

df_kim_investigation_1983 = df_kim_investigation_1983[cols]
summary_table(df_kim_investigation_1983)
df_jet_breakup = pd.concat([df_jet_breakup, df_kim_investigation_1983])

###########################################
# wu_atomizing_1983 / wu_measurement_1983 #
###########################################

# spray angle: wu_measurement_1983 fig. 9
# data table: wu_atomizing_1983 p. 116, tab. 2.5

# Issues with this data:
# 1. Number of data points is assumed to be in the "Data" column in table 2.5 (p. 116). This does not appear to be clarified in the text. The worst case scenario is 3. See p. 33: "At least three photographs were taken for each condition and the arithmetic average of the spray angles was calculated." This seems consistent with the "Data" column in table 2.5. Page A59 has 4 photos from this test series, so it appears there were 4 photos at least.
# 2. The rho_g/rho_l column appears to have a typo. wu_measurement_1983 suggests (caption for fig. 9) that it should read rho_g/rho_l \times 10^{-3}.

# Test with high pressure apparatus.

wu_atomizing_1983_csv = pd.read_csv('../data/wu_atomizing_1983/wu_atomizing_1983.csv', sep=',', header=0)

n_wu_atomizing_1983      = wu_atomizing_1983_csv['data'] # See issue #1.
nozzle_wu_atomizing_1983 = wu_atomizing_1983_csv['nozzle']
P_g_wu_atomizing_1983    = wu_atomizing_1983_csv['P_g (MPa)'] * 1e6 # Pa
rho_g_wu_atomizing_1983  = wu_atomizing_1983_csv['rho_g (kg/m^3)']
liquid_wu_atomizing_1983 = wu_atomizing_1983_csv['liquid']
dP_wu_atomizing_1983     = wu_atomizing_1983_csv['dP (MPa)'] * 1e6 # Pa
Ubar_0_wu_atomizing_1983 = wu_atomizing_1983_csv['Vbar_inj (m/s)']
Re_l0_wu_atomizing_1983  = wu_atomizing_1983_csv['Re_l0']
We_l0_wu_atomizing_1983  = wu_atomizing_1983_csv['We_l0']
theta_wu_atomizing_1983  = wu_atomizing_1983_csv['theta (deg.)'] * (np.pi / 180) # radians
rho_s_wu_atomizing_1983  = 1 / wu_atomizing_1983_csv['rho_g/rho_l']

# p. 114, tab. 2.3
mu_g = 1.70e-5 # kg/(m*s

i = 0
d_0_wu_atomizing_1983     = np.zeros(len(Ubar_0_wu_atomizing_1983))
L_0s_wu_atomizing_1983    = np.zeros(len(Ubar_0_wu_atomizing_1983))
Fr_0_wu_atomizing_1983    = np.zeros(len(Ubar_0_wu_atomizing_1983))
nu_s_wu_atomizing_1983    = np.zeros(len(Ubar_0_wu_atomizing_1983))
K_c_wu_atomizing_1983     = np.zeros(len(Ubar_0_wu_atomizing_1983))
Ma_g_wu_atomizing_1983    = np.zeros(len(Ubar_0_wu_atomizing_1983))
e_theta_wu_atomizing_1983 = np.zeros(len(Ubar_0_wu_atomizing_1983))
regime_turb_wu_atomizing_1983 = []
for nozzle in nozzle_wu_atomizing_1983:
   liquid = liquid_wu_atomizing_1983[i]
   Ubar_0 = Ubar_0_wu_atomizing_1983[i]
   Re_l0  = Re_l0_wu_atomizing_1983[i]
   We_l0  = We_l0_wu_atomizing_1983[i]
   rho_g  = rho_g_wu_atomizing_1983[i]
   P_g    = P_g_wu_atomizing_1983[i]
   
   if nozzle == '330-50':
      d_0_wu_atomizing_1983[i]  = 330e-6 # m
      d_0                       = d_0_wu_atomizing_1983[i]
      L_0s_wu_atomizing_1983[i] = 50
   else:
      raise ValueError('Invalid nozzle for wu_atomizing_1983:', nozzle)
   
   if liquid == 'n-hexane':
      # p. 115, tab. 2.4
      rho_l = 665     # kg/m^3
      mu_l  = 0.00032 # kg/(m*s)
      sigma = 0.01843 # N/m
      P_v   = 16.5e6  # Pa
   else:
      raise ValueError('Invalid liquid for wu_atomizing_1983:', liquid)
   
   nu_g = mu_g / rho_g
   rho_s_est = rho_l / rho_g
   
   # Check that the density ratio is reasonable due to the typo in there.
   #print rho_s_est, rho_s_wu_atomizing_1983[i], abs(rho_s_wu_atomizing_1983[i] - rho_s_est)
   assert(abs(rho_s_wu_atomizing_1983[i] - rho_s_est) < 1)
   
   Fr_0_wu_atomizing_1983[i]  = Ubar_0**2 / (g * d_0)
   nu_s_wu_atomizing_1983[i]  = nu_l / nu_g
   #K_c_wu_atomizing_1983[i]   = (P_g - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_wu_atomizing_1983[i]   = K_c(P_g, P_v, rho_l, Ubar_0, K_L_wellrounded, friction_factor_smooth(Re_l0), L_0s_wu_atomizing_1983[i])
   Ma_g_wu_atomizing_1983[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_N2)
   
   if Re_l0 > Re_turb:
      regime_turb_wu_atomizing_1983.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_wu_atomizing_1983.append('transitional')
   else:
      regime_turb_wu_atomizing_1983.append('laminar')
   
   # Worst case scenario. wu_atomizing_1983: > The standard deviation of each condition is normally less than 3% but in some cases it was as large as 7%.
   e_theta_wu_atomizing_1983[i] = max(z * 0.07 * theta_wu_atomizing_1983[i] / sqrt(n_wu_atomizing_1983[i]), 0.5 * (np.pi / 180)) # radians

I_0_wu_atomizing_1983 = I_fully_developed_array(Re_l0_wu_atomizing_1983, Re_turb)
Lambda_r_s_wu_atomizing_1983 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_wu_atomizing_1983, Re_turb), 'v')

df_wu_atomizing_1983 = pd.DataFrame({'We_l0': We_l0_wu_atomizing_1983})
df_wu_atomizing_1983['key']        = 'wu_atomizing_1983'
df_wu_atomizing_1983['alt key']    = 'wu_measurement_1983'
df_wu_atomizing_1983['liquid']     = 'n-hexane'
df_wu_atomizing_1983['gas']        = 'N2' # wu_atomizing_1983 p. 114, tab. 2.3
df_wu_atomizing_1983['degas']      = np.nan
df_wu_atomizing_1983['regulation'] = np.nan

df_wu_atomizing_1983['bend']              = False
df_wu_atomizing_1983['flow straightener'] = False
df_wu_atomizing_1983['contraction shape'] = 'smooth'  # wu_measurement_1983 p. 409R fig. 4; see also wu_atomizing_1983 p. 132 fig 2.4
df_wu_atomizing_1983['screen']            = False
df_wu_atomizing_1983['c']                 = 9.525e-3 / d_0_wu_atomizing_1983 # wu_measurement_1983 p. 409R fig. 4; see also wu_atomizing_1983 p. 132 fig 2.4
df_wu_atomizing_1983['trip']              = False
df_wu_atomizing_1983['L_r/d_0']           = 0
df_wu_atomizing_1983['eps_r/d_0']         = np.nan

df_wu_atomizing_1983['L_0/d_0']       = L_0s_wu_atomizing_1983
df_wu_atomizing_1983['roughness']     = 'smooth'
df_wu_atomizing_1983['f page fig']    = np.nan
df_wu_atomizing_1983['pipe material'] = '304 stainless steel' # wu_atomizing_1983 p. 23, wu_atomizing_1983 p. 113, tab. 2.2: > made from stainless steel hypodermic tubing
df_wu_atomizing_1983['est f']         = True
df_wu_atomizing_1983['check FD']      = np.nan

df_wu_atomizing_1983['t/d_0']       = (635e-6 - d_0_wu_atomizing_1983) / (2 * d_0_wu_atomizing_1983) # wu_atomizing_1983 p. 23
df_wu_atomizing_1983['L_tip/d_0']   = np.nan
df_wu_atomizing_1983['end checked'] = True # wu_atomizing_1983 p. 23: > Both ends of the tubing were then carefully cleaned under microscope

df_wu_atomizing_1983['orientation']         = np.nan # Unclear
df_wu_atomizing_1983['vibration isolation'] = np.nan
df_wu_atomizing_1983['d_chamber/d_0']       = 19e-2 / d_0_wu_atomizing_1983 # wu_atomizing_1983 p. 18

df_wu_atomizing_1983['regime photo'] = 'atomization'
df_wu_atomizing_1983['regime L_b']   = np.nan
df_wu_atomizing_1983['regime turb']  = regime_turb_wu_atomizing_1983
df_wu_atomizing_1983['est turb']     = True

#df_wu_atomizing_1983['We_l0']      = We_l0_wu_atomizing_1983
df_wu_atomizing_1983['Re_l0']      = Re_l0_wu_atomizing_1983
df_wu_atomizing_1983['I_0']        = I_0_wu_atomizing_1983
df_wu_atomizing_1983['Fr_0']       = Fr_0_wu_atomizing_1983
df_wu_atomizing_1983['rho_s']      = rho_s_wu_atomizing_1983
df_wu_atomizing_1983['nu_s']       = nu_s_wu_atomizing_1983
df_wu_atomizing_1983['K_c']        = K_c_wu_atomizing_1983
df_wu_atomizing_1983['Ma_g']       = Ma_g_wu_atomizing_1983
df_wu_atomizing_1983['Lambda_r_s'] = Lambda_r_s_wu_atomizing_1983
df_wu_atomizing_1983['est rho_s']  = False
df_wu_atomizing_1983['est nu_s']   = False # All data given in wu_atomizing_1983 p. 114, tab. 2.3.

df_wu_atomizing_1983['L_b/d_0']                      = np.nan
df_wu_atomizing_1983['L_b/d_0 stat error']           = np.nan
df_wu_atomizing_1983['L_b/d_0 resolution']           = np.nan
df_wu_atomizing_1983['L_b method']                   = np.nan
df_wu_atomizing_1983['L_b/d_0 page fig']             = np.nan
df_wu_atomizing_1983['L_b/d_0 transcription method'] = np.nan

df_wu_atomizing_1983['theta']                      = theta_wu_atomizing_1983
df_wu_atomizing_1983['theta stat error']           = e_theta_wu_atomizing_1983
df_wu_atomizing_1983['theta resolution']           = np.nan
df_wu_atomizing_1983['theta page fig']             = np.nan
df_wu_atomizing_1983['theta transcription method'] = np.nan

df_wu_atomizing_1983['D_10/d_0']                   = np.nan
df_wu_atomizing_1983['D_10/d_0 stat error']        = np.nan
df_wu_atomizing_1983['D_10/d_0 resolution']        = np.nan
df_wu_atomizing_1983['D_30/d_0']                   = np.nan
df_wu_atomizing_1983['D_30/d_0 stat error']        = np.nan
df_wu_atomizing_1983['D_30/d_0 resolution']        = np.nan
df_wu_atomizing_1983['D_32/d_0']                   = np.nan
df_wu_atomizing_1983['D_32/d_0 stat error']        = np.nan
df_wu_atomizing_1983['D_32/d_0 resolution']        = np.nan
df_wu_atomizing_1983['droplet x/d_0']              = np.nan
df_wu_atomizing_1983['D/d_0 page fig']             = np.nan
df_wu_atomizing_1983['D/d_0 transcription method'] = np.nan
df_wu_atomizing_1983['D/d_0 measurement method']   = np.nan
df_wu_atomizing_1983['v_d_bar/vp']                        = np.nan
df_wu_atomizing_1983['v_d_bar/vp stat error']             = np.nan
df_wu_atomizing_1983['v_d_bar/vp resolution']             = np.nan
df_wu_atomizing_1983['v_d_bar/vp page fig']               = np.nan
df_wu_atomizing_1983['v_d_bar/vp transcription method']   = np.nan
df_wu_atomizing_1983['e_p']                        = np.nan
df_wu_atomizing_1983['e_p page fig']               = np.nan
df_wu_atomizing_1983['e_p transcription method']   = np.nan

df_wu_atomizing_1983['x_trans/d_0']                      = np.nan
df_wu_atomizing_1983['x_trans/d_0 page fig']             = np.nan
df_wu_atomizing_1983['x_trans/d_0 transcription method'] = np.nan

df_wu_atomizing_1983['x_i/d_0']                      = np.nan
df_wu_atomizing_1983['x_i/d_0 stat error']           = np.nan
df_wu_atomizing_1983['x_i/d_0 resolution']           = np.nan
df_wu_atomizing_1983['x_i/d_0 page fig']             = np.nan
df_wu_atomizing_1983['x_i/d_0 transcription method'] = np.nan

df_wu_atomizing_1983['x_e/d_0']                      = np.nan
df_wu_atomizing_1983['x_e/d_0 stat error']           = np.nan
df_wu_atomizing_1983['x_e/d_0 resolution']           = np.nan
df_wu_atomizing_1983['x_e/d_0 page fig']             = np.nan
df_wu_atomizing_1983['x_e/d_0 transcription method'] = np.nan

# TODO: Add photos.
# photo: wu_atomizing_1983 p. A59
# details of photographic setup: wu_atomizing_1983 pp. 25-
df_wu_atomizing_1983['photo filename']   = np.nan
df_wu_atomizing_1983['exposure time']    = np.nan
df_wu_atomizing_1983['flash']            = np.nan
df_wu_atomizing_1983['x_low']            = np.nan
df_wu_atomizing_1983['x_mid']            = np.nan
df_wu_atomizing_1983['x_high']           = np.nan
df_wu_atomizing_1983['photo page fig']   = np.nan
df_wu_atomizing_1983['spectrum']         = np.nan
df_wu_atomizing_1983['lighting']         = np.nan
df_wu_atomizing_1983['background color'] = np.nan
df_wu_atomizing_1983['grid']             = np.nan
df_wu_atomizing_1983['distance']         = np.nan
df_wu_atomizing_1983['f-stop']           = np.nan
df_wu_atomizing_1983['focal length']     = np.nan
df_wu_atomizing_1983['camera model']     = np.nan
df_wu_atomizing_1983['sensor']           = np.nan

df_wu_atomizing_1983 = df_wu_atomizing_1983[cols]
summary_table(df_wu_atomizing_1983)
df_jet_breakup = pd.concat([df_jet_breakup, df_wu_atomizing_1983])

#####################################################
# shimizu_measurements_1984 / hiroyasu_breakup_1982 #
#####################################################

# Issues with this data:
# 1. In the spray angle plot, the lowest pressure case has a pressure of about 0.07 MPa (if you interpolate the position of the points like the data point extraction programs), but the text suggests the lowest pressure tested is 0.1 MPa. I used 0.1 MPa instead, and also rounded all the remaining pressures to be in line with what was discussed in the text.
# 2. The spray angle has an unclear definition. The plot says 2 \theta, which would imply that \theta is the half spray angle, but the results are most consistent \theta being the full spray angle when compared against other data. \theta here is assumed to be the full spray angle. It may have been labeled this way to increase the number of grid lines in the plot, which may have been placed automatically by software.

# TODO: Add additional information about the experimental setup from shimizu_measurements_1984.
# breakup length: figs. 8 and 13, p. 1712
# photos: fig. 15, p. 1713 (similar to hiroyasu_breakup_1982 but with scale; not the same photo)
# spray angle: fig. 22, p. 1714 (data already added)

# hiroyasu_breakup_1982:
# figs. 11 and 15
# Fig. 9 has a curve for breakup length but no data points. It can be used to determine L_b regime.
# Fig. 12 is likely duplicated by fig. 15, and it also does not seem to fully describe the differences between the fully developed data points.

shimizu_measurements_1984_csv = pd.read_csv('../data/shimizu_measurements_1984/shimizu_measurements_1984.csv', sep=',', header=0)

d_0_shimizu_measurements_1984            = shimizu_measurements_1984_csv['d_0 (mm)'] * 1e-3 # m
L_0s_shimizu_measurements_1984           = shimizu_measurements_1984_csv['L_0/d_0']
Ubar_0_shimizu_measurements_1984         = shimizu_measurements_1984_csv['Ubar_0 (m/s)']
P_atm_shimizu_measurements_1984          = shimizu_measurements_1984_csv['P_atm (MPa)'] * 1e6 # Pa
theta_shimizu_measurements_1984          = shimizu_measurements_1984_csv['theta (rad)'] / 2
xbavgs_shimizu_measurements_1984         = shimizu_measurements_1984_csv['xbavg (mm)'] * 1e-3 / d_0_shimizu_measurements_1984
regime_photo_shimizu_measurements_1984   = shimizu_measurements_1984_csv['regime photo']
regime_L_b_shimizu_measurements_1984     = shimizu_measurements_1984_csv['regime L_b']
photo_filename_shimizu_measurements_1984 = shimizu_measurements_1984_csv['photo filename']

i = 0

# water, unknown temperature
# gas is nitrogen, unknown temperature
rho_l = rho_water(T_std)
sigma = sigma_water(T_std)
nu_l  = nu_water(T_std)
P_v   = P_v_water(T_std)

We_l0_shimizu_measurements_1984      = np.zeros(len(Ubar_0_shimizu_measurements_1984))
Re_l0_shimizu_measurements_1984      = np.zeros(len(Ubar_0_shimizu_measurements_1984))
Fr_0_shimizu_measurements_1984       = np.zeros(len(Ubar_0_shimizu_measurements_1984))
rho_s_shimizu_measurements_1984      = np.zeros(len(Ubar_0_shimizu_measurements_1984))
nu_s_shimizu_measurements_1984       = np.zeros(len(Ubar_0_shimizu_measurements_1984))
K_c_shimizu_measurements_1984        = np.zeros(len(Ubar_0_shimizu_measurements_1984))
Ma_g_shimizu_measurements_1984       = np.zeros(len(Ubar_0_shimizu_measurements_1984))
x_low_shimizu_measurements_1984      = np.zeros(len(Ubar_0_shimizu_measurements_1984))
photo_page_fig_shimizu_measurements_1984   = []
spectrum_shimizu_measurements_1984         = []
background_color_shimizu_measurements_1984 = []
grid_shimizu_measurements_1984             = []
regime_turb_shimizu_measurements_1984      = []
for Ubar_0 in Ubar_0_shimizu_measurements_1984:
   d_0    = d_0_shimizu_measurements_1984[i]
   
   rho_g = rho_ideal_gas(P_atm_shimizu_measurements_1984[i], T_std, MW_N2)
   nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g # TODO: Change coefficients to be for N2.
   
   We_l0_shimizu_measurements_1984[i] = rho_l * Ubar_0**2 * d_0 / sigma
   Re_l0_shimizu_measurements_1984[i] = rho_l * Ubar_0 * d_0 / mu_l
   Fr_0_shimizu_measurements_1984[i]  = Ubar_0**2 / (g * d_0)
   rho_s_shimizu_measurements_1984[i] = rho_l / rho_g
   nu_s_shimizu_measurements_1984[i]  = nu_l / nu_g
   #K_c_shimizu_measurements_1984[i]   = (P_atm_shimizu_measurements_1984[i] - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_shimizu_measurements_1984[i]   = K_c(P_atm_shimizu_measurements_1984[i], P_v, rho_l, Ubar_0, K_L_sharpedged, friction_factor_smooth(Re_l0_shimizu_measurements_1984[i]), L_0s_shimizu_measurements_1984[i])
   Ma_g_shimizu_measurements_1984[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_N2)
   
   # TODO: Check that this is working correctly.
   if photo_filename_shimizu_measurements_1984[i] is np.nan:
      x_low_shimizu_measurements_1984[i] = np.nan
      photo_page_fig_shimizu_measurements_1984.append(np.nan)
      spectrum_shimizu_measurements_1984.append(np.nan)
      background_color_shimizu_measurements_1984.append(np.nan)
      grid_shimizu_measurements_1984.append(np.nan)
   else:
      x_low_shimizu_measurements_1984[i] = 0
      photo_page_fig_shimizu_measurements_1984.append('p. 71, fig. 6')
      spectrum_shimizu_measurements_1984.append('visible')
      background_color_shimizu_measurements_1984.append('light')
      grid_shimizu_measurements_1984.append(False)
   
   Re_l0 = Re_l0_shimizu_measurements_1984[i]
   if Re_l0 > Re_turb:
      regime_turb_shimizu_measurements_1984.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_shimizu_measurements_1984.append('transitional')
   else:
      regime_turb_shimizu_measurements_1984.append('laminar')
   
   i = i + 1

I_0_shimizu_measurements_1984 = I_fully_developed_array(Re_l0_shimizu_measurements_1984, Re_turb)
Lambda_r_s_shimizu_measurements_1984 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_shimizu_measurements_1984, Re_turb), 'v')

df_shimizu_measurements_1984 = pd.DataFrame({'We_l0': We_l0_shimizu_measurements_1984})
df_shimizu_measurements_1984['key']        = 'shimizu_measurements_1984'
df_shimizu_measurements_1984['alt key']    = 'hiroyasu_breakup_1982'
df_shimizu_measurements_1984['liquid']     = 'water' # p. 71L
df_shimizu_measurements_1984['gas']        = 'nitrogen' # p. 71L
df_shimizu_measurements_1984['degas']      = False
df_shimizu_measurements_1984['regulation'] = 'pressure regulator' # Based on shimizu_measurements_1984 fig. 2. Also see shimizu_measurements_1984 p. 1710R.

df_shimizu_measurements_1984['bend']              = False
df_shimizu_measurements_1984['flow straightener'] = False
df_shimizu_measurements_1984['contraction shape'] = 'sudden'
df_shimizu_measurements_1984['screen']            = False
df_shimizu_measurements_1984['c']                 = 10 # See fig. 3, p. 70. 3 mm / 0.3 mm
df_shimizu_measurements_1984['trip']              = False
df_shimizu_measurements_1984['L_r/d_0']           = 0
df_shimizu_measurements_1984['eps_r/d_0']         = np.nan

df_shimizu_measurements_1984['L_0/d_0']       = L_0s_shimizu_measurements_1984
df_shimizu_measurements_1984['roughness']     = 'smooth'
df_shimizu_measurements_1984['f page fig']    = np.nan
df_shimizu_measurements_1984['pipe material'] = np.nan
df_shimizu_measurements_1984['est f']         = True
df_shimizu_measurements_1984['check FD']      = False

df_shimizu_measurements_1984['t/d_0']       = np.nan # Can not be determined.
df_shimizu_measurements_1984['L_tip/d_0']   = np.nan
df_shimizu_measurements_1984['end checked'] = np.nan

df_shimizu_measurements_1984['orientation']         = 'down' # Possibly unclear; see fig. 2, p. 70. The view could be from the side or above. However, the bottom of p. 70L says "The screen wire detector could be moved up or down" which suggests the jet fires down.
df_shimizu_measurements_1984['vibration isolation'] = np.nan # Unclear.
df_shimizu_measurements_1984['d_chamber/d_0']       = 140e-3 / d_0_shimizu_measurements_1984 # shimizu_measurements_1984 p. 1710L

df_shimizu_measurements_1984['regime photo'] = regime_photo_shimizu_measurements_1984
df_shimizu_measurements_1984['regime L_b']   = regime_L_b_shimizu_measurements_1984
df_shimizu_measurements_1984['regime turb']  = regime_turb_shimizu_measurements_1984
df_shimizu_measurements_1984['est turb']     = True

#df_shimizu_measurements_1984['We_l0']      = We_l0_shimizu_measurements_1984
df_shimizu_measurements_1984['Re_l0']      = Re_l0_shimizu_measurements_1984
df_shimizu_measurements_1984['I_0']        = I_0_shimizu_measurements_1984
df_shimizu_measurements_1984['Fr_0']       = Fr_0_shimizu_measurements_1984
df_shimizu_measurements_1984['rho_s']      = rho_s_shimizu_measurements_1984
df_shimizu_measurements_1984['nu_s']       = nu_s_shimizu_measurements_1984
df_shimizu_measurements_1984['K_c']        = K_c_shimizu_measurements_1984
df_shimizu_measurements_1984['Ma_g']       = Ma_g_shimizu_measurements_1984
df_shimizu_measurements_1984['Lambda_r_s'] = Lambda_r_s_shimizu_measurements_1984
df_shimizu_measurements_1984['est rho_s']  = True
df_shimizu_measurements_1984['est nu_s']   = True

df_shimizu_measurements_1984['L_b/d_0']                      = xbavgs_shimizu_measurements_1984
df_shimizu_measurements_1984['L_b/d_0 stat error']           = np.nan
df_shimizu_measurements_1984['L_b/d_0 resolution']           = np.nan
df_shimizu_measurements_1984['L_b method']                   = 'electrical'
df_shimizu_measurements_1984['L_b/d_0 page fig']             = np.nan
df_shimizu_measurements_1984['L_b/d_0 transcription method'] = np.nan

df_shimizu_measurements_1984['theta']                      = theta_shimizu_measurements_1984
df_shimizu_measurements_1984['theta stat error']           = np.nan
df_shimizu_measurements_1984['theta resolution']           = np.nan
df_shimizu_measurements_1984['theta page fig']             = 'p. 73, fig. 15'
df_shimizu_measurements_1984['theta transcription method'] = 'figure'

df_shimizu_measurements_1984['D_10/d_0']                   = np.nan
df_shimizu_measurements_1984['D_10/d_0 stat error']        = np.nan
df_shimizu_measurements_1984['D_10/d_0 resolution']        = np.nan
df_shimizu_measurements_1984['D_30/d_0']                   = np.nan
df_shimizu_measurements_1984['D_30/d_0 stat error']        = np.nan
df_shimizu_measurements_1984['D_30/d_0 resolution']        = np.nan
df_shimizu_measurements_1984['D_32/d_0']                   = np.nan
df_shimizu_measurements_1984['D_32/d_0 stat error']        = np.nan
df_shimizu_measurements_1984['D_32/d_0 resolution']        = np.nan
df_shimizu_measurements_1984['droplet x/d_0']              = np.nan
df_shimizu_measurements_1984['D/d_0 page fig']             = np.nan
df_shimizu_measurements_1984['D/d_0 transcription method'] = np.nan
df_shimizu_measurements_1984['D/d_0 measurement method']   = np.nan
df_shimizu_measurements_1984['v_d_bar/vp']                        = np.nan
df_shimizu_measurements_1984['v_d_bar/vp stat error']             = np.nan
df_shimizu_measurements_1984['v_d_bar/vp resolution']             = np.nan
df_shimizu_measurements_1984['v_d_bar/vp page fig']               = np.nan
df_shimizu_measurements_1984['v_d_bar/vp transcription method']   = np.nan
df_shimizu_measurements_1984['e_p']                        = np.nan
df_shimizu_measurements_1984['e_p page fig']               = np.nan
df_shimizu_measurements_1984['e_p transcription method']   = np.nan

df_shimizu_measurements_1984['x_trans/d_0']                      = np.nan
df_shimizu_measurements_1984['x_trans/d_0 page fig']             = np.nan
df_shimizu_measurements_1984['x_trans/d_0 transcription method'] = np.nan

df_shimizu_measurements_1984['x_i/d_0']                      = np.nan
df_shimizu_measurements_1984['x_i/d_0 stat error']           = np.nan
df_shimizu_measurements_1984['x_i/d_0 resolution']           = np.nan
df_shimizu_measurements_1984['x_i/d_0 page fig']             = np.nan
df_shimizu_measurements_1984['x_i/d_0 transcription method'] = np.nan

df_shimizu_measurements_1984['x_e/d_0']                      = np.nan
df_shimizu_measurements_1984['x_e/d_0 stat error']           = np.nan
df_shimizu_measurements_1984['x_e/d_0 resolution']           = np.nan
df_shimizu_measurements_1984['x_e/d_0 page fig']             = np.nan
df_shimizu_measurements_1984['x_e/d_0 transcription method'] = np.nan

df_shimizu_measurements_1984['photo filename']   = photo_filename_shimizu_measurements_1984
df_shimizu_measurements_1984['exposure time']    = np.nan
df_shimizu_measurements_1984['flash']            = np.nan
df_shimizu_measurements_1984['x_low']            = x_low_shimizu_measurements_1984
df_shimizu_measurements_1984['x_mid']            = np.nan
df_shimizu_measurements_1984['x_high']           = np.nan
df_shimizu_measurements_1984['photo page fig']   = photo_page_fig_shimizu_measurements_1984
df_shimizu_measurements_1984['spectrum']         = spectrum_shimizu_measurements_1984
df_shimizu_measurements_1984['lighting']         = np.nan
df_shimizu_measurements_1984['background color'] = background_color_shimizu_measurements_1984
df_shimizu_measurements_1984['grid']             = grid_shimizu_measurements_1984
df_shimizu_measurements_1984['distance']         = np.nan
df_shimizu_measurements_1984['f-stop']           = np.nan
df_shimizu_measurements_1984['focal length']     = np.nan
df_shimizu_measurements_1984['camera model']     = np.nan
df_shimizu_measurements_1984['sensor']           = np.nan

df_shimizu_measurements_1984 = df_shimizu_measurements_1984[cols]
summary_table(df_shimizu_measurements_1984)
df_jet_breakup = pd.concat([df_jet_breakup, df_shimizu_measurements_1984])

######################
# arai_break-up_1985 #
######################

# DONE: No error analysis data provided. Making some assumptions.
# DONE: Add arai_breakup_1985 figs. 7 and 8 (Fig. 7 also reproduced by arai_similarity_1991 fig. 2.)

# TODO: figs. 7 and 8: Why is L/d = 40 different from L/d = 50? If this were caused by the density ratio/aerodynamic effects alone then I can't see how this would be possible. Also hard to believe it's turbulence too, as L/d = 40 should be fully developed, but maybe not in this case?

#rho_g = rho_ideal_gas(3e6, 30, MW_air)
#rho_l = rho_water(30)

#print rho_l / rho_g, (rho_l / rho_g - 1) / (rho_l / rho_g + 1)

# nozzle orifice diameter of 0.3 mm
# Liquid was water: p. IB/4/2: "( in this study, water was used )"
# p. IB/4/3: Re = 20000 is U = 53.4 m/s. This implies that nu = 8.01e-7, which returns almost exactly T = 30 C.
T_atm = 30 # C (assumed same temperature as the water)
rho_l_arai_break_up_1985 = rho_water(T_atm)
mu_l_arai_break_up_1985  = mu_water(T_atm)
nu_l_arai_break_up_1985  = nu_water(T_atm)
sigma_arai_break_up_1985 = sigma_water(T_atm)
d_0_arai_break_up_1985   = 0.3e-3; # m

arai_break_up_1985_csv = pd.read_csv('../data/arai_break-up_1985/arai_break-up_1985.csv', sep=',', header=0)

P_g_arai_break_up_1985        = arai_break_up_1985_csv['P_atm']
L_0s_arai_break_up_1985       = arai_break_up_1985_csv['L_0/d_0']
Re_l0_arai_break_up_1985      = arai_break_up_1985_csv['Re_l0']
L_bs_arai_break_up_1985       = arai_break_up_1985_csv['L_b/d_0']
theta_arai_break_up_1985      = arai_break_up_1985_csv['theta']
page_fig_arai_break_up_1985   = arai_break_up_1985_csv['page fig']
regime_L_b_arai_break_up_1985 = arai_break_up_1985_csv['regime L_b']

Oh_l0_arai_break_up_1985 = mu_l_arai_break_up_1985 / sqrt(rho_l_arai_break_up_1985 * sigma_arai_break_up_1985 * d_0_arai_break_up_1985)

We_l0_arai_break_up_1985 = (Re_l0_arai_break_up_1985 * Oh_l0_arai_break_up_1985)**2.
Fr_0_arai_break_up_1985  = ((Re_l0_arai_break_up_1985 * (mu_l_arai_break_up_1985 / rho_l_arai_break_up_1985))**2.) / (g * d_0_arai_break_up_1985**3)

I_0_arai_break_up_1985 = I_fully_developed_array(Re_l0_arai_break_up_1985, Re_turb)
Lambda_r_s_arai_break_up_1985 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_arai_break_up_1985, Re_turb), 'v')

i = 0
rho_s_arai_break_up_1985    = np.zeros(len(Re_l0_arai_break_up_1985))
nu_s_arai_break_up_1985     = np.zeros(len(Re_l0_arai_break_up_1985))
K_c_arai_break_up_1985      = np.zeros(len(Re_l0_arai_break_up_1985))
Ma_g_arai_break_up_1985     = np.zeros(len(Re_l0_arai_break_up_1985))
L_b_page_fig_arai_break_up_1985   = []
theta_page_fig_arai_break_up_1985 = []
regime_turb_arai_break_up_1985    = []
e_L_bs_arai_break_up_1985         = []
for P_g in P_g_arai_break_up_1985:
   Re_l0 = Re_l0_arai_break_up_1985[i]
   Ubar_0 = Re_l0 * nu_l_arai_break_up_1985 / d_0_arai_break_up_1985
   
   rho_g = rho_ideal_gas(P_g, T_std, MW_air)
   nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g
   
   rho_s_arai_break_up_1985[i] = rho_l / rho_g
   nu_s_arai_break_up_1985[i]  = nu_l / nu_g
   #K_c_arai_break_up_1985[i]   = (P_g - P_v_water(T_atm)) / (0.5 * rho_l * Ubar_0**2)
   K_c_arai_break_up_1985[i]   = K_c(P_g, P_v_water(T_atm), rho_l, Ubar_0, K_L_sharpedged, friction_factor_smooth(Re_l0), L_0s_arai_break_up_1985[i])
   Ma_g_arai_break_up_1985[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_atm), T_atm + T_zero, MW_air)
   
   if not(np.isnan(L_bs_arai_break_up_1985[i])):
      L_b_page_fig_arai_break_up_1985.append(page_fig_arai_break_up_1985[i])
      theta_page_fig_arai_break_up_1985.append(None)
      e_L_bs_arai_break_up_1985.append(0.2e-2 / d_0_arai_break_up_1985) # Following Phinney.
   elif not(np.isnan(theta_arai_break_up_1985[i])):
      L_b_page_fig_arai_break_up_1985.append(None)
      theta_page_fig_arai_break_up_1985.append(page_fig_arai_break_up_1985[i])
      e_L_bs_arai_break_up_1985.append(None)
   else:
      raise ValueError('Invalid data point for arai_break_up_1985. Needs to have either L_b or theta.')
   
   if Re_l0 > Re_turb:
      regime_turb_arai_break_up_1985.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_arai_break_up_1985.append('transitional')
   else:
      regime_turb_arai_break_up_1985.append('laminar')
   
   i = i + 1

df_arai_break_up_1985 = pd.DataFrame({'We_l0': We_l0_arai_break_up_1985})
df_arai_break_up_1985['key']        = 'arai_break-up_1985'
df_arai_break_up_1985['alt key']    = 'arai_similarity_1991' # TODO: Check this for any information that is missing from the 1985 paper.
df_arai_break_up_1985['liquid']     = 'water'
df_arai_break_up_1985['gas']        = 'air'
df_arai_break_up_1985['degas']      = False
df_arai_break_up_1985['regulation'] = 'pressure regulator' # Based on fig. 1.

df_arai_break_up_1985['bend']              = False
df_arai_break_up_1985['flow straightener'] = False
df_arai_break_up_1985['contraction shape'] = 'sudden'
df_arai_break_up_1985['screen']            = False
df_arai_break_up_1985['c']                 = 10 # See fig. 1.
df_arai_break_up_1985['trip']              = False
df_arai_break_up_1985['L_r/d_0']           = 0
df_arai_break_up_1985['eps_r/d_0']         = np.nan

df_arai_break_up_1985['L_0/d_0']       = L_0s_arai_break_up_1985
df_arai_break_up_1985['roughness']     = 'smooth'
df_arai_break_up_1985['f page fig']    = np.nan # TODO
df_arai_break_up_1985['pipe material'] = 'unknown, some acrylic' # p. IB/4/2: > Some of them were made of transparent acryl resin [...]
df_arai_break_up_1985['est f']         = True
df_arai_break_up_1985['check FD']      = False

df_arai_break_up_1985['t/d_0']       = np.nan # TODO: Figure this out based on the drawing of the nozzle?
df_arai_break_up_1985['L_tip/d_0']   = np.nan
df_arai_break_up_1985['end checked'] = np.nan

df_arai_break_up_1985['orientation']         = np.nan # Unclear; see fig. 1, p. IB/4/6. The view could be from the side or above.
df_arai_break_up_1985['vibration isolation'] = np.nan # Unclear.
df_arai_break_up_1985['d_chamber/d_0']       = np.nan # Unclear.

df_arai_break_up_1985['regime photo'] = np.nan
df_arai_break_up_1985['regime L_b']   = regime_L_b_arai_break_up_1985
df_arai_break_up_1985['regime turb']  = regime_turb_arai_break_up_1985
df_arai_break_up_1985['est turb']     = True

#df_arai_break_up_1985['We_l0']      = We_l0_arai_break_up_1985
df_arai_break_up_1985['Re_l0']      = Re_l0_arai_break_up_1985
df_arai_break_up_1985['I_0']        = I_0_arai_break_up_1985
df_arai_break_up_1985['Fr_0']       = Fr_0_arai_break_up_1985
df_arai_break_up_1985['rho_s']      = rho_s_arai_break_up_1985
df_arai_break_up_1985['nu_s']       = nu_s_arai_break_up_1985
df_arai_break_up_1985['K_c']        = K_c_arai_break_up_1985
df_arai_break_up_1985['Ma_g']       = Ma_g_arai_break_up_1985
df_arai_break_up_1985['Lambda_r_s'] = Lambda_r_s_arai_break_up_1985
df_arai_break_up_1985['est rho_s']  = True
df_arai_break_up_1985['est nu_s']   = True

df_arai_break_up_1985['L_b/d_0']                      = L_bs_arai_break_up_1985
df_arai_break_up_1985['L_b/d_0 stat error']           = np.nan
df_arai_break_up_1985['L_b/d_0 resolution']           = e_L_bs_arai_break_up_1985
df_arai_break_up_1985['L_b method']                   = 'electrical'
df_arai_break_up_1985['L_b/d_0 page fig']             = L_b_page_fig_arai_break_up_1985
df_arai_break_up_1985['L_b/d_0 transcription method'] = 'figure'

df_arai_break_up_1985['theta']                      = theta_arai_break_up_1985
df_arai_break_up_1985['theta stat error']           = np.nan # TODO
df_arai_break_up_1985['theta resolution']           = np.nan # TODO
df_arai_break_up_1985['theta page fig']             = theta_page_fig_arai_break_up_1985
df_arai_break_up_1985['theta transcription method'] = 'figure'

df_arai_break_up_1985['D_10/d_0']                   = np.nan
df_arai_break_up_1985['D_10/d_0 stat error']        = np.nan
df_arai_break_up_1985['D_10/d_0 resolution']        = np.nan
df_arai_break_up_1985['D_30/d_0']                   = np.nan
df_arai_break_up_1985['D_30/d_0 stat error']        = np.nan
df_arai_break_up_1985['D_30/d_0 resolution']        = np.nan
df_arai_break_up_1985['D_32/d_0']                   = np.nan
df_arai_break_up_1985['D_32/d_0 stat error']        = np.nan
df_arai_break_up_1985['D_32/d_0 resolution']        = np.nan
df_arai_break_up_1985['droplet x/d_0']              = np.nan
df_arai_break_up_1985['D/d_0 page fig']             = np.nan
df_arai_break_up_1985['D/d_0 transcription method'] = np.nan
df_arai_break_up_1985['D/d_0 measurement method']   = np.nan
df_arai_break_up_1985['v_d_bar/vp']                        = np.nan
df_arai_break_up_1985['v_d_bar/vp stat error']             = np.nan
df_arai_break_up_1985['v_d_bar/vp resolution']             = np.nan
df_arai_break_up_1985['v_d_bar/vp page fig']               = np.nan
df_arai_break_up_1985['v_d_bar/vp transcription method']   = np.nan
df_arai_break_up_1985['e_p']                        = np.nan
df_arai_break_up_1985['e_p page fig']               = np.nan
df_arai_break_up_1985['e_p transcription method']   = np.nan

df_arai_break_up_1985['x_trans/d_0']                      = np.nan
df_arai_break_up_1985['x_trans/d_0 page fig']             = np.nan
df_arai_break_up_1985['x_trans/d_0 transcription method'] = np.nan

df_arai_break_up_1985['x_i/d_0']                      = np.nan
df_arai_break_up_1985['x_i/d_0 stat error']           = np.nan
df_arai_break_up_1985['x_i/d_0 resolution']           = np.nan
df_arai_break_up_1985['x_i/d_0 page fig']             = np.nan
df_arai_break_up_1985['x_i/d_0 transcription method'] = np.nan

df_arai_break_up_1985['x_e/d_0']                      = np.nan
df_arai_break_up_1985['x_e/d_0 stat error']           = np.nan
df_arai_break_up_1985['x_e/d_0 resolution']           = np.nan
df_arai_break_up_1985['x_e/d_0 page fig']             = np.nan
df_arai_break_up_1985['x_e/d_0 transcription method'] = np.nan

df_arai_break_up_1985['photo filename']   = np.nan
df_arai_break_up_1985['exposure time']    = np.nan
df_arai_break_up_1985['flash']            = np.nan
df_arai_break_up_1985['x_low']            = np.nan
df_arai_break_up_1985['x_mid']            = np.nan
df_arai_break_up_1985['x_high']           = np.nan
df_arai_break_up_1985['photo page fig']   = np.nan
df_arai_break_up_1985['spectrum']         = np.nan
df_arai_break_up_1985['lighting']         = np.nan
df_arai_break_up_1985['background color'] = np.nan
df_arai_break_up_1985['grid']             = np.nan
df_arai_break_up_1985['distance']         = np.nan
df_arai_break_up_1985['f-stop']           = np.nan
df_arai_break_up_1985['focal length']     = np.nan
df_arai_break_up_1985['camera model']     = np.nan
df_arai_break_up_1985['sensor']           = np.nan

df_arai_break_up_1985 = df_arai_break_up_1985[cols]
summary_table(df_arai_break_up_1985)
df_jet_breakup = pd.concat([df_jet_breakup, df_arai_break_up_1985])

#########################
## debler_break-up_1988 #
#########################

# TODO: Add transitional regimes (R2F, F2S, S2A).

##load -ascii ../../data/debler_break-up_1988/ar_67_7.dat
##load -ascii ../../data/debler_break-up_1988/ar_78_9.dat
##load -ascii ../../data/debler_break-up_1988/ar_96.dat

#############################################
# ruff_structure_1990 / ruff_structure_1989 #
#############################################

# TODO: Get Tu at exit from paper.
# p. 39, tab. 2.4: test conditions table

ruff_structure_1990_csv = pd.read_csv('../data/ruff_structure_1990/ruff_structure_1990.csv', sep=',', header=0)

page_fig_ruff_structure_1990     = ruff_structure_1990_csv['page fig']
liquid_ruff_structure_1990       = ruff_structure_1990_csv['liquid']
d_0_ruff_structure_1990          = ruff_structure_1990_csv['d (mm)'] * 1e-3 # m
mdot_ruff_structure_1990         = ruff_structure_1990_csv['mdot (kg/s)']
Re_l0_ruff_structure_1990        = ruff_structure_1990_csv['Re_l0']
We_l0_ruff_structure_1990        = ruff_structure_1990_csv['We_l0']
We_g0_ruff_structure_1990        = ruff_structure_1990_csv['We_g0']
Oh_l0_ruff_structure_1990        = ruff_structure_1990_csv['Oh_l0']
regime_photo_ruff_structure_1990 = ruff_structure_1990_csv['regime photo']
theta_ruff_structure_1990        = ruff_structure_1990_csv['theta (deg)'] * (np.pi / 180)
P_atm_ruff_structure_1990        = ruff_structure_1990_csv['P_atm (kPa)'] * 1e3 # Pa
dP_ruff_structure_1990           = ruff_structure_1990_csv['dP (KPa)'] * 1e3 # Pa

#rho_s_ruff_structure_1990 = We_l0_ruff_structure_1990 / We_g0_ruff_structure_1990
#print rho_s_ruff_structure_1990

T_g_ruff_structure_1990 = 298 - T_zero # C, ruff_structure_1989 p. 903L, tab. 1

nu_l  = nu_water(T_g_ruff_structure_1990)
rho_l = rho_water(T_g_ruff_structure_1990)
rho_g = rho_ideal_gas(P_atm, T_g_ruff_structure_1990, MW_air)
nu_g  = mu_ideal_gas(T_g_ruff_structure_1990, mu_0_air, T_0_air, C_air) / rho_g

i = 0
Fr_0_ruff_structure_1990  = np.zeros(len(Re_l0_ruff_structure_1990))
rho_s_ruff_structure_1990 = np.zeros(len(Re_l0_ruff_structure_1990))
nu_s_ruff_structure_1990  = np.zeros(len(Re_l0_ruff_structure_1990))
K_c_ruff_structure_1990   = np.zeros(len(Re_l0_ruff_structure_1990))
Ma_g_ruff_structure_1990  = np.zeros(len(Re_l0_ruff_structure_1990))
regime_turb_ruff_structure_1990    = []
for Re_l0 in Re_l0_ruff_structure_1990:
   Re_l0 = Re_l0_ruff_structure_1990[i]
   #print Re_l0
   d_0   = d_0_ruff_structure_1990[i]
   P_atm = P_atm_ruff_structure_1990[i]
   
   Ubar_0 = Re_l0 * nu_l / d_0
   
   Fr_0_ruff_structure_1990[i]  = Ubar_0**2 / (g * d_0)
   rho_s_ruff_structure_1990[i] = rho_l / rho_g
   nu_s_ruff_structure_1990[i]  = nu_l / nu_g
   #K_c_ruff_structure_1990[i]   = (P_g - P_v_water(T_g_ruff_structure_1990)) / (0.5 * rho_l * Ubar_0**2)
   K_c_ruff_structure_1990[i]   = K_c(P_g, P_v_water(T_g_ruff_structure_1990), rho_l, Ubar_0, K_L_slightlyrounded, friction_factor_smooth(Re_l0), 41.)
   Ma_g_ruff_structure_1990[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_g_ruff_structure_1990), T_g_ruff_structure_1990 + T_zero, MW_air)
   
   if Re_l0 > Re_turb:
      regime_turb_ruff_structure_1990.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_ruff_structure_1990.append('transitional')
   else:
      regime_turb_ruff_structure_1990.append('laminar')
   
   i = i + 1

I_0_ruff_structure_1990 = I_fully_developed_array(Re_l0_ruff_structure_1990, Re_turb)
Lambda_r_s_ruff_structure_1990 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_ruff_structure_1990, Re_turb), 'v')

df_ruff_structure_1990 = pd.DataFrame({'We_l0': We_l0_ruff_structure_1990})
df_ruff_structure_1990['key']        = 'ruff_structure_1990'
df_ruff_structure_1990['alt key']    = 'sallam_properties_2002'
df_ruff_structure_1990['liquid']     = 'water'
df_ruff_structure_1990['gas']        = 'air'
df_ruff_structure_1990['degas']      = np.nan
df_ruff_structure_1990['regulation'] = np.nan

df_ruff_structure_1990['bend']              = np.nan
df_ruff_structure_1990['flow straightener'] = True
df_ruff_structure_1990['contraction shape'] = 'somewhat sudden' # p. 12, fig. 2.2
df_ruff_structure_1990['screen']            = False
df_ruff_structure_1990['c']                 = 35e-3 / d_0_ruff_structure_1990 # p. 11
df_ruff_structure_1990['trip']              = False
df_ruff_structure_1990['L_r/d_0']           = 0.
df_ruff_structure_1990['eps_r/d_0']         = np.nan

df_ruff_structure_1990['L_0/d_0']       = 41.
df_ruff_structure_1990['roughness']     = 'smooth'
df_ruff_structure_1990['f page fig']    = np.nan
df_ruff_structure_1990['pipe material'] = np.nan
df_ruff_structure_1990['est f']         = True # TODO: Use TI data given.
df_ruff_structure_1990['check FD']      = np.nan

df_ruff_structure_1990['t/d_0']       = np.nan
df_ruff_structure_1990['L_tip/d_0']   = np.nan
df_ruff_structure_1990['end checked'] = np.nan

df_ruff_structure_1990['orientation']         = 'down'
df_ruff_structure_1990['vibration isolation'] = np.nan # Unclear.
df_ruff_structure_1990['d_chamber/d_0']       = np.inf

df_ruff_structure_1990['regime photo'] = regime_photo_ruff_structure_1990
df_ruff_structure_1990['regime L_b']   = np.nan
df_ruff_structure_1990['regime turb']  = regime_turb_ruff_structure_1990
df_ruff_structure_1990['est turb']     = True

#df_ruff_structure_1990['We_l0']      = We_l0_ruff_structure_1990
df_ruff_structure_1990['Re_l0']      = Re_l0_ruff_structure_1990
df_ruff_structure_1990['I_0']        = I_0_ruff_structure_1990
df_ruff_structure_1990['Fr_0']       = Fr_0_ruff_structure_1990
df_ruff_structure_1990['rho_s']      = rho_s_ruff_structure_1990
df_ruff_structure_1990['nu_s']       = nu_s_ruff_structure_1990
df_ruff_structure_1990['K_c']        = K_c_ruff_structure_1990
df_ruff_structure_1990['Ma_g']       = Ma_g_ruff_structure_1990
df_ruff_structure_1990['Lambda_r_s'] = Lambda_r_s_ruff_structure_1990
df_ruff_structure_1990['est rho_s']  = True
df_ruff_structure_1990['est nu_s']   = True

df_ruff_structure_1990['L_b/d_0']                      = np.nan
df_ruff_structure_1990['L_b/d_0 stat error']           = np.nan
df_ruff_structure_1990['L_b/d_0 resolution']           = np.nan
df_ruff_structure_1990['L_b method']                   = np.nan
df_ruff_structure_1990['L_b/d_0 page fig']             = np.nan
df_ruff_structure_1990['L_b/d_0 transcription method'] = np.nan

df_ruff_structure_1990['theta']                      = theta_ruff_structure_1990
df_ruff_structure_1990['theta stat error']           = np.nan
df_ruff_structure_1990['theta resolution']           = np.nan
df_ruff_structure_1990['theta page fig']             = np.nan
df_ruff_structure_1990['theta transcription method'] = 'table'

df_ruff_structure_1990['D_10/d_0']                   = np.nan
df_ruff_structure_1990['D_10/d_0 stat error']        = np.nan
df_ruff_structure_1990['D_10/d_0 resolution']        = np.nan
df_ruff_structure_1990['D_30/d_0']                   = np.nan
df_ruff_structure_1990['D_30/d_0 stat error']        = np.nan
df_ruff_structure_1990['D_30/d_0 resolution']        = np.nan
df_ruff_structure_1990['D_32/d_0']                   = np.nan
df_ruff_structure_1990['D_32/d_0 stat error']        = np.nan
df_ruff_structure_1990['D_32/d_0 resolution']        = np.nan
df_ruff_structure_1990['droplet x/d_0']              = np.nan
df_ruff_structure_1990['D/d_0 page fig']             = np.nan
df_ruff_structure_1990['D/d_0 transcription method'] = np.nan
df_ruff_structure_1990['D/d_0 measurement method']   = np.nan
df_ruff_structure_1990['v_d_bar/vp']                        = np.nan
df_ruff_structure_1990['v_d_bar/vp stat error']             = np.nan
df_ruff_structure_1990['v_d_bar/vp resolution']             = np.nan
df_ruff_structure_1990['v_d_bar/vp page fig']               = np.nan
df_ruff_structure_1990['v_d_bar/vp transcription method']   = np.nan
df_ruff_structure_1990['e_p']                        = np.nan
df_ruff_structure_1990['e_p page fig']               = np.nan
df_ruff_structure_1990['e_p transcription method']   = np.nan

df_ruff_structure_1990['x_trans/d_0']                      = np.nan
df_ruff_structure_1990['x_trans/d_0 page fig']             = np.nan
df_ruff_structure_1990['x_trans/d_0 transcription method'] = np.nan

df_ruff_structure_1990['x_i/d_0']                      = np.nan
df_ruff_structure_1990['x_i/d_0 stat error']           = np.nan
df_ruff_structure_1990['x_i/d_0 resolution']           = np.nan
df_ruff_structure_1990['x_i/d_0 page fig']             = np.nan
df_ruff_structure_1990['x_i/d_0 transcription method'] = np.nan

df_ruff_structure_1990['x_e/d_0']                      = np.nan
df_ruff_structure_1990['x_e/d_0 stat error']           = np.nan
df_ruff_structure_1990['x_e/d_0 resolution']           = np.nan
df_ruff_structure_1990['x_e/d_0 page fig']             = np.nan
df_ruff_structure_1990['x_e/d_0 transcription method'] = np.nan

# TODO: Add photos.
df_ruff_structure_1990['photo filename']   = np.nan
df_ruff_structure_1990['exposure time']    = np.nan
df_ruff_structure_1990['flash']            = np.nan
df_ruff_structure_1990['x_low']            = np.nan
df_ruff_structure_1990['x_mid']            = np.nan
df_ruff_structure_1990['x_high']           = np.nan
df_ruff_structure_1990['photo page fig']   = np.nan
df_ruff_structure_1990['spectrum']         = np.nan
df_ruff_structure_1990['lighting']         = np.nan
df_ruff_structure_1990['background color'] = np.nan
df_ruff_structure_1990['grid']             = np.nan
df_ruff_structure_1990['distance']         = np.nan
df_ruff_structure_1990['f-stop']           = np.nan
df_ruff_structure_1990['focal length']     = np.nan
df_ruff_structure_1990['camera model']     = np.nan
df_ruff_structure_1990['sensor']           = np.nan

df_ruff_structure_1990 = df_ruff_structure_1990[cols]
summary_table(df_ruff_structure_1990)
df_jet_breakup = pd.concat([df_jet_breakup, df_ruff_structure_1990])

############################
# tseng_near-injector_1991 #
############################

# TODO: Add.
# p. 12: > Water was used as the test liquid, injected vertically downward into still air within a pressure vessel (Mc Daniel Tank Manufacturing Company) having a 1.5 m diameter and 4.5 m height
# p. 12: > centrifugal water pump
# p. 14: > honeycomb flow straightener (1.6 mm cells, 25 mm long), two screens to calm the flow (0.018 mm diameter wire, 16 x 16 square mesh), followed by a smooth converging section from the 35 mm diameter water supply tube to the 9.5 mm exit diameter.
# p. 15, fig. 2.2: flow straightener
# p. 35, tab. 2.3: test conditions
# p. 35: > still air at various pressures and $298\pm2$K; in atomization breakup regime

##################
# wu_liquid_1992 #
##################

# Issues with this data:
# 1. The 6.4 mm diameter appears to be labeled as 5.0 mm in table B.5 for water, glycerol, and n-heptane. I assume it's supposed to be 6.4 mm in table B.5 as it's that way in all other places this data is. sallam_properties_2002 p. 133, tab. B.4 seems to have the diameters labeled correctly, along with the data at seemingly higher precision. I assume this is transcribed from a plot and that the original data in the Wu dissertation is right. The numbers are similar to the original dissertation's aside from the diameter and integral scale.
# 2. Page 78 suggests they used \Lambda = d_0 / 8 to calculate \Lambda. Did they use the wrong \Lambda when scaling x_i/\Lambda for the 6.4 mm case? Check the SMD/\Lambda plot in wu_liquid_1992 and the x_i/d_0 plot in wu_onset_1995.
# 3. The We_l0 for the heptane data in wu_liquid_1992 seems inconsistent with the corresponding data in wu_onset_1995. It appears that wu_liquid_1992 has the correct Weber numbers.
# 4. For wu_onset_1995, the new data points do not have any ambient gas information. It is assumed these points are in air at atmospheric pressure and temperature.
# 5. No statistical/repeatability uncertainties are given for x_i, so the percentage error is assumed to be the same as for SMD.

# Uncertainty: p. 126 (pdf p. 143), p. 129 (pdf p. 146)
# p. 126 gives minimum 2.3% for all quantities (which presumably includes x_i)

# Tables: wu_liquid_1992 pp. 134-139 (pdf p. 151-139)

wu_liquid_1992_csv = pd.read_csv('../data/wu_liquid_1992/wu_liquid_1992.csv', sep=',', header=0)

key_wu_liquid_1992            = wu_liquid_1992_csv['key']
alt_key_wu_liquid_1992        = wu_liquid_1992_csv['alt key']
x_i_page_fig_wu_liquid_1992   = wu_liquid_1992_csv['x_i graph page fig']
SMD_i_page_fig_wu_liquid_1992 = wu_liquid_1992_csv['SMD_i graph page fig']
liquid_wu_liquid_1992         = wu_liquid_1992_csv['liquid']
gas_wu_liquid_1992            = wu_liquid_1992_csv['gas']
L_0s_wu_liquid_1992           = wu_liquid_1992_csv['L_0/d_0']
d_0_wu_liquid_1992            = wu_liquid_1992_csv['d (mm)'] * 1e-3 # m
P_g_wu_liquid_1992            = wu_liquid_1992_csv['P_atm (atm)'] * P_atm # Pa
rho_s_wu_liquid_1992          = wu_liquid_1992_csv['rho_l/rho_g']
rho_l_wu_liquid_1992          = wu_liquid_1992_csv['rho_l']
mu_l_wu_liquid_1992           = wu_liquid_1992_csv['mu_l']
sigma_wu_liquid_1992          = wu_liquid_1992_csv['sigma']
rho_g_wu_liquid_1992          = wu_liquid_1992_csv['rho_g']
Ubar_0_wu_liquid_1992         = wu_liquid_1992_csv['Ubar_0 (m/s)']
We_l0_wu_liquid_1992          = wu_liquid_1992_csv['We_l0']
x_is_wu_liquid_1992           = wu_liquid_1992_csv['x_i/d_0']
D_32_is_wu_liquid_1992        = wu_liquid_1992_csv['D_{32,i}/d_0']
v_d_s_wu_liquid_1992          = wu_liquid_1992_csv['vtilde_p/Ubar_0']
regime_photo_wu_liquid_1992   = wu_liquid_1992_csv['regime photo']
photo_filename_wu_liquid_1992 = wu_liquid_1992_csv['photo filename']

T_wu_liquid_1992 = 298 # K, \pm 3; see p. 34

i = 0
for We_l0, d_0, Ubar_0, rho_l, sigma in zip(We_l0_wu_liquid_1992, d_0_wu_liquid_1992, Ubar_0_wu_liquid_1992, rho_l_wu_liquid_1992, sigma_wu_liquid_1992):
   if np.isnan(We_l0):
      We_l0 = rho_l * Ubar_0**2. * d_0 / sigma
      We_l0_wu_liquid_1992[i] = We_l0

i = 0
Re_l0_wu_liquid_1992           = np.zeros(len(We_l0_wu_liquid_1992))
Fr_0_wu_liquid_1992            = np.zeros(len(We_l0_wu_liquid_1992))
rho_s_wu_liquid_1992           = np.zeros(len(We_l0_wu_liquid_1992))
nu_s_wu_liquid_1992            = np.zeros(len(We_l0_wu_liquid_1992))
K_c_wu_liquid_1992             = np.zeros(len(We_l0_wu_liquid_1992))
Ma_g_wu_liquid_1992            = np.zeros(len(We_l0_wu_liquid_1992))
d_chambers_wu_liquid_1992      = np.zeros(len(We_l0_wu_liquid_1992))
x_is_wu_liquid_1992_resolution = np.zeros(len(We_l0_wu_liquid_1992))
regime_turb_wu_liquid_1992     = []
for We_l0 in We_l0_wu_liquid_1992:
   rho_l  = rho_l_wu_liquid_1992[i]
   rho_g  = rho_g_wu_liquid_1992[i]
   mu_l   = mu_l_wu_liquid_1992[i]
   P_g    = P_g_wu_liquid_1992[i]
   d_0    = d_0_wu_liquid_1992[i]
   sigma  = sigma_wu_liquid_1992[i]
   
   gas = gas_wu_liquid_1992[i]
   
   if gas == 'air':
      mu_g = mu_ideal_gas(T_g_ruff_structure_1990, mu_0_air, T_0_air, C_air)
      nu_g = mu_g / rho_g
   elif gas == 'helium':
      mu_g = np.interp(T_wu_liquid_1992 - T_zero, np.array([0, 20, 50]), np.array([1.87, 1.96, 2.10])) * 1e-5
      nu_g = mu_g / rho_g
   elif gas == 'freon 12':
      nu_g = np.interp(T_wu_liquid_1992 - T_zero, np.array([0, 10, 20, 30]), np.array([0.214, 0.203, 0.198, 0.194])) * 1e-6
   else:
      raise ValueError('Invalid gas for wu_liquid_1992:', gas)
   
   nu_l = mu_l / rho_l
   
   Ubar_0 = sqrt((We_l0 * sigma) / (rho_l * d_0))
   #print Ubar_0, d_0
   
   # TODO: Generalize to other fluids in the future.
   P_v = P_v_water(T_wu_liquid_1992 - T_zero)
   
   Re_l0_wu_liquid_1992[i] = Ubar_0 * d_0 / nu_l
   Fr_0_wu_liquid_1992[i]  = Ubar_0**2 / (g * d_0)
   rho_s_wu_liquid_1992[i] = rho_l / rho_g
   nu_s_wu_liquid_1992[i]  = nu_l / nu_g
   #K_c_wu_liquid_1992[i]   = (P_g - P_v) / (0.5 * rho_l * Ubar_0**2)
   K_c_wu_liquid_1992[i]   = K_c(P_g, P_v, rho_l, Ubar_0, K_L_slightlyrounded, friction_factor_smooth(Re_l0_wu_liquid_1992[i]), L_0s_wu_liquid_1992[i])
   Ma_g_wu_liquid_1992[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_std), T_std + T_zero, MW_N2)
   
   if P_g == 1:
      d_chambers_wu_liquid_1992[i] = np.inf
   else:
      d_chambers_wu_liquid_1992[i] = 300e-3 / d_0 # p. 18
   
   if Re_l0 > Re_turb:
      regime_turb_wu_liquid_1992.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_wu_liquid_1992.append('transitional')
   else:
      regime_turb_wu_liquid_1992.append('laminar')
   
   # TODO: Add x_i errors.
   # p. 126 (pdf p. 143): > The high magnification of the system (110-320X) minimized this uncertainty: it is estimated to be 2.3%.
   # p. 127 (pdf p. 144): > measurements from the tip of the reference pin to the edge of the measuring volume have a maximum position uncertainty of ±70 $\mu$m in both x- and y-directions, with substantially lower uncertainties near the center of the measuring volume.
   # pp. 128-129 does not discuss location repeatability errors. These could be estimated under certain assumptions. For example, you could consider multiple data points as averages of 20 data points (wu_liquid_1992 p. 128 says each point is at least 20 data points, and they probably are mostly at the lower end). Then you can estimate the standard deviation for the repeatability based on the standard deviation of data points grouped closely together in We_l0 * I_0^3. This assumes that no other variables are involved.
   
   x_is_wu_liquid_1992_resolution[i] = 0.023 * x_is_wu_liquid_1992[i] + 70e-6 / d_0
   
   i = i + 1

#i = 0
#x_is_window_array   = []
#We_l0_window_center = 6e3
#We_l0_window_width  = 2e3
#for We_l0 in We_l0_wu_liquid_1992:
   #if (We_l0 < (We_l0_window_center + We_l0_window_width)) and (We_l0 > (We_l0_window_center - We_l0_window_width)):
      ##print We_l0, x_is_wu_liquid_1992[i]
      #x_is_window_array.append(x_is_wu_liquid_1992[i])
   
   #i = i + 1

## Assume that the ratio of standard deviation to mean is constant.
#x_is_std  = std(x_is_window_array)
#x_is_mean = mean(x_is_window_array)

##print x_is_std / x_is_mean, len(x_is_window_array)

i = 0
x_is_wu_liquid_1992_stat_error = np.zeros(len(We_l0_wu_liquid_1992))
for x_is in x_is_wu_liquid_1992:
   #x_is_wu_liquid_1992_stat_error[i] = (x_is_std / x_is_mean) * z * x_is
   #x_is_wu_liquid_1992_stat_error[i] = 0
   
   # Assumed to be the same as SMD.
   x_is_wu_liquid_1992_stat_error[i] = 0.33 * x_is
   
   i = i + 1

#print x_is_wu_liquid_1992_stat_error / x_is_wu_liquid_1992

I_0_wu_liquid_1992 = I_fully_developed_array(Re_l0_wu_liquid_1992, Re_turb)
Lambda_r_s_wu_liquid_1992 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_wu_liquid_1992, Re_turb), 'v')

df_wu_liquid_1992 = pd.DataFrame({'We_l0': We_l0_wu_liquid_1992})
df_wu_liquid_1992['key']        = key_wu_liquid_1992
df_wu_liquid_1992['alt key']    = alt_key_wu_liquid_1992
df_wu_liquid_1992['liquid']     = liquid_wu_liquid_1992
df_wu_liquid_1992['gas']        = gas_wu_liquid_1992
df_wu_liquid_1992['degas']      = np.nan
df_wu_liquid_1992['regulation'] = 'pressure regulator' # wu_liquid_1992 p. 50 fig. 4-2: pressure regulator

df_wu_liquid_1992['bend']              = False
df_wu_liquid_1992['flow straightener'] = False
df_wu_liquid_1992['contraction shape'] = 'mostly sudden'
df_wu_liquid_1992['screen']            = False
df_wu_liquid_1992['c']                 = np.nan # TODO: Is this written anywhere? wu_effects_1995 p. 178 says it is about 100.
df_wu_liquid_1992['trip']              = False
df_wu_liquid_1992['L_r/d_0']           = 0
df_wu_liquid_1992['eps_r/d_0']         = np.nan

df_wu_liquid_1992['L_0/d_0']       = L_0s_wu_liquid_1992
df_wu_liquid_1992['roughness']     = 'smooth'
df_wu_liquid_1992['f page fig']    = np.nan
df_wu_liquid_1992['pipe material'] = np.nan # Not listed?
df_wu_liquid_1992['est f']         = True
df_wu_liquid_1992['check FD']      = False

df_wu_liquid_1992['t/d_0']       = np.nan
df_wu_liquid_1992['L_tip/d_0']   = np.nan
df_wu_liquid_1992['end checked'] = np.nan

df_wu_liquid_1992['orientation']         = 'down' # wu_liquid_1992 p. 51 fig. 4-2: suggests the jet is directed down
df_wu_liquid_1992['vibration isolation'] = np.nan # TODO
df_wu_liquid_1992['d_chamber/d_0']       = d_chambers_wu_liquid_1992

df_wu_liquid_1992['regime photo'] = regime_photo_wu_liquid_1992
df_wu_liquid_1992['regime L_b']   = np.nan
df_wu_liquid_1992['regime turb']  = regime_turb_wu_liquid_1992
df_wu_liquid_1992['est turb']     = True

#df_wu_liquid_1992['We_l0']      = We_l0_wu_liquid_1992
df_wu_liquid_1992['Re_l0']      = Re_l0_wu_liquid_1992
df_wu_liquid_1992['I_0']        = I_0_wu_liquid_1992
df_wu_liquid_1992['Fr_0']       = Fr_0_wu_liquid_1992
df_wu_liquid_1992['rho_s']      = rho_s_wu_liquid_1992
df_wu_liquid_1992['nu_s']       = nu_s_wu_liquid_1992
df_wu_liquid_1992['K_c']        = K_c_wu_liquid_1992
df_wu_liquid_1992['Ma_g']       = Ma_g_wu_liquid_1992
df_wu_liquid_1992['Lambda_r_s'] = Lambda_r_s_wu_liquid_1992
df_wu_liquid_1992['est rho_s']  = False
df_wu_liquid_1992['est nu_s']   = True

df_wu_liquid_1992['L_b/d_0']                      = np.nan
df_wu_liquid_1992['L_b/d_0 stat error']           = np.nan
df_wu_liquid_1992['L_b/d_0 resolution']           = np.nan
df_wu_liquid_1992['L_b method']                   = np.nan
df_wu_liquid_1992['L_b/d_0 page fig']             = np.nan
df_wu_liquid_1992['L_b/d_0 transcription method'] = np.nan

df_wu_liquid_1992['theta']                      = np.nan
df_wu_liquid_1992['theta stat error']           = np.nan
df_wu_liquid_1992['theta resolution']           = np.nan
df_wu_liquid_1992['theta page fig']             = np.nan
df_wu_liquid_1992['theta transcription method'] = np.nan

df_wu_liquid_1992['D_10/d_0']                   = np.nan
df_wu_liquid_1992['D_10/d_0 stat error']        = np.nan
df_wu_liquid_1992['D_10/d_0 resolution']        = np.nan
df_wu_liquid_1992['D_30/d_0']                   = np.nan
df_wu_liquid_1992['D_30/d_0 stat error']        = np.nan
df_wu_liquid_1992['D_30/d_0 resolution']        = np.nan
df_wu_liquid_1992['D_32/d_0']                   = D_32_is_wu_liquid_1992
df_wu_liquid_1992['D_32/d_0 stat error']        = np.nan # TODO
df_wu_liquid_1992['D_32/d_0 resolution']        = np.nan # TODO
df_wu_liquid_1992['droplet x/d_0']              = np.nan
df_wu_liquid_1992['D/d_0 page fig']             = np.nan
df_wu_liquid_1992['D/d_0 transcription method'] = np.nan
df_wu_liquid_1992['D/d_0 measurement method']   = np.nan
df_wu_liquid_1992['v_d_bar/vp']                      = v_d_s_wu_liquid_1992
df_wu_liquid_1992['v_d_bar/vp stat error']           = np.nan # TODO
df_wu_liquid_1992['v_d_bar/vp resolution']           = np.nan # TODO
df_wu_liquid_1992['v_d_bar/vp page fig']             = np.nan # TODO
df_wu_liquid_1992['v_d_bar/vp transcription method'] = np.nan # TODO
df_wu_liquid_1992['e_p']                        = np.nan
df_wu_liquid_1992['e_p page fig']               = np.nan
df_wu_liquid_1992['e_p transcription method']   = np.nan

df_wu_liquid_1992['x_trans/d_0']                      = np.nan
df_wu_liquid_1992['x_trans/d_0 page fig']             = np.nan
df_wu_liquid_1992['x_trans/d_0 transcription method'] = np.nan

df_wu_liquid_1992['x_i/d_0']                      = x_is_wu_liquid_1992
df_wu_liquid_1992['x_i/d_0 stat error']           = x_is_wu_liquid_1992_stat_error
df_wu_liquid_1992['x_i/d_0 resolution']           = x_is_wu_liquid_1992_resolution
df_wu_liquid_1992['x_i/d_0 page fig']             = np.nan # TODO
df_wu_liquid_1992['x_i/d_0 transcription method'] = np.nan # TODO

df_wu_liquid_1992['x_e/d_0']                      = np.nan
df_wu_liquid_1992['x_e/d_0 stat error']           = np.nan
df_wu_liquid_1992['x_e/d_0 resolution']           = np.nan
df_wu_liquid_1992['x_e/d_0 page fig']             = np.nan
df_wu_liquid_1992['x_e/d_0 transcription method'] = np.nan

df_wu_liquid_1992['photo filename']   = photo_filename_wu_liquid_1992
df_wu_liquid_1992['exposure time']    = np.nan
df_wu_liquid_1992['flash']            = np.nan # TODO
df_wu_liquid_1992['x_low']            = np.nan # TODO
df_wu_liquid_1992['x_mid']            = np.nan
df_wu_liquid_1992['x_high']           = np.nan
df_wu_liquid_1992['photo page fig']   = np.nan # TODO
df_wu_liquid_1992['spectrum']         = np.nan
df_wu_liquid_1992['lighting']         = np.nan
df_wu_liquid_1992['background color'] = np.nan
df_wu_liquid_1992['grid']             = np.nan # TODO
df_wu_liquid_1992['distance']         = np.nan
df_wu_liquid_1992['f-stop']           = np.nan
df_wu_liquid_1992['focal length']     = np.nan
df_wu_liquid_1992['camera model']     = np.nan
df_wu_liquid_1992['sensor']           = np.nan

df_wu_liquid_1992 = df_wu_liquid_1992[cols]
summary_table(df_wu_liquid_1992)
df_jet_breakup = pd.concat([df_jet_breakup, df_wu_liquid_1992])

#######################
# mansour_effect_1994 #
#######################

# TODO: Add droplet size and turbulence intensity measurements.
# DONE: Add photos.
# mansour_effect_1994 pp. 132-133: > The nozzles were constructed from stainless steel hypodermic tubing, square cut on both ends and burr free. The injection orifices were examined under a 15x microscope and carefully deburred with a jeweler's broach, with particular care taken to avoid beveling corners. Tor the large capillary tube setup, the liquid entered a large cylinderical chamber followed by a 22:1 contraction. This was followed by an injector tube, D_1 = 3.051 mm diameter with a length to diameter ratio L_1/D_1 = 36. For the small capillary setup, the liquid entered a cylindrical chamber followed by a 4:1 contraction. This was followed by an injector tube, D_2 = 1.194 mm diamer with a length to diameter ratio L_2/D_2 = 255.

mansour_effect_1994_csv = pd.read_csv('../data/mansour_effects_1994/mansour_effects_1994_fig_7_2.csv', sep=',', header=0)

d_0_mansour_effect_1994          = mansour_effect_1994_csv['d_0 (mm)'] * 1e-3 # m
Ubar_0_mansour_effect_1994       = mansour_effect_1994_csv['Ubar_0 (m/s)']
L_bs_mansour_effect_1994         = mansour_effect_1994_csv['x_b/d_0']
regime_L_b_mansour_effect_1994   = mansour_effect_1994_csv['regime x_b']
regime_turb_mansour_effect_1994  = mansour_effect_1994_csv['regime turb']
regime_photo_mansour_effect_1994 = mansour_effect_1994_csv['regime photo']

# TODO: Find out if the liquid and gas information was given in the journal or conference papers.
rho_l = rho_water(T_std)
sigma = sigma_water(T_std)
nu_l  = nu_water(T_std)
P_v   = P_v_water(T_std)
rho_g = rho_ideal_gas(P_atm, T_std, MW_air)
nu_g  = mu_ideal_gas(T_std, mu_0_air, T_0_air, C_air) / rho_g

Re_l0_mansour_effect_1994      = Ubar_0_mansour_effect_1994 * d_0_mansour_effect_1994 / nu_l
We_l0_mansour_effect_1994      = rho_l * Ubar_0_mansour_effect_1994**2 * d_0_mansour_effect_1994 / sigma
I_0_mansour_effect_1994        = I_fully_developed_array(Re_l0_mansour_effect_1994, Re_turb) # TODO: Use the measured values.
Lambda_r_s_mansour_effect_1994 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_mansour_effect_1994, Re_turb), 'v')
Fr_0_mansour_effect_1994       = Ubar_0_mansour_effect_1994**2 / (g * d_0_mansour_effect_1994)
rho_s_mansour_effect_1994      = rho_l / rho_g
nu_s_mansour_effect_1994       = nu_l / nu_g
#K_c_mansour_effect_1994        = (P_atm - P_v) / (0.5 * rho_l * Ubar_0**2)
Ma_g_mansour_effect_1994       = Ubar_0_mansour_effect_1994 / c_ideal_gas(gamma_air(T_atm), T_atm + T_zero, MW_air)

i = 0
#regime_turb_mansour_effect_1994 = []
L_0s_mansour_effect_1994        = []
for Re_l0 in Re_l0_mansour_effect_1994:
   d_0 = d_0_mansour_effect_1994[i]
   
   if abs(d_0 - 3.051e-3) < 1e-8:
      L_0s_mansour_effect_1994.append(36)
   elif abs(d_0 - 1.194e-3) < 1e-8:
      L_0s_mansour_effect_1994.append(255)
   else:
      raise ValueError('Invalid nozzle diameter for mansour_effects_1994:', d_0)
   
   #if Re_l0 > Re_turb:
      #regime_turb_mansour_effect_1994.append('turbulent')
   #elif Re_l0 > Re_trans:
      #regime_turb_mansour_effect_1994.append('transitional')
   #else:
      #regime_turb_mansour_effect_1994.append('laminar')
   
   i = i + 1

K_c_mansour_effect_1994 = K_c(P_atm, P_v, rho_l, Ubar_0, K_L_guess, friction_factor_smooth_array(Re_l0_mansour_effect_1994), L_0s_mansour_effect_1994)

df_mansour_effect_1994 = pd.DataFrame({'We_l0': We_l0_mansour_effect_1994})
df_mansour_effect_1994['key']        = 'mansour_effect_1994'
df_mansour_effect_1994['alt key']    = 'mansour_effects_1994'
df_mansour_effect_1994['liquid']     = 'water' # TODO: Is this correct? I can't find it written anywhere in a quick look.
df_mansour_effect_1994['gas']        = 'air'
df_mansour_effect_1994['degas']      = np.nan
df_mansour_effect_1994['regulation'] = 'pressure regulator' # mansour_effect_1994 p. 50 fig. 4-2: pressure regulator

df_mansour_effect_1994['bend']              = np.nan
df_mansour_effect_1994['flow straightener'] = np.nan
df_mansour_effect_1994['contraction shape'] = np.nan
df_mansour_effect_1994['screen']            = np.nan
df_mansour_effect_1994['c']                 = np.nan
df_mansour_effect_1994['trip']              = np.nan
df_mansour_effect_1994['L_r/d_0']           = 0
df_mansour_effect_1994['eps_r/d_0']         = np.nan

df_mansour_effect_1994['L_0/d_0']       = L_0s_mansour_effect_1994
df_mansour_effect_1994['roughness']     = 'smooth'
df_mansour_effect_1994['f page fig']    = np.nan # TODO
df_mansour_effect_1994['pipe material'] = 'stainless steel'
df_mansour_effect_1994['est f']         = True # TODO: Add measured turbulence intensities.
df_mansour_effect_1994['check FD']      = False

df_mansour_effect_1994['t/d_0']       = np.nan
df_mansour_effect_1994['L_tip/d_0']   = np.nan
df_mansour_effect_1994['end checked'] = np.nan

df_mansour_effect_1994['orientation']         = 'down' # mansour_effect_1994 p. 51 fig. 4-2: suggests the jet is directed down
df_mansour_effect_1994['vibration isolation'] = True # Somewhere it says the jet is attached to a heavy weight.
df_mansour_effect_1994['d_chamber/d_0']       = np.inf

df_mansour_effect_1994['regime photo'] = regime_photo_mansour_effect_1994
df_mansour_effect_1994['regime L_b']   = regime_L_b_mansour_effect_1994
df_mansour_effect_1994['regime turb']  = regime_turb_mansour_effect_1994
df_mansour_effect_1994['est turb']     = True

#df_mansour_effect_1994['We_l0']      = We_l0_mansour_effect_1994
df_mansour_effect_1994['Re_l0']      = Re_l0_mansour_effect_1994
df_mansour_effect_1994['I_0']        = I_0_mansour_effect_1994
df_mansour_effect_1994['Fr_0']       = Fr_0_mansour_effect_1994
df_mansour_effect_1994['rho_s']      = rho_s_mansour_effect_1994
df_mansour_effect_1994['nu_s']       = nu_s_mansour_effect_1994
df_mansour_effect_1994['K_c']        = K_c_mansour_effect_1994
df_mansour_effect_1994['Ma_g']       = Ma_g_mansour_effect_1994
df_mansour_effect_1994['Lambda_r_s'] = Lambda_r_s_mansour_effect_1994
df_mansour_effect_1994['est rho_s']  = True
df_mansour_effect_1994['est nu_s']   = True

df_mansour_effect_1994['L_b/d_0']                      = L_bs_mansour_effect_1994
df_mansour_effect_1994['L_b/d_0 stat error']           = scipy.stats.t.ppf(1 - (1 - interval_probability_level) / 2, 20 - 1) * breakup_length_sigmas * L_bs_mansour_effect_1994 / sqrt(20) # mansour_effect_1994 p. 134: "All breakup length data reported an average over 20 frames or more." The standard deviation of the breakup length is estimated from phinney_breakup_1973.
df_mansour_effect_1994['L_b/d_0 resolution']           = 4.5 # p. 594: > The breakup length was very reproducible within one wavelength distance.
df_mansour_effect_1994['L_b method']                   = 'photographic'
df_mansour_effect_1994['L_b/d_0 page fig']             = 'p. 143, fig. 7-2'
df_mansour_effect_1994['L_b/d_0 transcription method'] = 'figure'

df_mansour_effect_1994['theta']                      = np.nan
df_mansour_effect_1994['theta stat error']           = np.nan
df_mansour_effect_1994['theta resolution']           = np.nan
df_mansour_effect_1994['theta page fig']             = np.nan
df_mansour_effect_1994['theta transcription method'] = np.nan

df_mansour_effect_1994['D_10/d_0']                   = np.nan
df_mansour_effect_1994['D_10/d_0 stat error']        = np.nan
df_mansour_effect_1994['D_10/d_0 resolution']        = np.nan
df_mansour_effect_1994['D_30/d_0']                   = np.nan
df_mansour_effect_1994['D_30/d_0 stat error']        = np.nan
df_mansour_effect_1994['D_30/d_0 resolution']        = np.nan
df_mansour_effect_1994['D_32/d_0']                   = np.nan
df_mansour_effect_1994['D_32/d_0 stat error']        = np.nan
df_mansour_effect_1994['D_32/d_0 resolution']        = np.nan
df_mansour_effect_1994['droplet x/d_0']              = np.nan
df_mansour_effect_1994['D/d_0 page fig']             = np.nan
df_mansour_effect_1994['D/d_0 transcription method'] = np.nan
df_mansour_effect_1994['D/d_0 measurement method']   = np.nan
df_mansour_effect_1994['v_d_bar/vp']                        = np.nan
df_mansour_effect_1994['v_d_bar/vp stat error']             = np.nan
df_mansour_effect_1994['v_d_bar/vp resolution']             = np.nan
df_mansour_effect_1994['v_d_bar/vp page fig']               = np.nan
df_mansour_effect_1994['v_d_bar/vp transcription method']   = np.nan
df_mansour_effect_1994['e_p']                        = np.nan
df_mansour_effect_1994['e_p page fig']               = np.nan
df_mansour_effect_1994['e_p transcription method']   = np.nan

df_mansour_effect_1994['x_trans/d_0']                      = np.nan
df_mansour_effect_1994['x_trans/d_0 page fig']             = np.nan
df_mansour_effect_1994['x_trans/d_0 transcription method'] = np.nan

df_mansour_effect_1994['x_i/d_0']                      = np.nan
df_mansour_effect_1994['x_i/d_0 stat error']           = np.nan
df_mansour_effect_1994['x_i/d_0 resolution']           = np.nan
df_mansour_effect_1994['x_i/d_0 page fig']             = np.nan
df_mansour_effect_1994['x_i/d_0 transcription method'] = np.nan

df_mansour_effect_1994['x_e/d_0']                      = np.nan
df_mansour_effect_1994['x_e/d_0 stat error']           = np.nan
df_mansour_effect_1994['x_e/d_0 resolution']           = np.nan
df_mansour_effect_1994['x_e/d_0 page fig']             = np.nan
df_mansour_effect_1994['x_e/d_0 transcription method'] = np.nan

df_mansour_effect_1994['photo filename']   = np.nan
df_mansour_effect_1994['exposure time']    = np.nan
df_mansour_effect_1994['flash']            = np.nan
df_mansour_effect_1994['x_low']            = np.nan
df_mansour_effect_1994['x_mid']            = np.nan
df_mansour_effect_1994['x_high']           = np.nan
df_mansour_effect_1994['photo page fig']   = np.nan
df_mansour_effect_1994['spectrum']         = np.nan
df_mansour_effect_1994['lighting']         = np.nan
df_mansour_effect_1994['background color'] = np.nan
df_mansour_effect_1994['grid']             = np.nan
df_mansour_effect_1994['distance']         = np.nan
df_mansour_effect_1994['f-stop']           = np.nan
df_mansour_effect_1994['focal length']     = np.nan
df_mansour_effect_1994['camera model']     = np.nan
df_mansour_effect_1994['sensor']           = np.nan

df_mansour_effect_1994 = df_mansour_effect_1994[cols]
summary_table(df_mansour_effect_1994)
df_jet_breakup = pd.concat([df_jet_breakup, df_mansour_effect_1994])

##########################
# sallam_properties_2002 #
##########################

# TODO: Add.
# theta (low pressure atmosphere), D_{32}, L_b, ligaments

# Issues with this data:
# 1. The Reynolds numbers I am computing differ from those reported in the paper and dissertation.

# pdf p. 165: data table

# sallam_properties_2002 p. 108 (pdf p. 128): breakup length standard deviation

# p. 26, tab. 2.3: test conditions table
# p. 26: still air at $99\pm0.5$ kPa and $297\pm0.5$ K ($\rho_\text{g} - 1.16$ kg/m$^3$ and $\nu_\text{g} = 15.9$ mm$^2$/s). Round injector with a rounded entry and a length-to-diameter ratio of 40:1.
T_sallam_properties_2002 = 297 - T_zero # C

sallam_properties_2002_csv = pd.read_csv('../data/sallam_properties_2002/sallam_properties_2002.csv', sep=',', header=0)

theta_page_fig_sallam_properties_2002 = sallam_properties_2002_csv['theta page fig']
x_i_page_fig_sallam_properties_2002   = sallam_properties_2002_csv['x_i page fig']
x_b_page_fig_sallam_properties_2002   = sallam_properties_2002_csv['xbavg page fig']
liquid_sallam_properties_2002         = sallam_properties_2002_csv['liquid']
d_0_sallam_properties_2002            = sallam_properties_2002_csv['d (mm)'] * 1e-3 # m
We_l0_sallam_properties_2002          = sallam_properties_2002_csv['We_l0']
theta_sallam_properties_2002          = sallam_properties_2002_csv['theta (deg)'] * (np.pi / 180)
x_is_sallam_properties_2002           = sallam_properties_2002_csv['x_i/d_0']
x_bs_sallam_properties_2002           = sallam_properties_2002_csv['xbavg/d_0']
regime_L_b_sallam_properties_2002     = sallam_properties_2002_csv['regime L_b']

nu_g  = 15.9e-6 # m^2/s
rho_g = 1.16    # kg/m^3
L_0s  = 40.

i = 0
Re_l0_sallam_properties_2002 = np.zeros(len(We_l0_sallam_properties_2002))
Fr_0_sallam_properties_2002  = np.zeros(len(We_l0_sallam_properties_2002))
rho_s_sallam_properties_2002 = np.zeros(len(We_l0_sallam_properties_2002))
nu_s_sallam_properties_2002  = np.zeros(len(We_l0_sallam_properties_2002))
K_c_sallam_properties_2002   = np.zeros(len(We_l0_sallam_properties_2002))
Ma_g_sallam_properties_2002  = np.zeros(len(We_l0_sallam_properties_2002))
regime_turb_sallam_properties_2002    = []
for We_l0, d_0, liquid in zip(We_l0_sallam_properties_2002, d_0_sallam_properties_2002, liquid_sallam_properties_2002):
   
   # p. 26 (pdf p. 46)
   if liquid == 'water':
      mu_l  = 8.94e-4 # kg/(m*s)
      sigma = 70.8e-3 # N/m
      rho_l = 997.    # kg/m^3
   elif liquid == 'ethanol':
      mu_l  = 16.0e-4 # kg/(m*s)
      sigma = 24.0e-3 # N/m
      rho_l = 800.    # kg/m^3
   else:
      sys.exit('Invalid liquid for sallam_properties_2002.')
   
   nu_l  = mu_l / rho_l
   
   Ubar_0 = sqrt(We_l0 * sigma / (rho_l * d_0))
   
   Re_l0_sallam_properties_2002[i] = Ubar_0 * d_0 / nu_l
   Fr_0_sallam_properties_2002[i]  = Ubar_0**2 / (g * d_0)
   rho_s_sallam_properties_2002[i] = rho_l / rho_g
   nu_s_sallam_properties_2002[i]  = nu_l / nu_g
   #K_c_sallam_properties_2002[i]   = (P_g - P_v_water(T_sallam_properties_2002)) / (0.5 * rho_l * Ubar_0**2)
   K_c_sallam_properties_2002[i]   = K_c(P_g, P_v_water(T_sallam_properties_2002), rho_l, Ubar_0, K_L_wellrounded, friction_factor_smooth(Re_l0_sallam_properties_2002[i]), L_0s)
   Ma_g_sallam_properties_2002[i]  = Ubar_0 / c_ideal_gas(gamma_air(T_sallam_properties_2002), T_sallam_properties_2002 + T_zero, MW_air)
   
   if Re_l0 > Re_turb:
      regime_turb_sallam_properties_2002.append('turbulent')
   elif Re_l0 > Re_trans:
      regime_turb_sallam_properties_2002.append('transitional')
   else:
      regime_turb_sallam_properties_2002.append('laminar')
   
   i = i + 1

I_0_sallam_properties_2002 = I_fully_developed_array(Re_l0_sallam_properties_2002, Re_turb)
Lambda_r_s_sallam_properties_2002 = Lambda_s_array('smooth', friction_factor_smooth_array(Re_l0_sallam_properties_2002, Re_turb), 'v')

df_sallam_properties_2002 = pd.DataFrame({'We_l0': We_l0_sallam_properties_2002})
df_sallam_properties_2002['key']        = 'sallam_properties_2002'
df_sallam_properties_2002['alt key']    = 'sallam_liquid_2002'
df_sallam_properties_2002['liquid']     = 'water'
df_sallam_properties_2002['gas']        = 'air'
df_sallam_properties_2002['degas']      = np.nan
df_sallam_properties_2002['regulation'] = np.nan

df_sallam_properties_2002['bend']              = False
df_sallam_properties_2002['flow straightener'] = False
df_sallam_properties_2002['contraction shape'] = 'smooth'
df_sallam_properties_2002['screen']            = False
df_sallam_properties_2002['c']                 = 165.e-3 / d_0_sallam_properties_2002 # p. 16
df_sallam_properties_2002['trip']              = False
df_sallam_properties_2002['L_r/d_0']           = 0
df_sallam_properties_2002['eps_r/d_0']         = np.nan

df_sallam_properties_2002['L_0/d_0']       = L_0s
df_sallam_properties_2002['roughness']     = 'smooth'
df_sallam_properties_2002['f page fig']    = np.nan
df_sallam_properties_2002['pipe material'] = np.nan
df_sallam_properties_2002['est f']         = True # TODO: Use TI data given.
df_sallam_properties_2002['check FD']      = np.nan

df_sallam_properties_2002['t/d_0']       = np.nan
df_sallam_properties_2002['L_tip/d_0']   = np.nan
df_sallam_properties_2002['end checked'] = np.nan

df_sallam_properties_2002['orientation']         = 'down'
df_sallam_properties_2002['vibration isolation'] = np.nan # Unclear.
df_sallam_properties_2002['d_chamber/d_0']       = np.inf

df_sallam_properties_2002['regime photo'] = np.nan
df_sallam_properties_2002['regime L_b']   = regime_L_b_sallam_properties_2002
df_sallam_properties_2002['regime turb']  = regime_turb_sallam_properties_2002
df_sallam_properties_2002['est turb']     = True

#df_sallam_properties_2002['We_l0']      = We_l0_sallam_properties_2002
df_sallam_properties_2002['Re_l0']      = Re_l0_sallam_properties_2002
df_sallam_properties_2002['I_0']        = I_0_sallam_properties_2002
df_sallam_properties_2002['Fr_0']       = Fr_0_sallam_properties_2002
df_sallam_properties_2002['rho_s']      = rho_s_sallam_properties_2002
df_sallam_properties_2002['nu_s']       = nu_s_sallam_properties_2002
df_sallam_properties_2002['K_c']        = K_c_sallam_properties_2002
df_sallam_properties_2002['Ma_g']       = Ma_g_sallam_properties_2002
df_sallam_properties_2002['Lambda_r_s'] = Lambda_r_s_sallam_properties_2002
df_sallam_properties_2002['est rho_s']  = False
df_sallam_properties_2002['est nu_s']   = False

df_sallam_properties_2002['L_b/d_0']                      = x_bs_sallam_properties_2002
df_sallam_properties_2002['L_b/d_0 stat error']           = 0.1 * x_bs_sallam_properties_2002 # sallam_liquid_2002 p. 430: > Several pulsed shadowgraph images were averaged in order to find mean liquid column breakup lengths with experimental uncertainties (95% confidence) less than 10%.
df_sallam_properties_2002['L_b/d_0 resolution']           = 0.5e-3 / d_0_sallam_properties_2002 # sallam_liquid_2002 p. 430: > The nozzle assembly could be traversed vertically with an accuracy of 0.5 mm
df_sallam_properties_2002['L_b method']                   = 'photographic'
df_sallam_properties_2002['L_b/d_0 page fig']             = x_b_page_fig_sallam_properties_2002
df_sallam_properties_2002['L_b/d_0 transcription method'] = 'table'

df_sallam_properties_2002['theta']                      = theta_sallam_properties_2002
df_sallam_properties_2002['theta stat error']           = np.nan
df_sallam_properties_2002['theta resolution']           = np.nan
df_sallam_properties_2002['theta page fig']             = np.nan
df_sallam_properties_2002['theta transcription method'] = 'table'

df_sallam_properties_2002['D_10/d_0']                   = np.nan
df_sallam_properties_2002['D_10/d_0 stat error']        = np.nan
df_sallam_properties_2002['D_10/d_0 resolution']        = np.nan
df_sallam_properties_2002['D_30/d_0']                   = np.nan
df_sallam_properties_2002['D_30/d_0 stat error']        = np.nan
df_sallam_properties_2002['D_30/d_0 resolution']        = np.nan
df_sallam_properties_2002['D_32/d_0']                   = np.nan
df_sallam_properties_2002['D_32/d_0 stat error']        = np.nan
df_sallam_properties_2002['D_32/d_0 resolution']        = np.nan
df_sallam_properties_2002['droplet x/d_0']              = np.nan
df_sallam_properties_2002['D/d_0 page fig']             = np.nan
df_sallam_properties_2002['D/d_0 transcription method'] = np.nan
df_sallam_properties_2002['D/d_0 measurement method']   = np.nan
df_sallam_properties_2002['v_d_bar/vp']                        = np.nan
df_sallam_properties_2002['v_d_bar/vp stat error']             = np.nan
df_sallam_properties_2002['v_d_bar/vp resolution']             = np.nan
df_sallam_properties_2002['v_d_bar/vp page fig']               = np.nan
df_sallam_properties_2002['v_d_bar/vp transcription method']   = np.nan
df_sallam_properties_2002['e_p']                        = np.nan
df_sallam_properties_2002['e_p page fig']               = np.nan
df_sallam_properties_2002['e_p transcription method']   = np.nan

df_sallam_properties_2002['x_trans/d_0']                      = np.nan
df_sallam_properties_2002['x_trans/d_0 page fig']             = np.nan
df_sallam_properties_2002['x_trans/d_0 transcription method'] = np.nan

df_sallam_properties_2002['x_i/d_0']                      = x_is_sallam_properties_2002
df_sallam_properties_2002['x_i/d_0 stat error']           = 0.33 * x_is_sallam_properties_2002 # Assumed same as you did for wu_liquid_1992.
df_sallam_properties_2002['x_i/d_0 resolution']           = x_is_sallam_properties_2002 * 0.023 # p. 125 gives minimum 2.3%
df_sallam_properties_2002['x_i/d_0 page fig']             = x_i_page_fig_sallam_properties_2002
df_sallam_properties_2002['x_i/d_0 transcription method'] = 'table'

df_sallam_properties_2002['x_e/d_0']                      = np.nan
df_sallam_properties_2002['x_e/d_0 stat error']           = np.nan
df_sallam_properties_2002['x_e/d_0 resolution']           = np.nan
df_sallam_properties_2002['x_e/d_0 page fig']             = np.nan
df_sallam_properties_2002['x_e/d_0 transcription method'] = np.nan

# TODO: Add photos.
df_sallam_properties_2002['photo filename']   = np.nan
df_sallam_properties_2002['exposure time']    = np.nan
df_sallam_properties_2002['flash']            = np.nan
df_sallam_properties_2002['x_low']            = np.nan
df_sallam_properties_2002['x_mid']            = np.nan
df_sallam_properties_2002['x_high']           = np.nan
df_sallam_properties_2002['photo page fig']   = np.nan
df_sallam_properties_2002['spectrum']         = np.nan
df_sallam_properties_2002['lighting']         = np.nan
df_sallam_properties_2002['background color'] = np.nan
df_sallam_properties_2002['grid']             = np.nan
df_sallam_properties_2002['distance']         = np.nan
df_sallam_properties_2002['f-stop']           = np.nan
df_sallam_properties_2002['focal length']     = np.nan
df_sallam_properties_2002['camera model']     = np.nan
df_sallam_properties_2002['sensor']           = np.nan

df_sallam_properties_2002 = df_sallam_properties_2002[cols]
summary_table(df_sallam_properties_2002)
df_jet_breakup = pd.concat([df_jet_breakup, df_sallam_properties_2002])

##########################
# mayer_atomization_2004 #
##########################

# spray angle: p. 538, fig. 31
# TODO: Add.

######################################
# osta_effect_2010 / osta_study_2012 #
######################################

# TODO: Add.
# theta, ligaments
# theta table on p. 52 (pdf p. 70)
# see table B.1 on p. 109 (pdf p. 127) for velocity estimate function

#########################
# trettel_modeling_2020 #
#########################

# TODO: Add t/d_0.

trettel_modeling_2020_csv = pd.read_csv('../../../experiments/consolidated-trajectory-experiments.csv', sep=',', header=0)

frame_rate = 29.97 # frames/s

photo_filename_trettel_modeling_2020  = trettel_modeling_2020_csv['video filename'] # C
temp_trettel_modeling_2020_weather    = (5./9.) * (trettel_modeling_2020_csv['temp. (F, weather service)'] - 32.0) # C
temp_trettel_modeling_2020_anemometer = trettel_modeling_2020_csv['temp. (C, anemometer)']
P_atm_trettel_modeling_2020_temp      = trettel_modeling_2020_csv['P_atm (inHg)'] * 3386.39 # Pa
d_0_trettel_modeling_2020             = in_to_m(trettel_modeling_2020_csv['d_0 (in)']) # m
L_0_trettel_modeling_2020             = in_to_m(trettel_modeling_2020_csv['L_0 (in)']) # m
roughness_trettel_modeling_2020       = trettel_modeling_2020_csv['roughness']
Vol_trettel_modeling_2020             = trettel_modeling_2020_csv['Vol. (L)'] * 0.001 # m^3
T_trettel_modeling_2020               = (trettel_modeling_2020_csv['end frame']- trettel_modeling_2020_csv['start frame']) / frame_rate # s
regime_photo_trettel_modeling_2020    = trettel_modeling_2020_csv['regime photo']

temp_trettel_modeling_2020 = np.array([])
for temp_weather, temp_anemometer in zip(temp_trettel_modeling_2020_weather, temp_trettel_modeling_2020_anemometer):
   if temp_weather > 0.: # not sure why not(temp_weather is np.nan) doesn't work
      temp_trettel_modeling_2020 = np.append(temp_trettel_modeling_2020, temp_weather)
   else:
      temp_trettel_modeling_2020 = np.append(temp_trettel_modeling_2020, temp_anemometer)

P_atm_trettel_modeling_2020 = np.array([])
for P_atm_i in P_atm_trettel_modeling_2020_temp:
   if P_atm_i > 0.:
      P_atm_trettel_modeling_2020 = np.append(P_atm_trettel_modeling_2020, P_atm_i)
   else:
      P_atm_trettel_modeling_2020 = np.append(P_atm_trettel_modeling_2020, P_atm)

Q_trettel_modeling_2020      = Vol_trettel_modeling_2020 / T_trettel_modeling_2020
A_0_trettel_modeling_2020    = (np.pi / 4.) * d_0_trettel_modeling_2020**2.
Ubar_0_trettel_modeling_2020 = Q_trettel_modeling_2020 / A_0_trettel_modeling_2020

rho_g_trettel_modeling_2020 = rho_ideal_gas(P_atm_trettel_modeling_2020, temp_trettel_modeling_2020, MW_air)
nu_g_trettel_modeling_2020  = mu_ideal_gas(temp_trettel_modeling_2020, mu_0_air, T_0_air, C_air) / rho_g_trettel_modeling_2020

rho_l_trettel_modeling_2020 = rho_water(temp_trettel_modeling_2020)
sigma_trettel_modeling_2020 = sigma_water(temp_trettel_modeling_2020)
nu_l_trettel_modeling_2020  = nu_water(temp_trettel_modeling_2020)
P_v_trettel_modeling_2020   = P_v_water(temp_trettel_modeling_2020)

We_l0_trettel_modeling_2020      = rho_l_trettel_modeling_2020 * Ubar_0_trettel_modeling_2020**2. * d_0_trettel_modeling_2020 / sigma_trettel_modeling_2020
Re_l0_trettel_modeling_2020      = Ubar_0_trettel_modeling_2020 * d_0_trettel_modeling_2020 / nu_l_trettel_modeling_2020
f_trettel_modeling_2020 = friction_factor_smooth_array(Re_l0_trettel_modeling_2020)
Fr_0_trettel_modeling_2020       = Ubar_0_trettel_modeling_2020**2. / (d_0_trettel_modeling_2020 * g)
rho_s_trettel_modeling_2020      = rho_l_trettel_modeling_2020 / rho_g_trettel_modeling_2020
nu_s_trettel_modeling_2020       = nu_l_trettel_modeling_2020 / nu_g_trettel_modeling_2020
K_c_trettel_modeling_2020        = K_c(P_atm_trettel_modeling_2020, P_v_trettel_modeling_2020, rho_l_trettel_modeling_2020, Ubar_0_trettel_modeling_2020, K_L_sharpedged, friction_factor_smooth_array(Re_l0_trettel_modeling_2020), L_0_trettel_modeling_2020/d_0_trettel_modeling_2020)

Ma_g_trettel_modeling_2020        = []
Lambda_r_s_trettel_modeling_2020  = []
regime_turb_trettel_modeling_2020 = []
I_0_trettel_modeling_2020         = []
for Ubar_0, temp, roughness, f, Re_l0 in zip(Ubar_0_trettel_modeling_2020, temp_trettel_modeling_2020, roughness_trettel_modeling_2020, f_trettel_modeling_2020, Re_l0_trettel_modeling_2020):
   Ma_g_trettel_modeling_2020.append(Ubar_0 / c_ideal_gas(gamma_air(temp), temp + T_zero, MW_air))
   
   if Re_l0 > Re_turb:
      regime_turb_trettel_modeling_2020.append('turbulent')
      Lambda_r_s_trettel_modeling_2020.append(Lambda_s(roughness, f, 'v'))
      I_0_trettel_modeling_2020.append(I_fully_developed(friction_factor_smooth(Re_l0)))
   elif Re_l0 > Re_trans:
      regime_turb_trettel_modeling_2020.append('transitional')
      Lambda_r_s_trettel_modeling_2020.append(np.nan)
      I_0_trettel_modeling_2020.append(np.nan)
   else:
      regime_turb_trettel_modeling_2020.append('laminar')
      Lambda_r_s_trettel_modeling_2020.append(np.nan)
      I_0_trettel_modeling_2020.append(np.nan)

df_trettel_modeling_2020 = pd.DataFrame({'key': 'trettel_modeling_2020',
                                       'We_l0': We_l0_trettel_modeling_2020})
df_trettel_modeling_2020['liquid']     = 'water'
df_trettel_modeling_2020['alt key']    = np.nan
df_trettel_modeling_2020['gas']        = 'air'
df_trettel_modeling_2020['degas']      = False
df_trettel_modeling_2020['regulation'] = True

df_trettel_modeling_2020['bend']              = True
df_trettel_modeling_2020['flow straightener'] = True
df_trettel_modeling_2020['screen']            = True
df_trettel_modeling_2020['contraction shape'] = 'rough'
df_trettel_modeling_2020['c']                 = 0.824 * 2.54e-2 / d_0_trettel_modeling_2020
df_trettel_modeling_2020['trip']              = False
df_trettel_modeling_2020['L_r/d_0']           = 0
df_trettel_modeling_2020['eps_r/d_0']         = np.nan

df_trettel_modeling_2020['L_0/d_0']       = L_0_trettel_modeling_2020 / d_0_trettel_modeling_2020
df_trettel_modeling_2020['roughness']     = roughness_trettel_modeling_2020
df_trettel_modeling_2020['f page fig']    = np.nan
df_trettel_modeling_2020['pipe material'] = 'brass'
df_trettel_modeling_2020['est f']         = True
df_trettel_modeling_2020['check FD']      = False

df_trettel_modeling_2020['t/d_0']       = np.nan
df_trettel_modeling_2020['L_tip/d_0']   = np.nan
df_trettel_modeling_2020['end checked'] = np.nan

df_trettel_modeling_2020['orientation']         = 'angled'
df_trettel_modeling_2020['vibration isolation'] = False
df_trettel_modeling_2020['d_chamber/d_0']       = np.nan

#df_trettel_modeling_2020['We_l0']      = We_l0_trettel_modeling_2020
df_trettel_modeling_2020['Re_l0']      = Re_l0_trettel_modeling_2020
df_trettel_modeling_2020['I_0']        = I_0_trettel_modeling_2020
df_trettel_modeling_2020['Fr_0']       = Fr_0_trettel_modeling_2020
df_trettel_modeling_2020['rho_s']      = rho_s_trettel_modeling_2020
df_trettel_modeling_2020['nu_s']       = nu_s_trettel_modeling_2020
df_trettel_modeling_2020['K_c']        = K_c_trettel_modeling_2020
df_trettel_modeling_2020['Ma_g']       = Ma_g_trettel_modeling_2020
df_trettel_modeling_2020['Lambda_r_s'] = Lambda_r_s_trettel_modeling_2020
df_trettel_modeling_2020['est rho_s']  = False
df_trettel_modeling_2020['est nu_s']   = True

df_trettel_modeling_2020['regime photo'] = regime_photo_trettel_modeling_2020
df_trettel_modeling_2020['regime L_b']   = np.nan
df_trettel_modeling_2020['regime turb']  = regime_turb_trettel_modeling_2020
df_trettel_modeling_2020['est turb']     = True

df_trettel_modeling_2020['L_b/d_0']                      = np.nan
df_trettel_modeling_2020['L_b/d_0 stat error']           = np.nan
df_trettel_modeling_2020['L_b/d_0 resolution']           = np.nan
df_trettel_modeling_2020['L_b method']                   = np.nan
df_trettel_modeling_2020['L_b/d_0 page fig']             = np.nan
df_trettel_modeling_2020['L_b/d_0 transcription method'] = np.nan

df_trettel_modeling_2020['theta']                      = np.nan
df_trettel_modeling_2020['theta stat error']           = np.nan
df_trettel_modeling_2020['theta resolution']           = np.nan
df_trettel_modeling_2020['theta page fig']             = np.nan
df_trettel_modeling_2020['theta transcription method'] = np.nan

df_trettel_modeling_2020['D_10/d_0']                        = np.nan
df_trettel_modeling_2020['D_10/d_0 stat error']             = np.nan
df_trettel_modeling_2020['D_10/d_0 resolution']             = np.nan
df_trettel_modeling_2020['D_32/d_0']                        = np.nan
df_trettel_modeling_2020['D_32/d_0 stat error']             = np.nan
df_trettel_modeling_2020['D_32/d_0 resolution']             = np.nan
df_trettel_modeling_2020['D_30/d_0']                        = np.nan
df_trettel_modeling_2020['D_30/d_0 stat error']             = np.nan
df_trettel_modeling_2020['D_30/d_0 resolution']             = np.nan
df_trettel_modeling_2020['droplet x/d_0']                   = np.nan
df_trettel_modeling_2020['D/d_0 page fig']                  = np.nan
df_trettel_modeling_2020['D/d_0 transcription method']      = np.nan
df_trettel_modeling_2020['D/d_0 measurement method']        = np.nan
df_trettel_modeling_2020['v_d_bar/vp']                      = np.nan
df_trettel_modeling_2020['v_d_bar/vp stat error']           = np.nan
df_trettel_modeling_2020['v_d_bar/vp resolution']           = np.nan
df_trettel_modeling_2020['v_d_bar/vp page fig']             = np.nan
df_trettel_modeling_2020['v_d_bar/vp transcription method'] = np.nan
df_trettel_modeling_2020['e_p']                             = np.nan
df_trettel_modeling_2020['e_p page fig']                    = np.nan
df_trettel_modeling_2020['e_p transcription method']        = np.nan

df_trettel_modeling_2020['x_trans/d_0']                      = np.nan
df_trettel_modeling_2020['x_trans/d_0 page fig']             = np.nan
df_trettel_modeling_2020['x_trans/d_0 transcription method'] = np.nan

df_trettel_modeling_2020['x_i/d_0']                      = np.nan
df_trettel_modeling_2020['x_i/d_0 stat error']           = np.nan
df_trettel_modeling_2020['x_i/d_0 resolution']           = np.nan
df_trettel_modeling_2020['x_i/d_0 page fig']             = np.nan
df_trettel_modeling_2020['x_i/d_0 transcription method'] = np.nan

df_trettel_modeling_2020['x_e/d_0']                      = np.nan
df_trettel_modeling_2020['x_e/d_0 stat error']           = np.nan
df_trettel_modeling_2020['x_e/d_0 resolution']           = np.nan
df_trettel_modeling_2020['x_e/d_0 page fig']             = np.nan
df_trettel_modeling_2020['x_e/d_0 transcription method'] = np.nan

df_trettel_modeling_2020['photo filename']   = photo_filename_trettel_modeling_2020
df_trettel_modeling_2020['exposure time']    = np.nan
df_trettel_modeling_2020['flash']            = False
df_trettel_modeling_2020['x_low']            = 0.
df_trettel_modeling_2020['x_mid']            = np.nan
df_trettel_modeling_2020['x_high']           = np.nan
df_trettel_modeling_2020['photo page fig']   = np.nan
df_trettel_modeling_2020['spectrum']         = 'visible'
df_trettel_modeling_2020['lighting']         = 'ambient, outdoors'
df_trettel_modeling_2020['background color'] = np.nan
df_trettel_modeling_2020['grid']             = False
df_trettel_modeling_2020['distance']         = np.nan
df_trettel_modeling_2020['f-stop']           = np.nan
df_trettel_modeling_2020['focal length']     = np.nan
df_trettel_modeling_2020['camera model']     = 'Contour ROAM model 1600'
df_trettel_modeling_2020['sensor']           = np.nan

df_trettel_modeling_2020 = df_trettel_modeling_2020[cols]
summary_table(df_trettel_modeling_2020)
df_jet_breakup = pd.concat([df_jet_breakup, df_trettel_modeling_2020])

#######
# END #
#######

summary_table(df_jet_breakup)

with open('../outputs/'+data_file+'.pickle', 'w') as f:
   pickle.dump([df_jet_breakup, metadata], f)

macros_breakdown = open('../outputs/macros/breakdown.tex', 'w')

macros_breakdown.write(r'\newcommand{\compilationcitations}{\cite{')

key_array = []
for key in df_jet_breakup['key']:
   if key == 'eisenklam_flow_1958':
      # removed due to high Re_crit
      continue
   elif key == 'chen_mechanics_1962':
      key = 'chen_disintegration_1964'
   elif key == 'skrebkov_turbulentnyye_1963':
      key = 'skrebkov_turbulent_1966'
   elif key == 'grant_newtonian_1965':
      key = 'grant_newtonian_1966'
   elif key == 'kim_investigation_1983':
      # removed due to not actually measuring \xiavg
      continue
   elif key == 'wu_atomizing_1983':
      key = 'wu_measurements_1983'
   elif key == 'wu_liquid_1992':
      key = 'wu_primary_1992'
   elif key == 'mansour_effect_1994':
      key = 'mansour_effects_1994'
   elif key == 'sallam_properties_2002':
      key = 'sallam_liquid_2002'
   
   if not(key in key_array):
      key_array.append(key)
      if len(key_array) == 1:
         macros_breakdown.write(key)
      else:
         macros_breakdown.write(','+key)

macros_breakdown.write('}}\n')

#print key_array

macros_breakdown.close()
