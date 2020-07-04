# TODO

- Update to Python 3
- Convert variable names to be consistent with papers (e.g., L_b ==> x_b, I_0 ==> Tubar_0, v_d_bar ==> vdavg)
- Change regime names to be consistent with papers and not archaic
- Add more data validation checks (see section below)
- Add type hinting (both function and variable annotation, with latest syntax). Some links:
  - https://towardsdatascience.com/down-with-technical-debt-clean-python-for-data-scientists-aa7592eff7fc
  - https://www.bernat.tech/the-state-of-type-hints-in-python/
  - https://stackoverflow.com/questions/32557920/what-are-type-hints-in-python-3-5
- Make many loops use zip rather than indexing
- Add data and script for pipe turbulence statistics
- Add dimension checks via a Python package that does dimensional analysis (via [Pint](https://github.com/hgrecco/pint)); use for unit conversions as well, e.g., make the raw CSV spreadsheets in the units as printed and convert rather than putting the raw CSV spreadsheets in SI units
- Use the [uncertainties](https://github.com/lebigot/uncertainties) package to track uncertainties.
- Require quantified uncertainty for all quantities
- Use one file for each study to make editing the code easier.
- Don't make queries directly to the database and instead use OOP or functions. This would allow you to switch between pandas and SQLite if you want to.
- "codebook" per http://datacolada.org/69
- standardized input CSV file form?
- Have Gnumeric spreadsheets for every CSV file, and use [data validation](http://www.fifi.org/doc/gnumeric-doc/html/C/protectionandvalidation.html) in the spreadsheets
- Change x_low, etc. to x_s_low as this is actually x_low/d_0.
- Regression considering uncertainty.
- Bayesian experimental design
- Bayesian model selection
- Add more citations listed in the "Data sources to add" section below to the bib file. Stopped at kitamura_influence_1978.
- https://liw.fi/hackingdoc/
- https://liw.fi/free-software-website/
- https://vincent.bernat.ch/en/blog/2019-sustainable-python-script
- Test all functions used. Start with assertions on all inputs and outputs (test bounds and dimensions), [Hypothesis for property testing](https://datascience.blog.wzb.eu/2019/11/08/property-based-testing-for-scientific-code-in-python/) and [doctest](https://en.wikipedia.org/wiki/Doctest).
- [Require keyword arguments](http://www.jefftk.com/p/require-keyword-arguments)
- Start using Pylint regularly.

# Output data format

## pandas

### Advantages

- Already used (though probably poorly)
- Integrates with Pint
- Simple syntax in most cases eases data exploration

### Disadvantages

- Python-only
- Data file (Pickle file) is also Python only
- Pickle files are not a recognized archival standard (Indeed, I have had trouble loading Pickle files after upgrading Python)
- Slow, though not at the size of the current database
- New: started in 2008
- Constraints seem irritating to add

## SQLite

### Advantages

- [A LOC recommended archival standard](https://www.loc.gov/preservation/resources/rfs/data.html)
- Fast
- Can readily be used in other programming languages
- [Wide support](https://en.wikipedia.org/wiki/Sqlite)
- Older than pandas: started in 2000
- [Constraints](https://www.tutorialspoint.com/sqlite/sqlite_constraints.htm) are easy

### Disadvantages 

- Doesn't support booleans (but 0 and 1 are okay)
- Harder run queries than pandas in most cases

# Data sources to add

- Add photos from trettel_turbulent_2020.
- Obtain more data for liquid-liquid systems, liquid jets in cross flow, liquid jets with co-flows
- Get more spray angles from photos.
- Add more data from Phinney for various atmospheric densities. (Not sure what I meant by this. I'll check if all of Phinney's data has been added.)
- Re-add data from wu_liquid_1992. Type raw data in and process to produce the Weber number, e.g., table B.5 has velocities not entered in the spreadsheet. Compare this against the one in there right now, which I believe was produced from figures? Also check that the Weber and Reynolds numbers are actually higher than 1e6 as this is the only data I have that high.
- duffie_factors_1953/duffie_factors_1951
- richardson_mechanism_1954
- keith_liquid-liquid_1955/keith_drop_1951
- christiansen_breakup_1957/christiansen_influence_1957
- Add atomization regime markers from palmer_water_1962 fig. 2 (also see report by same author).
- ryley_construction_1963
- Replace Grant photos with higher quality versions from grant_newtonian_1965-photos.pdf
- smirnov_effect_1965
- meister_formation_1966/scheele_drop_1968/scheele_drop_1968-1/meister_drop_1969/meister_prediction_1969
- rao_drop_1966 ("capillary" but what is its length?)
- dotson_study_1967
- newman_preliminary_1967
- kroesser_stability_1968, kroesser_viscoelastic_1969
- takahashi_effect_1969, takahashi_laminar_1971/kitamura_stability_1986, takahashi_stability_1972/takahashi_breakup_1972/takahashi_laminar_1971, kitamura_stability_1982
- Add fenn_newtonian_1969/fenn_ambient_1968 (Should be useful to examine density ratio effects on LR-DT boundary.)
- rutland_theoretical_1970 (see p. 1692R: 4 mm diameter, 30 cm length nozzle)
- skelland_dispersed_1971/minhas_dispersed_1969
- gordon_non-newtonian_1973 (both Newtonian and non-Newtonian fluids)
- parkin_production_1973 (p. 67: nozzle not long enough)
- skelland_jet_1974/johnson_jet_1973
- van_de_sande_air_1974 fig III.13 (pdf p. 41), fig. III.15 (turbulent Rayleigh, pdf p. 43; also see van_de_sande_jet_1976 p. 221L), fig. III.17 (can't disambiguate, pdf p. 46, also van_de_sande_surface_1973 fig. 3); van_de_sande_air_1974 pdf p. "each value plotted is the arithmetic mean of at least five photographs"
- dunn_high_1975/dunn_stability_1974
- lafrance_capillary_1975 (most experimental results are forced, but p. 78 has some earlier unforced data), lafrance_drop_1974 (droplet size and breakup length in laminar and turbulent Rayleigh regimes)
- vliem_influence_1975/sterling_mechanisms_1981 (droplet size and breakup length for turbulent Rayleigh)
- skelland_dispersed_1977
- kitamura_influence_1978
  - fig. 8: \rhol/\rhog \approx 1.3, highest Re \approx 3000
  - fig. 9: \rhol/\rhog \approx 0.8, highest Re \approx 1500
- van_den_akker_spontaneous_1980 (long enough for laminar fully developed flow?)
- abbott_analysis_1982
- eroglu_coaxial_1991/eroglu_initial_1991/eroglu_wave_1991/farago_morphological_1992/farago_parametric_1990
- Add breakup lengths from sauerwein_theoretische_1992 pp. 121--122? Not clear without translation if this is fully developed
- mayer_zur_1993/mayer_coaxial_1994/mayer_rocket_1995
- olinger_lock-states_1993 (forced dripping)
- hardalupas_characteristics_1994
- engelbert_breakup_1995
- Add Ruff and Tseng data, which apparently has low \rhol/\rhog photographic regimes according to wu_effects_1995 p. 189, fig. 7 (pdf p. 15) / faeth_structure_1995 p. 117, fig. 1.5.
- Based on the same figure as the previous line's TODO, so does wu_effects_1995.
- Add amagai_frequency_1997/arai_surface_1999
- leroux_break-up_1997 (ask Prof. Dumouchel for data, redundant with leroux_stability_1996?)
- mun_effects_1998 (L_0/d_0 = 25.6. Not long enough for laminar jets?)
- rhim_measurement_1998 (breakup lengths)
- Add tamaki_effects_1998
- Add clanet_transition_1999
- Add blaisot_determination_2000, godelle_phase_2000
- Add godelle_symbolic_2000
- Add malot_experimental_2001 (has a L_0/d_0 = 50 case; ask Prof. Dumouchel for data)
- blaisot_instabilities_2003
- tang_laminar_2003/tang_cylindrical_2004
- aalburg_primary_2005/aalburg_deformation_2002 (has data with no cross-flow: fig. 4, figs. 7--9)
- badens_laminar_2005
- Add spray angle from osta_effect_2010 pdf p. 70
- Add salyers_spray_2010
- umemura_two-valued_2011
- Add walker_effect_2012
- moallemi_breakup_2016 (see pdf p. 10)
- omocea_breakup_2016 (unclear on length of nozzle)
- rajendran_experimental_2017/rajendran_experimental_2012
- kiaoulias_evaluation_2019 (has pressure drop so I can better estimate f!)
- patrascu_dominant-wavelength_2019 (unclear on length of nozzle)
- torregrosa_study_2020 (See section 2.5. This is DNS, but does not report any breakup quantities.)
- felis_experimental_2020

# Data validation

- bounds on as many variables as possible
- SQLite constraints
- if known liquid, check that provided fluid properties are near what is expected
- Can K_c be negative?
- Check that all basic quantities (e.g., Re, Fr, etc.) are defined
- Reject all laminar data that isn't fully developed.
- Require increasing variables like pressure, velocity, Reynolds number, Weber number, etc. in each series. Do this for certain series, resetting when a new series is encountered?
- Require a "series" for every data point. The velocity (or Reynolds number, etc.) for each successive data point here is required to be incrementing.

# Database fields

- Add D_mode (most common droplet size) to database
- chow_experimental_1980 doesn't have many data points but discusses many things done to have a high quality experiment. Add some more database fields inspired by this.
- Add motion compensation key for photos. thorne_effect_1978, taylor_water_1983 (I later read that Hoyt and Taylor had two different setups with different methods of motion compensation?)
- Add transcription error. Break into error from image resolution, error from not knowing the center of the point due to its odd shape, and error from imprecise placement of the points.
- Add Knudsen number. https://en.wikipedia.org/wiki/Knudsen_number
- Look into papers on scale effects in hydraulic structures for other ideas.
- Add whether the flow is cavitating or not.
- Add field to mark as true when data for the wrong fluid was used due to the correct data being unavailable or otherwise.
- Add duration of jet in time.
- Add Jacob number for evaporation/flashing. See Incropera p. 376, kitamura_critical_1986, cleary_flashing_2007. Alternative: how close the temperature is to boiling. nezgada_effect_1970, nezgada_effect_1970-1
- Check whether spray angles are for the near or far fields.
- Characterize velocity profile. logan_flow_1963
- Add index of refraction. corey_droplet_1969 pdf p. 37, balewski_experimental_2008
- MAYBE: Add information about cavitation inception. History of water? keller_influence_1980, ooi_scale_1985 (dissolved air content in ppm), brennen_cavitation_2013 p. 22
- MAYBE: Add other Mach number: U_0/c_l.
- MAYBE: Add distance to flat plate divided by nozzle diameter for breakup measurements. See takahashi_stability_1972 p. 75. li_study_2007 seems to be on the same topic, but I don't know if they come to the same conclusion.
- MAYBE: birkhoff_jets_1957 p. 330: > Also, dust in the air, the chemistry of the surrounding gas, and electrification may affect the final behavior decisively, by determining whether colliding drops amalgamate or rebound.
- MAYBE: Electrically ground the nozzle?
