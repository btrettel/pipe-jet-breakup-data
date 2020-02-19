rm -rfv ../outputs/
python prepare_data.py
python plot_data.py
python atomization.py # needs to run before plot_data_2.py and low_atm_density_regimes.py
python plot_data_2.py
python low_atm_density_regimes.py
python initially_laminar.py
python validation_problems_example.py
