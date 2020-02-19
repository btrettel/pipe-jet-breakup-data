# pipe-jet-breakup-data
Compilation of jet breakup data where the nozzle is a long pipe.

Feel free to file issues, pull requests, or contact me with any improvements or comments. Email address: ![email address](http://trettel.org/email-address.png)

# Requirements

- Python 2
- numpy
- scipy
- matplotlib
- GitPython
- Pandas

# How to use

Assuming you are using Linux:

Clone and navigate to the repository in a terminal. Then `cd scripts/`. Initially the data has not been compiled and lives in CSV spreadsheets in the `data/` directory. You can compile the data into the `outputs/data/` directory (which will be automatically created if it does not exist) by running `python prepare_data.py`.

The other scripts then read the `pipe-jet-breakup-data.pickle` file to produce plots and perform other analyses. These scripts will create files in the `outputs/` directory. You can write your own script to do your own analysis too.

# Statistics

- 27 studies
- 1350 data points in total

- liquid Weber number range: 2.0e0 to 1.9e6
- liquid Reynolds number range: 1.5e1 to 1.0e6
- turbulence intensity range: 4.4% to 12.7%
- liquid-to-gas density ratio range: 1.9e1 to 1.4e5
- nozzle length-to-diameter ratio range: 16.8 to 2207.9

- 1105 xbavg/d_0 data points
- 64 theta data points
- 48 D data points
- 90 x_i/d_0 data points
- 17 v_d_bar/vp data points
- 114 photographs
- 120 photographic regimes
- 1100 breakup length regimes
- 1194 regimes in total

# Data compilation philosophy

Trettel, B. (2019, March 22). Improving the validation of turbulent jet breakup models. https://doi.org/10.31224/osf.io/k2fnm

# TODO

- Update to Python 3
- Convert variable names to be consistent with papers
- Change regime names to be consistent with papers and not archaic
- Add more data
- Add more data validation checks
- Make many loops use zip rather than indexing
- Add data and script for pipe turbulence statistics
