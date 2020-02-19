# pipe-jet-breakup-data
Compilation of jet breakup data where the nozzle is a long pipe.

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

# TODO

- Update to Python 3
- Convert variable names to be consistent with papers
- Change regime names to be consistent with papers and not archaic
- Add more data
- Add more data validation checks
- Make many loops use zip rather than indexing
- Add data and script for pipe turbulence statistics
