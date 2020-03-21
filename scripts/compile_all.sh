#!/usr/bin/env bash
# from https://github.com/anordal/shellharden/blob/master/how_to_do_things_safely_in_bash.md
if test "$BASH" = "" || "$BASH" -uc "a=();true \"\${a[@]}\"" 2>/dev/null; then
   # Bash 4.4, Zsh
   set -euo pipefail
else
   # Bash 4.3 and older chokes on empty arrays with set -u.
   set -eo pipefail
fi
shopt -s nullglob globstar

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# pipe-jet-breakup-data is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

rm -rfv ../outputs/
python generate_database.py
python plot_data.py
python atomization.py # needs to run before plot_data_2.py and low_atm_density_regimes.py
python plot_data_2.py
python initially_laminar.py # needs to run before low_atm_density_regimes.py
python turbulent_rayleigh.py # needs to run before low_atm_density_regimes.py
python low_atm_density_regimes.py
python validation_problems_example.py
python model_comparisons.py
