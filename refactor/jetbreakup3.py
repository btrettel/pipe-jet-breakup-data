#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# pipe-jet-breakup-data is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

import csv
from datetime import datetime
from git import Repo
import numpy as np
import matplotlib.pyplot as plt
import pint
import os
from scipy import interpolate

# Configuration
db_file = 'pipe-jet-breakup-data.sqlite'

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
lastchangeddate     = datetime.utcfromtimestamp(headcommit.authored_date).strftime('%Y-%m-%dT%H:%M:%SZ')

# Configure Pint

ureg = pint.UnitRegistry(system='mks', auto_reduce_dimensions=True)
ureg.setup_matplotlib()

# constants
# TODO: Add uncertainties and units
# TODO: For g and P_atm, it would be wise to add uncertainties based on how likely these are to be near the standard value. I'd expect both to vary a bit.
# TODO: https://en.wikipedia.org/wiki/Theoretical_gravity#Somigliana_equation
g      = 9.80665 # m/s, <https://en.wikipedia.org/wiki/Standard_gravity>
MW_air = 28.97   # kg/kmol
MW_N2  = 2*14.007  # kg/kmol
P_atm  = 101325.  # Pa
Rbar  = 8.31446261815324e3  # J/(kmol*K), <https://en.wikipedia.org/wiki/Gas_constant>
T_zero = 273.15  # K
