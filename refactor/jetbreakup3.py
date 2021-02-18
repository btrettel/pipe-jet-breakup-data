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
lastchangeddate     = datetime.utcfromtimestamp(headcommit.authored_date).strftime('%Y-%m-%dT%H:%M:%SZ')

# Configure Pint

ureg = pint.UnitRegistry(system='mks',  auto_reduce_dimensions=True)
ureg.setup_matplotlib()

# TODO: CSV reader function which has a similar user interface to Pandas (i.e., df['name']) but in this case also handles units automatically obtained from the header.
