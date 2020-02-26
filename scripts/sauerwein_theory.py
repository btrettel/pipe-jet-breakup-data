#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of pipe-jet-breakup-data.
# 
# pipe-jet-breakup-data is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# pipe-jet-breakup-data is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with pipe-jet-breakup-data. If not, see <https://www.gnu.org/licenses/>.

from jetbreakup import *

# Assuming homogeneous isotropic turbulence. Theory from sauerwein_theoretische_1992.

# configuration
We_l0 = 1.e2

Tubar_0_max = np.sqrt((5./6.) / We_l0)
print 'Tubar_0_max =', Tubar_0_max
assert(Tubar_0_max > 0.001)
Tubar_0_arr = np.linspace(0.0001, 0.5 * Tubar_0_max, 1e2)
S_arr = (6./5.) * We_l0 * Tubar_0_arr**2
assert(np.max(S_arr) < 1.)

# breakup length
C_TR_arr   = np.log((2./3.) * Tubar_0_arr**(-2)) # Not part of Sauerwein's theory
xbavgs_arr = C_TR_arr * We_l0**(1./2.) / (1. - S_arr)

plt.plot(100. * Tubar_0_arr, xbavgs_arr)
plt.xlabel(r'$\overline{\mathrm{Tu}}_0$, \%')
plt.ylabel(r'$\langle x_\text{b} \rangle/d_0$')
plt.grid()
plt.savefig('../outputs/figures/sauerwein_breakup_length.png')
plt.close()

# droplet size
lambda_ms_arr = np.pi / np.sqrt(0.5 * (1. - S_arr))
D_s_arr = (3. * lambda_ms_arr / 2.)**(1./3.)

plt.plot(100. * Tubar_0_arr, D_s_arr)
plt.xlabel(r'$\overline{\mathrm{Tu}}_0$, \%')
plt.ylabel(r'$D/d_0$')
plt.grid()
plt.savefig('../outputs/figures/sauerwein_droplet_size.png')
plt.close()
