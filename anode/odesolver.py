#*
# @file odesolver.py 
# This file is part of ANODE library.
#
# ANODE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ANODE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ANODE.  If not, see <http://www.gnu.org/licenses/>.
#*
from .scheme import Euler, RK2, Bosh3, RK4, RK4_alt, Dopri5


def odesolver(func, z0, options = None):
    if options == None:
        Nt = 2
        t0 = 0
    else:
        Nt = options['Nt']
        t0 = options['t0']
    print(z0.size())
    if (options['method'] == 'Euler' or options['method'] == 'euler'):
        solver = Euler(func, z0, Nt = Nt)
    elif (options['method'] == 'RK2' or options['method'] == 'rk2'):
        solver = RK2(func, z0, Nt = Nt)
    elif (options['method'] == 'fixed_bosh3'):
        solver = Bosh3(func, z0, Nt = Nt)
    elif (options['method'] == 'RK4' or options['method'] == 'rk4'):
        solver = RK4(func, z0, Nt = Nt)
    elif (options['method'] == 'RK4_alt' or options['method'] == 'rk4_alt'):
        solver = RK4_alt(func, z0, Nt = Nt)
    elif (options['method'] == 'fixed_dopri5'):
        solver = Dopri5(func, z0, Nt = Nt)
    else:
        print('error unsupported method passed')
        return
    z1 = solver.integrate(z0,t0)

    return z1
