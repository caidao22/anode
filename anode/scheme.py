#*
# @file scheme.py 
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
from .time_stepper import Time_Stepper

class Euler(Time_Stepper):
    def step(self, func, t, dt, y):
        out = y + dt * func(t, y)
        return out

class RK2(Time_Stepper):
    def step(self, func, t, dt, y):
        k1 = dt * func(t, y)
        k2 = dt * func(t + dt / 2.0, y + 1.0 / 2.0 * k1)
        out = y + k2
        return out

class Bosh3(Time_Stepper):
    def step(self, func, t, dt, y):
        k1 = dt * func(t, y)
        k2 = dt * func(t + dt/2.0, y + 1.0/2.0*k1)
        k3 = dt * func(t + dt*3.0/ 4.0, y + 3.0/4.0*k2)
        out = y + 2.0/9.0*k1 + 1.0/3.0*k2 + 4.0/9.0*k3
        return out

class RK4(Time_Stepper):
    def step(self, func, t, dt, y):
        k1 = dt * func(t, y)
        k2 = dt * func(t + dt / 2.0, y + 1.0 / 2.0 * k1)
        k3 = dt * func(t + dt / 2.0, y + 1.0 / 2.0 * k2)
        k4 = dt * func(t + dt, y + k3)
        out = y + 1.0 / 6.0 * k1 + 1.0 / 3.0 * k2 + 1.0 / 3.0 * k3 + 1.0 / 6.0 * k4
        return out

class RK4_alt(Time_Stepper):
    def step(self, func, t, dt, y):
        k1 = dt * func(t,y)
        k2 = dt * func(t + dt / 3.0, y + k1 / 3.0 )
        k3 = dt * func(t + dt * 2.0 / 3.0, y + (k1 / -3.0 + k2) )
        k4 = dt * func(t + dt, y +  (k1 - k2 + k3) )
        out = y + (k1 + 3.0*k2 + 3.0*k3 + k4) / 8.0
        return out

class Dopri5(Time_Stepper):
    alpha=[1./5., 3./10., 4./5., 8./9., 1., 1.]
    beta=[
         [1./5.],
         [3./40., 9./40.],
         [44./45., -56./15., 32./9.],
         [19372./6561., -25360./2187., 64448./6561., -212./729.],
         [9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.],
         [35./384., 0, 500./1113., 125./192., -2187./6784., 11./84.],
         ]
    c_sol=[35./384., 0, 500./1113., 125./192., -2187./6784., 11./84., 0]
    c_error=[
            35./384. - 1951./21600.,
            0,
            500./1113. - 22642./50085.,
            125./192. - 451./720.,
            -2187./6784. - -12231./42400.,
            11./84. - 649./6300.,
            -1. / 60.,
            ]

    def step(self, func, t, dt, y):
        k = [dt * func(t,y)]
        stage = 6
        for i in range(stage):
            ti = t + self.alpha[i] * dt
            yi = y
            for j in range(i+1):
                yi = yi + self.beta[i][j]*k[j]
            if i < 5:
                k = k + [ dt * func(ti, yi) ]
        out = yi
        #for i in range(stage):
        #    out = out + self.c_sol[i]*k[i]
        return out
