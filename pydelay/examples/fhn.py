import numpy as np
import pylab as pl

# import the solver
from pydelay import dde23

# define the FHN equations as in 
# Dahlem, M. A. , Hiller, G., Panchuk, A.  and Schöll, E. , Dynamics of delay-coupled excitable neural systems, Int. J. Bifur. Chaos 19, 745 (2009)
eqns = { 
        'x1': '(x1 - pow(x1,3)/3.0 - y1 + C*(x2(t-tau) - x1))/eps',
        'y1': 'x1 + a',
        'x2': '(x2 - pow(x2,3)/3.0 - y2 + C*(x1(t-tau) - x2))/eps',
        'y2': 'x2 + a'
        }

# set the parameters and the delay
params = { 
        'a': 1.3,
        'eps': 0.01,
        'C': 0.5,
        'tau': 3.0
        }


# initalise the solver
dde = dde23(eqns=eqns, params=params)

# set the simulation parameters
dde.set_sim_params(tfinal=200)

# When nothing is specified, the history for all variables 
# is initialized to 0.
#
dde.hist_from_funcs({'x1': lambda t: 1.0})

# run the simulation
dde.run()

# sample the solution with sample size dt=0.01 between 170 and 200
sol = dde.sample(170, 200, 0.01)

# plot the solution
x1 = sol['x1']
y1 = sol['y1']
x2 = sol['x2']
y2 = sol['y2']
t = sol['t']

pl.subplot(221)
pl.plot(t, x1, 'r')
pl.plot(t, y1, 'g')
pl.xlabel('$t$')
pl.ylabel('$x_1, y_1$')

pl.subplot(222)
pl.plot(x1, x2, 'r')
pl.xlabel('$x_1$')
pl.ylabel('$x_2$')

pl.subplot(223)
pl.plot(t, x2, 'r')
pl.plot(t, y2, 'g')
pl.xlabel('$t$')
pl.ylabel('$x_2, y_2$')

pl.subplot(224)
pl.plot(y1, y2, 'g')
pl.xlabel('$y_2$')
pl.ylabel('$y_2$')

pl.show()
