import numpy as np
import pylab as pl

# import the solver
from pydelay import dde23

eqns = { 
        'x': '-x + k*x(t-10) + A* f(w,t)'
        }

params = { 
        'k': 0.1,
        'w': 2.0,
        'A': 0.5
        }

# We can define a c function to be used in the equations
mycode = """
    double f(double w, double t) {
        return sin(w*t);
    }
    """

# initalise the solver
dde = dde23(eqns=eqns, params=params, supportcode=mycode)

# set the simulation parameters
dde.set_sim_params(tfinal=40)

# we can define the history as a python function
def myhist(t):
    return 0.01*t**2

dde.hist_from_funcs({'x': myhist})

# run the simulation
dde.run()

sol = dde.sample(0.01)
t = sol['t']
x = sol['x']

pl.plot(t, x)
pl.xlim((0,40))
pl.xlabel('$t$')
pl.ylabel('$x$')
pl.ylim((x.min(), x.max()))
pl.show()
