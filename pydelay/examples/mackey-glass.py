import numpy as np
import pylab as pl
from pydelay import dde23

eqns = { 'x' : '0.25 * x(t-tau) / (1.0 + pow(x(t-tau),10.0)) -0.1*x' }

dde = dde23(eqns=eqns, params={'tau': 15})
dde.set_sim_params(tfinal=1000, dtmax=1.0, AbsTol=10**-6, RelTol=10**-3)

histfunc = {'x': lambda t: 0.5 } 
dde.hist_from_funcs(histfunc, 51)
dde.run()

sol1 = dde.sample(515, 1000, 0.1)
x1 = sol1['x']
sol2 = dde.sample(500, 1000-15, 0.1)
x2 = sol2['x']

pl.plot(x1, x2)
pl.xlabel('$x(t)$')
pl.ylabel('$x(t-15)$')

pl.show()
