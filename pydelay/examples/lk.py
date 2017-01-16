import numpy as np
import pylab as pl
from pydelay import dde23
import time

eqns = { 'E:C': '0.5*(1.0+ii*a)*E*n + K*E(t-tau)',
         'n'  : '(p - n - (1.0 +n) * pow(abs(E),2))/T'}

params = { 'a'  : 4.0, 
           'p'  : 1.0, 
           'T'  : 1000.0, 
           'K'  : 0.1, 
           'tau': 1000,
           'nu' : 10**-5,
           'n0' : 10.0
         }

noise = { 'E': 'sqrt(0.5*nu*(n+n0)) * (gwn() + ii*gwn())' }

dde = dde23(eqns=eqns, params=params, noise=noise)

tfinal = 20000
dde.set_sim_params(tfinal=tfinal)

# use a dictionary to set the history
thist = np.linspace(0, 1000, 10000)
Ehist = np.sin(0.01*thist)*np.exp(1.0j*0.001*thist)
nhist = np.sin(0.01*thist)-1

dic = {'t' : thist, 'E': Ehist, 'n': nhist}
dde.hist_from_arrays(dic)

dde.run()

t = dde.sol['t']
E = dde.sol['E']
n = dde.sol['n']

spl = dde.sample(-1000, 20000, 0.1)

pl.plot(t, abs(E), '.', label='calculated points')
pl.plot(spl['t'], abs(spl['E']), 'g', label='spline interpolation')
pl.plot(t[:-1], t[1:] - t[:-1], 'k', label='step size')
pl.legend()

pl.xlim((-1000, tfinal))
pl.show()
