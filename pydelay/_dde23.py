#!/usr/bin/env python
from __future__ import division

# MIT/X Consortium License
#
# pydelay, A tool for solving delay differential equations.
# Copyright (C) 2009  Valentin Flunkert
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import warnings
warnings.filterwarnings(action = 'ignore',
                        message = 'BaseException.message has been deprecated')

import numpy as np
from scipy import weave
from scipy.interpolate import splrep, splev, spalde
import math
import sys
import hashlib
import time
import re
from helper import _parenthesis_balancedQ, _symbols_allowedQ,\
        gen_disconts, gwn_code, assure_array, isrealnum, isnum

class dde23:
    """This class translates a DDE to C and solves it using the Bogacki-Shampine method.

    *Attributes of class instances:*

    **For user relevant attributes:**

    `self.sol`
        Dictionary storing the solution (when the simulation has finished).
        The keys are the variable names and ``'t'`` corresponding to the 
        sampled times.

    `self.discont` 
        List of discontinuity times. This is generated from the occurring
        delays by propagating the discontinuity at ``t=0``. The solver will step on these 
        discontinuities. If you want the solver to step onto certain times they can
        be inserted here. 

    `self.rseed`  
        Can be set to initialise the random number generator with a specific
        seed. If not set it is initialised with the time.

    `self.hist`   
        Dictionary with the history. Don't manipulate the history arrays directly!
        Use the provided functions to set the history.

    `self.Vhist`  
        Dictionary with the time derivatives of the history.


    **For user less relevant attributes:**

    `self.delays` 
        List of the delays occurring in the equations.

    `self.chunk`  
        When arrays become to small they are grown by this number.

    `self.spline_tck` 
        Dictionary which stores the tck spline representation of the solutions.
        (see ``scipy.interpolate``)

    `self.eqns`   
        Stores the eqn dictionary.

    `self.params` 
        Stores the parameter dictionary.

    `self.simul`  
        Dictionary of the simulation parameters.

    `self.noise`  
        Stores the noise dictionary.

    `self.debug`  
        Stores the debug flag.

    `self.delayhashs` 
        List of hashs for each delay (this is used in the generated C-code).

    `self.vars`    
        List of variables extracted from the eqn dictionary keys.

    `self.types`  
        Dictionary of C-type names of each variable.

    `self.nptypes` 
        Dictionary of numpy-type names of each variable.
    """

    def __init__(self, eqns, params=None, noise=None, 
            supportcode = '', debug=False):
        """Initialise the solver.

        `eqns`  
            Dictionary defining for each variable the derivative. 
            Delays are written as as ``(t-...)``
            example::

                eqns = {
                    'y1': '- y1 * y2(t-tau1) + y2(t-tau2)',
                    'y2': 'a * y1 * y2(t-tau1) - y2',
                    'y3': 'y2 - y2(t-tau2)'
                    }

            You can also directly use numbers or combination of parameters as delays::

                eqns = {
                    'x1': '-a*x1 + x1(t - 1.0)',
                    'x2': 'x2-b*x1(t-2.0*a+b)
                    }

            At the moment only constant delays are supported.

            The string defining the equation has to be a valid C expression, i.e.,
            use ``pow(a,b)`` instead of ``a**b`` etc. (this might change in the future)::

                eqns = {'y': '-2.0 * sin(t) * pow(y(t-tau), 2)'}

            Complex variable can be defined using ``:C`` or ``:c`` in the variable name.
            The imaginary unit can be used through ``ii`` in the equations::

                eqns = {'z:C': '(-la + ii * w0) * z' }

        `params` 
            Dictionary defining the parameters (including delays) used in eqns.
            example::

                 params = {
                    'a'   : 1.0, 
                    'tau1': 1.0, 
                    'tau2': 10.0
                    }

        `noise` 
            Dictionary for noise terms. The function ``gwn()`` can be accessed in 
            the noise string and provides a Gaussian white noise term of unit variance.
            example::

                 noise = {'x': '0.01*gwn()'}

        `debug` 
            If set to ``True`` the solver gives verbose output to stdout while running.

        """

        self.eqns = eqns     # the equations
        if params == None:
            params = {}
        self.params = params # the parameters for the equations
        if noise == None:
            noise = {}
        self.noise = noise   
        self.supportcode = supportcode
        self.debug = debug
        self.rseed = int(str(int(time.time()*10))[3:]) # generate a random seed from the time
        self.delayhashs = [] # the md5 hash values of the delay terms
        self.code = ''       # the generated c code
        self.vars = []
        self.types = {}      # the types of the variables: e.g. {'x': 'double', 'z': 'Complex'}
        self.nptypes = {}    # the numpy types e.g. {'x': 'NPY_DOUBLE', 'z': 'NPY_CDOUBLE'}
        self.delays = []
        self.discont = []
        self.simul = {}      # the simulation parameters "tfinal", "AbsTol",...
        self.sol = {}        # the solutions
        self.hist = {}       # the history
        self.Vhist = {}      # the time derivative of the history
        self.spline_tck = {} # stores the tck spline representation of the solutions
        self.chunk = 10000   # grow arrays by this number if too small
        self._maxdelay = 1
        self._update()
        self.__hist_set = False #history has been set for at least one variable
        # initialise the histories to zero
        maxdelay = self._maxdelay
        nn = 101
        self.hist['t'] = np.linspace(-maxdelay, 0, nn, endpoint=True)
        self._initzeros(nn)

    def _initzeros(self, nn):
        """initialize the history of all variables to zero arrays of len nn.
        The time history is not set."""
        for var in self.vars:
            self.hist[var]  = np.zeros(nn) 
            self.Vhist[var] = np.zeros(nn)
            if self.types[var] == 'Complex':
                self.hist[var] = self.hist[var] + 0.0j
                self.Vhist[var] = self.Vhist[var] + 0.0j 
        self._nstart = nn-1

    def _update(self):
        self.vars = []
        self.types = {}
        self.nptypes = {}
        self.delays = []
        self.delayhashs = []
        for vardef in self.eqns:
            if ':' in vardef:
                try:
                    vname, vartype = vardef.split(':')
                except:
                    raise AssertionError,\
                            '%s is not a valid variable definition. Use "varname" or "varname:C"'%vardef
                if vartype in ['C', 'c']:
                    self.nptypes[vname] = 'NPY_CDOUBLE'
                    self.types[vname] = 'Complex'
                else:
                    raise AssertionError,\
                            'Unknown type in variable definition "%s"'%vardef
                self.vars.append(vname)
            else:
                self.types[vardef] = 'double'
                self.nptypes[vardef] = 'NPY_DOUBLE'
                self.vars.append(vardef)

        for vardef in self.eqns:
            var = vardef.split(':')[0]
            eq = self.eqns[vardef]
            assert _symbols_allowedQ(eq), \
                    r"""Forbidden symbol in eqn for %s:

                    %s

                    Forbidden symbols are: ** ^ { } [ ] & # ! @ %% ? ; > < | \
                    """%(var, eq)
            assert _parenthesis_balancedQ(eq), \
                    """Eqn for %s has unbalanced parenthesis:

                    %s
                    """%(var, eq)

        self._check_simul_pars()
        self._check_params()
        self._check_eqns()
        self._generate_code()

    def _check_eqns(self):
        for vardef in self.eqns:
            var = vardef.split(':')[0]
            eqn = self.eqns[vardef]
            funcsymbols = ['sin(', 'cos(', 'tan(', 'asin(', 'acos(', 'atan(', 'atan2(', 'cosh(',
                    'sinh(', 'exp(', 'abs(', 'fabs(', 'tan(', 'conj(', 'real(', 'imag(', 'tanh(', 
                    'log(', 'pow(', 'sqrt(', 'ceil(', 'floor(', 'log10(']    

            symbols = ['t','ii','Heavi('] + self.params.keys() + self.vars + funcsymbols + ['%s('%var for var in self.vars]
            symbpattern = r'([a-zA-Z]+\w*[\(]{0,1})'
            res = re.findall(symbpattern, eqn)
            for match in res:
                if self.supportcode == '':
                    assert match in symbols, "Undefined string '%s' in eqn for '%s':\n\n%s\n"%(match, var, eqn)
                else:
                    if match not in symbols and match not in self.supportcode:
                        sys.stderr.write("Warning: Suspicious string '%s' in eqn for '%s':\n\n%s\nIs this a function defined in your supportcode?\n"%(match, var, eqn))

    def _check_simul_pars(self):
        for key, value in zip(self.simul, self.simul.values()):
            assert key in ['AbsTol', 'RelTol', 'dtmin', 'dtmax', 'dt0', 'MaxIter', 'tfinal'],\
                    "'%s' is not a valid simulation parameter"%key
            assert isrealnum(value), \
                 "'%s' is set to '%s', has to be a positive number"%(key, value)
            assert value > 0, "'%s' is set to '%s', has to be a positive number"%(key, value)

    def _check_params(self):
        for key, value in zip(self.params, self.params.values()):
            assert isnum(value), \
                 "The parameter '%s' is set to '%s'. All parameters have to be numbers."%(key, value)

    def set_sim_params(self, tfinal=100, AbsTol=1e-6, RelTol=1e-3, 
            dtmin=1e-6, dtmax=None, dt0=None, MaxIter=1e9):
        """
        `tfinal` 
            End time of the simulation (the simulation always starts at ``t=0``).

        `AbsTol`, `RelTol` 
            The relative and absolute error tolerance. If the estimated error `e` 
            for a variable `y` obeys `e <= AbsTol + RelTol*|y|` then the step is accepted.  
            Otherwise the step will be repeated with a smaller step size.

        `dtmin`, `dtmax` 
            Minimum and maximum step size used.

        `dt0` 
            initial step size

        `MaxIter` 
            maximum number of steps. The simulation stops if this is reached.
        """
        if dtmax != None and dtmax < 0.9*self._mindelay:
            self.simul['dtmax']  = dtmax
        else:
            self.simul['dtmax']  = min([tfinal/50.0, self._mindelay*0.9])
        if dt0 != None:
            self.simul['dt0']  = dt0
        else:
            self.simul['dt0']  = dtmin
        self.simul['dtmin']  = dtmin
        self.simul['AbsTol'] = AbsTol
        self.simul['RelTol'] = RelTol
        order = 3
        for eq in self.eqns.values():
            if "'" in eq:
                try:
                    order = int(math.ceil(tfinal/self._mindelay))
                except:
                    pass
        discont = gen_disconts(0, tfinal, self.delays, order=order)[1:] #remove 0 from discont
        if tfinal not in discont:
            discont.append(tfinal)
        self.discont = discont

        self.simul['tfinal']   = tfinal
        self.simul['MaxIter'] = MaxIter

        self._update()

    def hist_from_dict(self, *args):
        """This function has been replaced by hist_from_arrays (with an s) 
        which takes a dictionary of arrays."""
        assert False, """This function has been replaced by hist_from_arrays (with an s) \
                which takes a dictionary of arrays."""


    def hist_from_funcs(self, dic, nn=101):
        """
        Initialise the histories with the functions stored in the dictionary `dic`.
        The keys are the variable names.  The function will be called as ``f(t)`` 
        for ``t`` in ``[-maxdelay, 0]`` on `nn` samples in the interval.

        This function provides the simplest way to set the history.
        It is often convenient to use python ``lambda`` functions for ``f``.
        This way you can define the history function in place.

        If any variable names are missing in the dictionaries, the history of these
        variables is set to zero and a warning is printed. If the dictionary contains 
        keys not matching any variables these entries are ignored and a warning is 
        printed.

        Example: Initialise the history of the variables ``x`` and ``y`` with 
        ``cos`` and ``sin`` functions using a finer sampling resolution::

            from math import sin, cos

            histdic = {
                'x': lambda t: cos(0.2*t),
                'y': lambda t: sin(0.2*t)
            }

            dde.hist_from_funcs(histdic, 500)
        """
        #TODO assert dic and not function
        maxdelay = self._maxdelay
        t = np.linspace(-maxdelay, 0, nn, endpoint=True) 
        
        for key in dic:
            if key not in self.vars:
                sys.stderr.write('Warning: The key %s in the history dictionary is not a variable name.\
                        It will be ignored.\n'%key)

        arraydic = {'t': t}
        for key in dic:
            vecf = np.vectorize(dic[key])
            arraydic[key] = vecf(t)

        self.hist_from_arrays(arraydic)

    def hist_from_func(self, *args):
        """This function has been replaced by hist_from_funcs (with an s) 
        which takes a dictionary of functions."""
        assert False, """This function has been replaced by hist_from_funcs (with an s)\
                which takes a dictionary of functions."""


    def hist_from_array(self, *args):
        """This function has been replaced by hist_from_arrays (with an s) which takes a dictionary
        and is more secure."""
        assert False, """This function has been replaced by hist_from_arrays which takes a dictionary of arrays.
        Please use this function."""

    def hist_from_arrays(self, dic, useend=True):
        """Initialise the history using a dictionary of arrays with variable names as keys.  
        Additionally a time array can be given corresponding to the key ``t``. 
        All arrays in `dic` have to have the same lengths.

        If an array for ``t`` is given the history is interpreted as points 
        ``(t,var)``. Otherwise the arrays will be evenly spaced out over the interval 
        ``[-maxdelay, 0]``.

        If useend is True the time array is shifted such that the end time is
        zero. This is useful if you want to use the result of a previous simulation 
        as the history.

        If any variable names are missing in the dictionaries, the history of these
        variables is set to zero and a warning is printed. 
        If the dictionary contains keys not matching any variables (or ``'t'``) these
        entries are ignored and a warning is printed.

        Example:::
            
            t = numpy.linspace(0, 1, 500)
            x = numpy.cos(0.2*t)
            y = numpy.sin(0.2*t)

            histdic = {
                't': t,
                'x': x,
                'y': y
            }
            dde.hist_from_arrays(histdic)

        """
        if self.__hist_set:
            sys.stderr.write("Warning: The history has already been set. It will be overwritten.\n")
        self.__hist_set = True

        maxdelay = self._maxdelay

        for key in dic:
            dic[key] = assure_array(dic[key])
        laengen = map(len, dic.values())
        laengen.sort()
        assert laengen[0] == laengen[-1],\
                'Error: not all arrays in dic have the same length'

        if 't' in dic:
            t = assure_array(dic['t']).copy()
            if useend:
                t = t-t[-1]
        else:
            nn = laengen[0]
            t = np.linspace(-maxdelay, 0, nn, endpoint=True)
        assert t[0] <= -maxdelay,\
            ('Error: The history has to reach back at least to the maximum delay %s.'%maxdelay +
             'The time array you gave spans [%s, %s]'%(t[0], t[-1]))
        assert t[-1] >= 0, ('Error: The history has to reach up to t=0.' +
            'The time array you gave spans [%s, %s]'%(t[0], t[-1]))

        # only copy the data that are needed
        delta = 0.001
        indices = (-maxdelay <= t) * (t<=0)
        t = t[indices]
        nn = len(t)

        tsample = np.zeros(nn+1)
        tsample[1:] = t
        tsample[0] = -maxdelay - delta * (t[1]-t[0])
        self.hist['t'] = tsample

        for key in dic:
            if key != 't' and key not in self.vars:
                sys.stderr.write('Warning: The key %s in the history dictionary is not a variable name.\
                            It will be ignored.\n'%key)

        for var in self.vars:
            if var in dic:
                if self.types[var] == 'Complex':
                    y = assure_array(dic[var])[indices]+0.0j
                    tck_real = splrep(t, y.real)
                    tck_imag = splrep(t, y.imag)
                    self.hist[var] = splev(tsample, tck_real) + 1.0j*splev(tsample, tck_imag)
                    self.Vhist[var] = splev(tsample, tck_real, der=1) + 1.0j*splev(tsample, tck_imag, der=1)
                else:
                    y = assure_array(dic[var])[indices]
                    tck = splrep(t, y)
                    self.hist[var] = splev(tsample, tck)
                    self.Vhist[var] = splev(tsample, tck, der=1)
            else:
                sys.stderr.write('History for %s not given. Will be set to zero.\n'%var)
                if self.types[var] == 'Complex':
                    self.hist[var] = np.zeros(nn+1) +0.0j
                    self.Vhist[var] = np.zeros(nn+1) +0.0j
                else:
                    self.hist[var] = np.zeros(nn+1)
                    self.Vhist[var] = np.zeros(nn+1)
        self._nstart = nn

    def _generate_code(self):
        """Generate the c code.
        The code is stored in 
           
        self.code and self.includes

        """

        ## These two functions are called by re.sub:
        # replace the delay terms by the appropriate function call
        def repldelay(matchobj):
            var = matchobj.group(1)
            delay = matchobj.group(2)
            s = delay
            for par in self.params:
                s=s.replace('PAR%s'%par, str(self.params[par]))
            try:
                delval = -eval(s)
            except:
                raise AssertionError, """Error: The delay term "%s" is an undefined expression. 
                (State and time dependent delays are not yet supported.)"""%delay
            assert delval >= 0, 'Error: The delay term  "%s" is in the future.'%delay

            if delval > 0:
                dhash = hashlib.md5(delay).hexdigest()
                self.delayhashs.append(dhash)
                self.delays.append(delval)
                return 'interp_%(var)s(t %(delay)s, &n_%(dhash)s) '%locals()
            else:
                return ' %(var)s '%locals()

        #replace the derivative delay terms ...
        def Vrepldelay(matchobj):
            var = matchobj.group(1)
            delay = matchobj.group(2)
            dhash = hashlib.md5(delay).hexdigest()
            self.delayhashs.append(dhash)
            s = delay
            for par in self.params:
                s=s.replace('PAR%s'%par, str(self.params[par]))
            try:
                delval = -eval(s)
            except:
                raise AssertionError, """Error: The delay term "%s" is an undefined expression. 
                (State and time dependent delays are not yet supported.)"""%delay
            assert delval > 0, """Error: The delay term  "%s" is in the future. 
            Or the eqn involves a derivative term without delay."""%delay
            self.delays.append(delval)
            dhash = hashlib.md5(delay).hexdigest()
            self.delayhashs.append(dhash)
            return 'interp_%(var)s(t %(delay)s, &n_%(dhash)s)'%locals()

        # generate the C equations

        # TODO: to speed things up 
        # - first collect all delay terms, then check the delays and build replace patterns 
        #   for each delay and go through

        eqnkeys = [key.split(':')[0] for key in self.eqns]
        eqnvals = self.eqns.values()

        CEQNSk1line = ';'.join(eqnvals)
        CEQNSk2line = ';'.join(eqnvals)
        CEQNSk3line = ';'.join(eqnvals)
        CEQNSk4line = ';'.join(eqnvals)

        noisekeys = [key.split(':')[0] for key in self.noise]
        noisevals = self.noise.values()

        CNOISEline = ';'.join(noisevals)

        for prm in self.params:
            compiled_pattern = re.compile(r'\b%s\b'%prm)
            CEQNSk1line = compiled_pattern.sub(r' PAR%s '%prm, CEQNSk1line)
            CEQNSk2line = compiled_pattern.sub(r' PAR%s '%prm, CEQNSk2line)
            CEQNSk3line = compiled_pattern.sub(r' PAR%s '%prm, CEQNSk3line)
            CEQNSk4line = compiled_pattern.sub(r' PAR%s '%prm, CEQNSk4line)
            CNOISEline  = compiled_pattern.sub(r' PAR%s '%prm, CNOISEline)

        for v in self.vars:
            compiled_pattern = re.compile(r'\b(%s)\s*\(\s*t\s*(.*?)\s*\)'%v)
            CEQNSk1line = compiled_pattern.sub(repldelay, CEQNSk1line)
            CEQNSk2line = compiled_pattern.sub(repldelay, CEQNSk2line)
            CEQNSk3line = compiled_pattern.sub(repldelay, CEQNSk3line)
            CEQNSk4line = compiled_pattern.sub(repldelay, CEQNSk4line)
            CNOISEline  = compiled_pattern.sub(repldelay, CNOISEline)

#            compiled_pattern = re.compile(r"\b(%s)'\s*\(\s*t\s*(.*?)\s*\)"%v)
#            CEQNSk1line = compiled_pattern.sub(Vrepldelay, CEQNSk1line)
#            CEQNSk2line = compiled_pattern.sub(Vrepldelay, CEQNSk2line)
#            CEQNSk3line = compiled_pattern.sub(Vrepldelay, CEQNSk3line)
#            CEQNSk4line = compiled_pattern.sub(Vrepldelay, CEQNSk4line)
#            CNOISEline  = compiled_pattern.sub(Vrepldelay, CNOISEline)
            
            compiled_pattern = re.compile(r'\b%s\b'%v)
            CEQNSk1line = compiled_pattern.sub(r' pt_%s_ar[SIM_n] '%v, CEQNSk1line)
            CEQNSk2line = compiled_pattern.sub(r' (pt_%s_ar[SIM_n] + 0.5 * dt * k1%s) '%(v,v), CEQNSk2line)
            CEQNSk3line = compiled_pattern.sub(r' (pt_%s_ar[SIM_n] + 0.75 * dt * k2%s) '%(v,v), CEQNSk3line)
            CEQNSk4line = compiled_pattern.sub(r' TEMP%s '%v, CEQNSk4line)
            CNOISEline  = compiled_pattern.sub(r' pt_%s_ar[SIM_n] '%v, CNOISEline)

        CEQNSk1 = dict(zip(eqnkeys, CEQNSk1line.split(';')))
        CEQNSk2 = dict(zip(eqnkeys, CEQNSk2line.split(';')))
        CEQNSk3 = dict(zip(eqnkeys, CEQNSk3line.split(';')))
        CEQNSk4 = dict(zip(eqnkeys, CEQNSk4line.split(';')))
        CNOISE  = dict(zip(noisekeys, CNOISEline.split(';')))

        self.delayhashs = list(set(self.delayhashs))
        self.delays = list(set(self.delays))
        try:
            self._maxdelay = max(self.delays)
            self._mindelay = min(self.delays)
        except:
            self._maxdelay = 5.0
            self._mindelay = 5.0


        code =\
        """
#ifdef MANUAL
int main()
{
double dtmin = 0.0001;
double dtmax = 0.1;
double tfinal = 100.0;
double RelTol = 1.0E-3;
double AbsTol = 1.0E-6;
int chunk = 10000;
int nstart = 101;
double dt0 = 0.01;
double maxdelay = ...; // Set the maximum delay here !!!
long unsigned int MaxIter = 10000000;
int NumOfDiscont = 4;
double discont[4] = {maxdelay, 2*maxdelay, 3*maxdelay, tfinal};
"""

        for par in self.params:
            if 'complex' in str(type(self.params[par])):
                code += "Complex PAR%s = %s;\n"%(par, self.params[par])
            else:
                code += "double PAR%s = %s;\n"%(par, self.params[par])


        code +=\
        """
#endif

t = 0.0; 
long unsigned int i;
long unsigned int NumberOfMinSteps = 0;
int nextdsc = 0;
int hitdsc = false;
double dist;
int TakingMinStep = 0;
SIM_n = nstart; 
SIM_size = SIM_n + chunk;
double RelErr;
double thresh = AbsTol/RelTol;
double dt = dt0;
srand((unsigned)RSEED);
"""

        for dhash in self.delayhashs:
            code += 'n_%s = SIM_n;\n'%dhash

        for var in self.vars:
            code +=\
            """

pt_%(var)s_ar  = (%(type)s*) malloc((SIM_n+chunk) * sizeof(%(type)s));
pt_%(var)s_Var = (%(type)s*) malloc((SIM_n+chunk) * sizeof(%(type)s));

            """%{'var': var, 'type': self.types[var]}

        code +=\
        """
pt_t_ar = (double *) malloc((SIM_n+chunk) * sizeof(double));

#ifndef MANUAL
        """



        for var in self.vars:
            code +=\
            """
for(i = 0; i < nstart+1; i++) {
    pt_%(var)s_ar[i]  = hist%(var)s_ar[i];
    pt_%(var)s_Var[i] = Vhist%(var)s_ar[i];
}
            """%{'var': var, 'type': self.types[var]}

        code +="""
for(i = 0; i < nstart+1; i++) 
    pt_t_ar[i] = Thist_ar[i];
#endif

#ifdef MANUAL
for(i = 0; i < nstart+1; i++) 
    pt_t_ar[i]  = -maxdelay*(nstart-i)/nstart; 

/* set the history here when running generated code directly */
"""

        for var in self.vars:
            code +=\
            """ 
for(i = 0; i < nstart+1; i++) {
    pt_%(var)s_ar[i]  = 0.2; // history value for the variable
    pt_%(var)s_Var[i] = 0.0; // history of the derivatives (0.0 if constant history)
}
"""%{'var': var, 'type': self.types[var]}

        code +=\
        """
#endif

        """

        # initalise some variables
        for var in self.vars:
            code += '\n%s %s;\n'%(self.types[var],\
                    ', '.join(['k%s%s'%(i, var) for i in [1,2,3,4]]))
            code += '%s TEMP%s;\n'%(self.types[var], var)
            code += '%s ERROR%s;\n'%(self.types[var], var)

        code += '//k1 need to be calculated only for the first step\n' +\
                '//due to the FSAL property k1(n+1)=k4(n)\n'
        for var in self.vars: 
            code += '\tk1%s'%var + ' = ' + CEQNSk1[var] + ';\n'

        code += """
while((t <= tfinal || hitdsc) && SIM_n-nstart <= MaxIter) {

"""

        # evaluate the equations in each step
        # and set the time correctly (for non-autonomous equations)
        code += '\tt += dt * 0.5;\n'
        for var in self.vars: 
            code += '\tk2%s'%var + ' = ' + CEQNSk2[var] + ';\n'
        code += '\tt += dt * 0.25;\n'
        for var in self.vars: 
            code += '\tk3%s'%var + ' = ' + CEQNSk3[var] + ';\n'
        for var in self.vars:
            code += '\tTEMP%(var)s = pt_%(var)s_ar[SIM_n] + dt * 1.0/9.0 * (2.0*k1%(var)s + 3.0*k2%(var)s + 4.0*k3%(var)s);\n'%{'var': var}
        code += '\tt += dt * 0.25;\n'
        for var in self.vars:
            code += '\tk4%s'%var + ' = '  + CEQNSk4[var] + ';\n'
            code += '\tERROR%(var)s = dt/72.0 * (-5.0*k1%(var)s + 6.0*k2%(var)s + 8.0*k3%(var)s - 9.0*k4%(var)s);\n'%{'var': var}

            if self.debug:
                code += 'std::cout << "t = " << t << "\tdt = " << dt << "\tk1 = " << k1%(var)s << "\tk2 = " << k2%(var)s << "\tk3 = " << k3%(var)s << "\t k4 = " << k4%(var)s << std::endl;\n'%{'var': var}
       
        code += '\tRelErr = 0.0;\n'
        for var in self.vars:
            code += '\tERROR%(var)s = ERROR%(var)s/MAX( MAX(ABS%(var)s(pt_%(var)s_ar[SIM_n]), ABS%(var)s(TEMP%(var)s)), thresh);\n\n'%{'var': var}

        for var in self.vars:
            code += '\tRelErr = MAX(RelErr, ABS%(var)s(ERROR%(var)s));\n'%{'var': var}
            if self.debug:
                code += 'std::cout << "\t RelErr = " << RelErr << "\t RelTol = " << RelTol << std::endl;\n'

        code += '\tif(RelErr <= RelTol || TakingMinStep ) {\n'

        for var in self.vars:
            if var in self.noise:
                code += '\t\tpt_%(var)s_ar[SIM_n+1] = TEMP%(var)s + sqrt(dt) * ( %(noise)s );\n'%{'var': var, 'noise': CNOISE[var]}+\
                        '\t\tpt_%(var)s_Var[SIM_n+1] = k4%(var)s;\n'%{'var': var}
            else:
                code += '\t\tpt_%(var)s_ar[SIM_n+1] = TEMP%(var)s;\n\t\tpt_%(var)s_Var[SIM_n+1] = k4%(var)s;\n'%{'var': var}
            code += '\t\tk1%(var)s = k4%(var)s; //FSAL\n'%{'var': var}

        code +=\
        """
        pt_t_ar[SIM_n+1] = t;
        SIM_n++;
        if(SIM_n - nstart > MaxIter) {
            std::cerr << "Warning: MaxIter reached! EndTime of simulation: " << t << std::endl << std::flush;
        }
        //std::cout << "pass: " << pow(RelTol/RelErr, 0.3333333333333) << std::endl;
        dt = dt * MAX(0.5, 0.8*pow(RelTol/RelErr, 0.3333333333333));

        // hit discontinuities
        hitdsc = false;
        if(nextdsc<NumOfDiscont) {
            dist = discont[nextdsc] - t;
            //std::cout << t << "\t" << discont[nextdsc] << "\t" << dist << std::endl;
            if(dist <= MIN(1.1*dt, dtmax)) {
                dt = dist;
                nextdsc++;
                hitdsc=true;
            }
            else if(dist <= 2*dt) {
                dt = 0.5 * dist;
            }
        }

    } else {
        //std::cout << "not passed" << std::endl;
        t-= dt;
        dt = 0.5*dt;
    }
    
    if(dt < dtmin && !hitdsc)  {
        //TODO: fix this    
        //if(NumberOfMinSteps == 0)
        //    std::cerr << "Warning: step size very small" << std::endl;
        NumberOfMinSteps++;
        TakingMinStep = true;
        dt = dtmin;
    } else {
        TakingMinStep = false;
    }
    
    if(dt > dtmax) 
        dt = dtmax;

            """

        # grow arrays if end is reached
        code +=\
            """
    //grow arrays if they are too small
    if(SIM_n+1 == SIM_size)
    {
        SIM_size += chunk;

        double *p;
        p = (double *) realloc(pt_t_ar, (SIM_size) * sizeof(double));
        if (!p) {
            std::cout << "realloc fail, there is probably not enough memory" << std::endl;
            exit(1);
        } 
        pt_t_ar = p;
            """

        for var in self.vars:
            code +=\
            """
        %(type)s *p%(var)s;
        p%(var)s = (%(type)s *) realloc(pt_%(var)s_ar, (SIM_size) * sizeof(%(type)s));
        if (!p%(var)s) {
            std::cout << "realloc fail, there is probably not enough memory" << std::endl;
            exit(1);
        } 
        pt_%(var)s_ar = p%(var)s;

        p%(var)s = (%(type)s *) realloc(pt_%(var)s_Var, (SIM_size) * sizeof(%(type)s));
        if (!p%(var)s) {
            std::cout << "realloc fail, there is probably not enough memory" << std::endl;
            exit(1);
        } 
        pt_%(var)s_Var = p%(var)s;
            """%{'var': var, 'type': self.types[var]}

        code += """
    }
}
double *p;
p = (double *) realloc(pt_t_ar, (SIM_n) * sizeof(double));
if (!p) 
    std::cout << "realloc fail when shrinking arrays" << std::endl;
else
    pt_t_ar = p;

"""

        for var in self.vars:
            code +=\
            """
%(type)s *p%(var)s;
p%(var)s = (%(type)s *) realloc(pt_%(var)s_ar, (SIM_n) * sizeof(%(type)s));
if (!p%(var)s) 
    std::cout << "realloc fail when shrinking arrays" << std::endl;
else
    pt_%(var)s_ar = p%(var)s;

p%(var)s = (%(type)s *) realloc(pt_%(var)s_Var, (SIM_n) * sizeof(%(type)s));
if (!p%(var)s) 
    std::cout << "realloc fail when shrinking arrays" << std::endl;
else
    pt_%(var)s_Var = p%(var)s;

            """%{'var': var, 'type': self.types[var]}

        code +=\
        """
#ifndef MANUAL
PyObject *erg;
erg = PyDict_New();

PyObject *myarray;
double *array_buf;
npy_intp dim = SIM_n;

myarray = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);

array_buf = (double *) PyArray_DATA(myarray);

for(i = 0; i < SIM_n; i++) 
    *array_buf++ = pt_t_ar[i];
free(pt_t_ar);

PyDict_SetItem(erg, PyString_InternFromString("t"), myarray);
Py_DECREF(myarray);

        """

        for var in self.vars:
            code +=\
            """
PyObject *myarray_%(var)s;
//PyObject *myarray_V%(var)s;
%(type)s *array_buf%(var)s; 
//%(type)s *array_bufV%(var)s;

myarray_%(var)s = PyArray_SimpleNew(1, &dim, %(nptype)s);
//myarray_V%(var)s = PyArray_SimpleNew(1, &dim, %(nptype)s);

array_buf%(var)s = (%(type)s *) PyArray_DATA(myarray_%(var)s);
//array_bufV%(var)s = (%(type)s *) PyArray_DATA(myarray_V%(var)s);

for(i = 0; i < SIM_n; i++) {
    *array_buf%(var)s++ = pt_%(var)s_ar[i];
    //*array_bufV%(var)s++ = pt_%(var)s_Var[i];
}
free(pt_%(var)s_ar);
free(pt_%(var)s_Var);

PyDict_SetItem(erg, PyString_InternFromString("%(var)s"), myarray_%(var)s);
//PyDict_SetItem(erg, PyString_InternFromString("V%(var)s"), myarray_V%(var)s);
Py_DECREF(myarray_%(var)s);
//Py_DECREF(myarray_V%(var)s);

            """ % {'var': var, 'type': self.types[var], 'nptype': self.nptypes[var]}
        code += """
if(NumberOfMinSteps>1)
    std::cerr << "Number of Minimum steps taken: " << NumberOfMinSteps << std::endl;

return_val =  erg;
Py_DECREF(erg);
#endif
"""
        code +=\
            """ 
#ifdef MANUAL
for(i = 0; i < SIM_n; i++) 
    std:: cout << pt_t_ar[i] << "\t" << %s << std::endl;
    return 0;
}
#endif
"""%(' << "\t" << '.join(['pt_%s_ar[i]'%var for var in self.vars]))

        includes=\
        """
#include<complex>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<iostream>

/****** uncomment the following line to run generated source code directly ***************/
//#define MANUAL

#define false 0;
#define true 1;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#define SQUARE(X) ((X) * (X))

double Heavi(double tin) {
    if(tin>=0)
        return 1.0;
    else
        return 0.0;
}

typedef std::complex<double> Complex;
const Complex ii(0.0, 1.0);
double t;
double *pt_t_ar;
long unsigned int SIM_n;
long unsigned int SIM_size;
"""

        for var in self.vars: 
            if self.types[var] == 'double':
                absfun = 'fabs'
            else:
                absfun = 'abs'
            includes += '#define ABS%s(X) (%s(X))\n'%(var, absfun)

        includes += '\n/* positions of lags in the arrays */\n'
        for dhash in self.delayhashs:
            includes += 'long unsigned int n_%s;\n'%dhash
        includes += '/* arrays holding the variabls */\n'
        for var in self.vars:
            includes += '%s *pt_%s_ar;\n'%(self.types[var], var)
            includes += '%s *pt_%s_Var;\n'%(self.types[var], var)

        for var in self.vars:
            includes +=\
            """
inline %(vartype)s hermite_%(var)s(const double &, const double &, const %(vartype)s &, const %(vartype)s &, 
                                                   const double &, const %(vartype)s &, const %(vartype)s &);
inline %(vartype)s dt_hermite_%(var)s(const double &, const double &, const %(vartype)s &, const %(vartype)s &, 
                                                   const double &, const %(vartype)s &, const %(vartype)s &);\n

            """%{'var': var, 'vartype': self.types[var]}

        tempcode=\
        """
%(vartype)s interp_%(var)s(double t, long unsigned int *n0) 
{
    while(pt_t_ar[*n0+1] < t)
        (*n0)++;
    while(pt_t_ar[*n0] > t) 
        (*n0)--;

    return hermite_%(var)s(t, pt_t_ar[*n0],   pt_%(var)s_ar[*n0],   pt_%(var)s_Var[*n0], 
                              pt_t_ar[*n0+1], pt_%(var)s_ar[*n0+1], pt_%(var)s_Var[*n0+1]);
}

//interpolation of d(%(var)s)/dt
%(vartype)s dt_interp_%(var)s(double t, long unsigned int *n0) 
{
    while(pt_t_ar[*n0+1] < t)
        (*n0)++;
    while(pt_t_ar[*n0] > t) 
        (*n0)--;

    return dt_hermite_%(var)s(t, pt_t_ar[*n0],   pt_%(var)s_ar[*n0],   pt_%(var)s_Var[*n0], 
                              pt_t_ar[*n0+1], pt_%(var)s_ar[*n0+1], pt_%(var)s_Var[*n0+1]);
}

//hermite interpolation
inline %(vartype)s hermite_%(var)s(const double &t, const double &tn, const %(vartype)s &Xn, const %(vartype)s &Vn, 
                                     const double &tnp1, const %(vartype)s &Xnp1, const %(vartype)s &Vnp1) 
{
    double h = tnp1 - tn;
    double s = (t - tn) / h;

    return (1.0 + 2 * s) * SQUARE(s-1.0) * Xn + (3.0 - 2 * s) * SQUARE(s) * Xnp1 + h * s * SQUARE(s-1.0) *Vn + h * (s - 1) * SQUARE(s) *Vnp1;
}

inline %(vartype)s dt_hermite_%(var)s(const double &t, const double &tn, const %(vartype)s &Xn, const %(vartype)s &Vn, 
                                     const double &tnp1, const %(vartype)s &Xnp1, const %(vartype)s &Vnp1) 
{
    double h = tnp1 - tn;
    double s = (t - tn) / h;

    return (1.0-4.0*s+3.0*SQUARE(s))*Vn + s*(3.0*s-2.0)*Vnp1 + 6.0*(s-1.0)*s*(Xn-Xnp1)/h;
}

        """

        for var in self.vars:
            includes += tempcode%{'var': var, 'vartype': self.types[var]}
        includes += self.supportcode

        #include source for Gaussian white noise:
        gwn = False
        for line in self.noise.values():
            if 'gwn' in line:
                gwn = True
        if gwn:
            includes += gwn_code

        self.includes = includes
        code += '\n// This hash is included to detect a change in\n'+\
                '// the support code and recompile in this case\n'+\
                '// %s\n'%hashlib.md5(includes).hexdigest()
        self.code = code

    def output_ccode(self):
        txt = ''
        txt += self.includes
        txt += "\n\n"
        txt += self.code
        return txt

    def run(self):
        """run the simulation"""
        # get the parameters for the simulation into the local variables
        # so that they can be given to the c extension
        if not self.__hist_set: 
            sys.stderr.write("""Warning: The history has not been set at all or not been set 
            with the functions hist_from_arrays or hist_from_funcs. Please don't access the 
            history arrays directly.\n""")

        nstart = self._nstart

        for prm in self.params:
            exec('PAR%s=%s'%(prm, self.params[prm]))

        for prm in self.simul:
            exec('%s=%s'%(prm, self.simul[prm]))
        
        for var in self.vars:
            exec('hist%s_ar  = self.hist["%s"]'%(var, var))
            exec('Vhist%s_ar = self.Vhist["%s"]'%(var, var))

        Thist_ar = self.hist['t']
        chunk = self.chunk

        discont = list(set(self.discont))
        discont.sort()

        discont = np.array(discont)
        NumOfDiscont = len(discont)
        RSEED = self.rseed

        self.sol = weave.inline(self.code,
                            ['PAR%s'%p for p in self.params] +\
                            ['hist%s_ar'%var for var in self.vars] +\
                            ['Vhist%s_ar'%var for var in self.vars] +\
                            ['Thist_ar', 'dtmin', 'dtmax',\
                            'dt0', 'tfinal', 'RelTol', 'discont', 'NumOfDiscont',\
                            'AbsTol', 'nstart', 'chunk', 'MaxIter', 'RSEED'],
                   support_code=self.includes,
                   verbose=0,
                   compiler = 'gcc')

        for var in self.vars:
            if self.types[var] == 'Complex':
                try:
                    self.spline_tck[var+'real'] = splrep(self.sol['t'], self.sol[var].real)
                    self.spline_tck[var+'imag'] = splrep(self.sol['t'], self.sol[var].imag)
                except:
                    sys.stderr.write('Error: Scipy could not calculate spline for variable %s.\n'%var)
                    sys.stderr.flush()
            else:
                try:
                    self.spline_tck[var] = splrep(self.sol['t'], self.sol[var])
                except:
                    sys.stderr.write('Error: Scipy could not calculate spline for variable %s.\n'%var)
                    sys.stderr.flush()

    def sol_spl(self, t): 
        """Sample the solutions at times `t`.
        
        `t`
            Array of time points on which to sample the solution.

        Returns a dictionary with the sampled arrays. The keys are the 
        variable names. The key ``'t'`` corresponds to the sampling times.
        """
        erg = {'t': t}
        for var in self.vars:
            if self.types[var] == 'Complex':
                erg[var] = splev(t, self.spline_tck[var+'real']) + 1.0j*splev(t, self.spline_tck[var+'imag']) 
            else:
                erg[var] = splev(t, self.spline_tck[var])
        return erg

    def sample(self, tstart=None, tfinal=None, dt=None):
        """Sample the solution with `dt` steps between `tstart` and `tfinal`.
        
        `tstart`, `tfinal`
            Start and end value of the interval to sample.
            If nothing is specified `tstart` is set to zero
            and `tfinal` is set to the simulation end time.

        `dt` 
            Sampling size used. If nothing is specified a reasonable 
            value is calculated.

        Returns a dictionary with the sampled arrays. The keys are the 
        variable names. The key ``'t'`` corresponds to the sampling times.
        """
        if tstart == None:
            tstart = 0
        if tfinal == None:
            tfinal = self.sol['t'].max()
        if dt == None:
            dt = 0.5*(self.sol['t'][-1]-self.sol['t'][0])/len(self.sol['t'])
        t = np.arange(tstart, tfinal, dt)
        erg = self.sol_spl(t)
        return erg

if __name__ == '__main__':
    eqns = {'y1': '- y1 * y2(t-1.0) + y2(t-tau2)',
            'y2': 'y1 * y2(t-tau1) - y2',
            'y3': 'y2 - y2(t-tau2)'}

    params = {'tau1': 1.0, 'tau2': 10.0}

    dde = dde23(eqns=eqns, params=params, supportcode='')

    dde.set_sim_params(tfinal=100, AbsTol=10**-6, RelTol=10**-3)

    initfuncs = {
        'y1': lambda x: 5.0,
        'y2': lambda x: 0.1,
        'y3': lambda x: 1.1
        }

    dde.hist_from_funcs(initfuncs)
    dde.run()

    import pylab as pl

    t = dde.sol['t']
    y1 = dde.sol['y1']
    y2 = dde.sol['y2']
    y3 = dde.sol['y3']
    pl.plot(t,y1, 'rx')
    pl.plot(t,y2, 'gx')
    pl.plot(t,y3, 'bx')
    
    tspl=np.linspace(0,100, 10000)

    sol = dde.sample(dt=0.2)

    pl.plot(sol['t'], sol['y1'], 'ro')
    pl.plot(sol['t'], sol['y2'], 'go')
    pl.plot(sol['t'], sol['y3'], 'bo')

    pl.show()
