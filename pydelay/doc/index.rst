.. pydelay documentation master file, created by
   sphinx-quickstart on Mon Nov  2 10:52:40 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================================================
`pydelay <http://pydelay.sourceforge.net/>`_ v |version|
========================================================

.. toctree::
   :maxdepth: 2

Introduction
============
pydelay is a program which translates a system of delay differential equations
(DDEs) into simulation C-code and compiles and runs the code (using scipy
weave). This way it is easy to quickly implement a system of DDEs but you still
have the speed of C. The Homepage can be found here:

        http://pydelay.sourceforge.net/ 

It is largely inspired by `PyDSTool <http://www.cam.cornell.edu/~rclewley/cgi-bin/moin.cgi/>`_.

The algorithm used is based on the Bogacki-Shampine method [1]_
which is also implemented in `Matlab's dde23` [2]_.

We also want to mention `PyDDE <http://users.ox.ac.uk/~clme1073/python/PyDDE/>`_
-- a different python program for solving DDEs.

**License**

pydelay is licensed under the MIT License.

.. htmlonly::
    **Citing**

    You can cite pydelay through the `arXiv <http://arxiv.org>`_.::
    
        @unpublished{FLU09a,
            title = {pydelay -- a python tool for solving delay differential equations},
            author = {Flunkert, V. and Sch{\"o}ll, E. },
            year = {2009},
            eprint={nlin.CD/0911.1633}
        }

Installation and requirements
-----------------------------
**Unix:**

You need `python <http://www.python.org/>`_, `numpy and scipy <http://numpy.scipy.org/>`_ and the `gcc`-compiler.
To plot the solutions and run the examples you also need `matplotlib <http://matplotlib.sourceforge.net/index.html>`_.

To install pydelay grab the latest `tar.gz` from the website and install the package in the usual way::

        cd pydelay-$version
        python setup.py install

When the package is installed, you can get some info about the functions and the usage with::

        pydoc pydelay


**Windows:**

The solver has not been tested on a `windows` machine. It could perhaps work
under `cygwin <http://www.cygwin.com/>`_.

An example 
----------

The following example shows the basic usage. It solves the Mackey-Glass
equations [3]_ for initial conditions which lead to a periodic orbit (see [4]_
for this example).

.. plot:: pyplots/mackey-glass.py
     :include-source: 

Usage
=====

Defining the equations, delays and parameters
---------------------------------------------

Equations are defined using a python dictionary.
The keys are the variable names and the entry is the right hand side of the 
differential equation.
The string defining the equation has to be a valid C expression, i.e.,
use ``pow(a,b)`` instead of ``a**b`` etc.

Delays are written as ``(t-delay)``, where ``delay`` can be 
some expression involving parameters and numbers but not (yet) involving 
the time ``t`` or the dynamic variables::

    eqns = {
        'y1': '- y1 * y2(t-tau) + y2(t-1.0)',
        'y2': 'a * y1 * y2(t-2*tau) - y2',
        'y3': 'y2 - y2(t-(tau+1))'
      }


Complex variables can be defined by adding ``':c'`` or
``':C'`` in the eqn-dictionary.
The imaginary unit can be used through ``'ii'`` in the equations::

    eqns = {
        'z:c': '(la + ii*w0 + g*pow(abs(z),2) )*z + b*(z(t-tau) - z(t))',
    }

Parameters are defined in a separate dictionary where the keys are 
the parameter names, i.e.,::

    params = {
        'a'  : 0.2,
        'tau': 1.0
    }

Setting the history
-------------------

The history of the variables is stored in the dictionary ``dde23.hist``.
The keys are the variable names and there is an additional key ``'t'`` for
the time array of the history. 

There is a second dictionary ``dde23.Vhist``
where the time derivatives of the history is stored (this is needed for the
solver). When the solver is initialized, i.e.,::

        dde = dde23(eqns, params)

the history of all variables (defined in ``eqns``) is initialized to an
array of length ``nn=101`` filled with zeros. The time array is evenly
spaced in the interval ``[-maxdelay, 0]``.

It is possible to manipulate these arrays directly, however this is not
recommended since one easily ends up with an ill-defined history resulting for
example in segfaults or false results.

Instead use the following methods to set the history.

.. automethod:: pydelay.dde23.hist_from_funcs

.. automethod:: pydelay.dde23.hist_from_arrays

Note that the previously used methods ``hist_from_dict``, ``hist_from_array``
and ``hist_from_func`` (the last two without ``s``) have been removed, since it
was too easy to make mistakes with them.

The solution
------------
After the solver has run, the solution (including the history) is stored 
in the dictionary ``dde23.sol``. The keys are again the variable names 
and the time ``'t'``. Since the solver uses an adaptive step size method,
the solution is not sampled at regular times. 

To sample the solutions at regular (or other custom spaced) times there 
are two functions.

.. automethod:: pydelay.dde23.sample

.. automethod:: pydelay.dde23.sol_spl

These functions use a cubic spline interpolation of the solution data.

Noise
-----
Noise can be included in the simulations. Note however, that the method used is
quite crude (an Euler method will be added which is better suited for noise
dominated dynamics). The deterministic terms are calculated with the usual
Runge-Kutta method and then the noise term is added with the proper scaling of
``$\sqrt{dt}$`` at the final step. To get accurate results one should use small
time steps, i.e., ``dtmax`` should be set small enough.

The noise is defined in a separate dictionary. The function ``gwn()`` can
be accessed in the noise string and is a Gaussian white noise term of unit
variance. The following code specifies an Ornstein-Uhlenbeck process.::

    eqns = { 'x': '-x' }
    noise = { 'x': 'D * gwn()'}
    params = { 'D': 0.00001 }

    dde = dde23(eqns=eqns, params=params, noise=noise)

You can also use noise terms of other forms by specifying an appropriate
C-function (see the section on custom C-code).

Custom C-code
-----------------
You can access custom C-functions in your equations by 
adding the definition as ``supportcode`` for the solver. 
In the following example a function ``f(w,t)`` is defined through 
C-code and accessed in the eqn string.::

    # define the eqn f is the C-function defined below
    eqns = { 'x': '- x + k*x(t-tau) + A*f(w,t)' }
    params = { 
        'k'  : 0.1, 
        'w'  : 2.0, 
        'A'  : 0.5, 
        'tau': 10.0
    }

    mycode = """
    double f(double t, double w) {
        return sin(w * t);
    }
    """

    dde = dde23(eqns=eqns, params=params, supportcode=mycode)

When defining custom code you have to be careful with the types.
The type of complex variables in the C-code is ``Complex``.
Note in the above example that ``w`` has to be given as an input to 
the function, because the parameters can only be accessed from the eqns string 
and not inside the supportcode. (Should this be changed?)

Using custom C-code is often useful for switching terms on and off.
For example the Heaviside function may be defined and used as follows.::

    # define the eqn f is the C-function defined below
    eqns = { 'z:c': '(la+ii*w)*z - Heavi(t-t0)* K*(z-z(t-tau))' }
    params = { 
        'K'  : 0.1 , 
        'w'  : 1.0,
        'la' : 0.1, 
        'tau': pi, 
        't0' : 2*pi 
    }

    mycode = """
    double Heavi(double t) {
        if(t>=0)
            return 1.0;
        else 
            return 0.0;
    }
    """
    dde = dde23(eqns=eqns, params=params, supportcode=mycode)

This code would switch a control term on when ``t>t0``.
Note that ``Heavi(t-t0)`` does not get translated to a delay 
term, because ``Heavi`` is not a system variable.

Since this scenario occurs so frequent the Heaviside function (as defined above)
is included by default in the source code.

Use and modify generated code
-----------------------------
The compilation of the generated code is done with ``scipy.weave``. 
Instead of using weave to run the code you can directly access the generated
code via the function ``dde23.output_ccode()``. This function returns 
the generated code as a string which you can then store in a source file.

To run the generated code manually you have to set the precompiler flag\\
``#define MANUAL`` (uncomment the line in the source file) to exclude the
python related parts and include some other parts making the code a valid stand
alone source file. After this the code should compile with 
``g++ -lm -o prog source.cpp`` and you can run the program manually.

You can specify the history of all variables in the source file
by setting the ``for`` loops after the comment\\
``/* set the history here ...  */``.

Running the code manually can help you debug, if some problem 
occurs and also allows you to extend the code easily.

Another example
---------------

The following example shows some of the things discussed above. 
The code simulates the Lang-Kobayashi laser equations [5]_

.. math::

        E'(t) &= \frac{1}{2}(1+i\alpha) n E + K E(t-\tau)\\
        T n'(t) &= p - n - (1+n) | E|^2

.. plot:: pyplots/lk.py
     :include-source: 

Module Reference
================

.. automethod:: pydelay.dde23.__init__
.. automethod:: pydelay.dde23.set_sim_params
.. automethod:: pydelay.dde23.hist_from_arrays
.. automethod:: pydelay.dde23.hist_from_funcs
.. automethod:: pydelay.dde23.output_ccode
.. automethod:: pydelay.dde23.run
.. autoclass:: pydelay.dde23
        
.. :members: __init__, set_sim_params, hist_from_arrays, hist_from_funcs, output_ccode, run 

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. [1] Bogacki, P. and Shampine, L. F., A 3(2) pair of Runge - Kutta formulas, Applied Mathematics Letters 2, 4, 321 ISSN 0893-9659, (1989).
.. [2] Shampine, L. F. and Thompson, S., Solving DDEs in Matlab, Appl. Num. Math. 37, 4, 441 ( 2001)
.. [3] Mackey, M. C. and Glass, L. (1977). Pathological physiological conditions resulting from instabilities in physiological control system. Science, 197(4300):287-289.
.. [4] http://www.scholarpedia.org/article/Mackey-Glass_equation  
.. [5] Lang, R.  and Kobayashi, K. , External optical feedback effects on semiconductor injection laser properties, IEEE J. Quantum Electron. 16, 347 (1980)
