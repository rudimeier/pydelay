# Copyright (C) 2009 Valentin Flunkert

import re
import numbers
import numpy as np

def gen_disconts(t0, t1, delays, initdisconts=None, order=3, rounddigits=5):
    """Generate a list of all possible discontinuities in the range [t0, t1]
    up to 'order' with initial discontinuities given by 'initdisconts'.
    The adaptive step size methods should step on the discontinuity points.

    'delays' can be a dictionary or a list of numbers.

    >>> [0, 2, 4, 6, 8, 10, 12] == gen_disconts(0, 100, [2.0, 4.0], [-10, 0], 3)
    True

    >>> [0, 3, 5, 6, 8, 10] == gen_disconts(0, 100, [3.0, 5.0], [0], 2)
    True
    """
    if initdisconts == None:
        initdisconts = [0]
    if isinstance(delays, dict):
        delays = delays.values()
    newdis = initdisconts
    alldis = []
    for o in range(order):
        alldis = alldis + newdis
        tempdis = newdis
        newdis = []
        for dis in tempdis:
            for delay in delays:
                newdis.append(dis + delay)
    alldis = alldis + newdis
    alldis = [round(dis,rounddigits) for dis in alldis if t0 <= dis and dis <= t1]
    erg = list(set(alldis))
    erg.sort()
    return erg

def _symbols_allowedQ(eqn):
    """Takes an eqn string and returns True if all symbols
    are allowed and False otherwise
    
    >>> map(_symbols_allowedQ, ['2.0*sin(x)* 3/2 + 1-2', 'a**b', 'a^b', 'a{b', 'a}b'])
    [True, False, False, False, False]
    >>> map(_symbols_allowedQ, ['a[b', 'a]b', 'a&c', 'a#b', 'a!b', 'a@b', 'a$b', 'a%b'])
    [False, False, False, False, False, False, False, False]
    >>> map(_symbols_allowedQ, ['a?b', 'a;b', 'a>b', '|', "\\\\"])
    [False, False, False, False, False]
    """
    ForbiddenPattern =\
            r'(\*\*|\^|\{|\}|\[|\]|\&|\#|\!|\@|\$|\%|\?|\;|\>|\<|\||\\)'
    res = re.search(ForbiddenPattern, eqn)
    if res:
        return False
    else:
        return True

def _parenthesis_balancedQ(eqn):
    """Takes an eqn string and return True if parenthesis
    are balanced and False otherwise

    >>> map(_parenthesis_balancedQ,['(fdjd)*d((2-1)+x*2)', 'fs*(1-(x*2*(a+b))', 'dfs * (x-2) + b)'])
    [True, False, False]
    """
    # 'opened_parenthesis' is increased, with '(' and decreased with ')'
    # and should be zero at the end (and never negative)
    opened_parenthesis = 0
    for ch in eqn:
        if ch == '(':
            opened_parenthesis += 1
        elif ch in ')':
            if opened_parenthesis == 0:
                return False
            else:
                opened_parenthesis -= 1
    if opened_parenthesis == 0:
        return True
    else:
        return False

def assure_array(array_candidate):
    """If array_candidate is a tuple, list or array cast it to a numpy array 
    and return it, else assert an error"""
    assert isinstance(array_candidate, (tuple, list, np.ndarray)),\
        "Error: this should be a list, tuple or array"
    return np.array(array_candidate)

def isrealnum(real_candidate):
    """
    >>> map(isrealnum, [0.1, 0.2, 0.1+0.2j, "a", (1.0, 2.0)])
    [True, True, False, False, False]
    """
    if not isinstance(real_candidate, numbers.Number):
        return False
    if not hasattr(real_candidate, 'imag'):
        return True
    return real_candidate.imag == 0

def isnum(num_candidate):
    """
    >>> map(isnum, [0.1, 0.2, 0.1+0.2j, "a", (1.0, 2.0)])
    [True, True, True, False, False]
    """
    return isinstance(num_candidate, numbers.Number)

gwn_code = """
double ra01()
{  
    return(double(rand())/RAND_MAX);
}

double gwn()
{
    static int iset=0;
    static double gset;

    double fac,rsq,v1,v2;

    if(iset==0) 
    {
        do 
        {
            v1=2.0*ra01()-1.;
            v2=2.0*ra01()-1.;
            rsq=v1*v1+v2*v2;
        } 
        while(rsq>1.0 || rsq==0);

        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } 
    else 
    {
        iset=0;
        return gset;
    }
}
"""

if __name__== '__main__':
    import doctest
    doctest.testmod()
