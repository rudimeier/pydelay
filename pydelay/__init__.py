"""
pydelay, A tool for solving delay differential equations.
Copyright (C) 2009  Valentin Flunkert

Author: Valentin Flunkert <flunkert@gmail.com>

Last update: 23.10.2009
"""
from _dde23 import dde23 
#from constantStepper import dde3
__all__ = ['dde23', 'gen_disconts']
__version__ = '0.1.1'
