# Jozef Skakala, PML, 2016.

# This short module contains some generic parameter values.


import numpy as np


def inc_momenta():
    return np.arange(0.2, 3.2, 0.2)

def flux_momenta():
    return np.arange(0.05, 3.25, 0.2)

def max_scale():
    return 90.0

def min_scale():
    return 3.0

def scale_coeff():
    return 1.1

def ratio_bound():
    return 0.99

def PDF_argument():
    return np.exp(np.linspace(-5,5, num=1000))
    
