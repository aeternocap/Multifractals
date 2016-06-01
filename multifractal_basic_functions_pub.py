#Jozef Skakala, PML, 2016.

# In this module are some quite trivial functions that can be used for the purpose of multifractal analysis. The number of the functions can be easily extended, on the other hand there already are some functions that might not be of much use..  

import numpy as np
import math
from scipy.integrate import quad
import multifractal_parameter_values_pub as mpv

# This rounds up / down float variables (haven't found this between basic functions, strange...).

def round_up(x): 
    if round(x)<x:
        return round(x)+1 
    else: 
        return round(x)

def round_down(x): 
    if round(x)<=x:
        return round(x)
    else:
        return round(x)-1



def signature(x):
    if x<0: 
        return -1 
    else: 
        return 1


def Heaviside(x):
    if x<0:
        return 0 
    else: 
        return 1


def st_gauss(x):
    return 2*np.exp(-x**2/2.0)/sqrt(2*3.1416)


# Simple algorithm to find the norm of function

def norm(function, argument):

    n_steps = len(argument)
    norm = 0.0    

    for step in range(0, n_steps-1):
   
        delta = argument[step+1]-argument[step]
        norm += delta*(function[step]+function[step+1])/2

    return norm


# This is used in the extrapolation of the fluxes (relating uniform PDF for which the generator is present to the real UM PDF).

def integrate_up_to_value(function, argument, value):

    size_function = len(function)
    index = 0
    integral = 0

    while ((integral < value) & (index < size_function)):
      
        step_size = argument[index+1]-argument[index]
        integral += ((function[index+1]+function[index])/2.0)*step_size
        index += 1

    return argument[index]


# These two are used to calculate inverse Mellin transform to obtain the UM PDF. Haven't found anything better...


def inverse_mellin_UM_integrand(z, x_o, C, alpha, scale):

    q = 1.5 + z*1j

    function = x_o**(-q)*scale**(C/(alpha-1.0)*((q-1)**alpha-q+1))/(2*np.pi)

    return function.real


def inverse_mellin_UM(parameters, scale):

    H = parameters[0]
    C = parameters[2]
    alpha = parameters[1]

    x = mpv.PDF_argument()
    function_transformed = np.zeros((len(x)))
    index = 0

    for x_o in x:

        function_transformed[index] = quad(inverse_mellin_UM_integrand, -np.inf, np.inf, args=(x_o, C, alpha, scale))[0]
        index +=1 

   # function_transformed = function_transformed/norm(function_transformed, x)  # Although PDF is normalized by definition, this is to correct numerical errors.

    return function_transformed

    
