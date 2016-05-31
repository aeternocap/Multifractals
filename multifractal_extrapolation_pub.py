# Author: Jozef Skakala, PML, 2016 
# This is the core function for the extrapolation. Plug in field at larger scales ('field') and obtain returned field at lower scales (determined by the iteraion exponent: `n_iterations').

import numpy as np
import multifractal_basic_functions_pub as mbf
import multifractal_extrapolate_functions_pub as mef
from multifractal_class_pub import multifractals
import multifractal_parameter_values_pub as mpv

def field_extrapolation(field, n_iterations):

    field_multifractal = multifractals(field)    # Calculate the multifractal scaling of the field
    parameters = field_multifractal.UM_parameters()  # Extract the UM parameters
    fluxes = field_multifractal.fluxes()  # Extract the fluxes
    ratio_bound = mpv.ratio_bound() # Upper bound on the fraction of the space where fluctuations are re-distributed
    factor = parameters[4]*(1/2.0**n_iterations)**parameters[0]   # the factor that relates fluctuations to fluxes
    

    fluxes_extrapolated = mef.fluxes_extrapolate(fluxes, n_iterations, parameters)   # Extrapolate fluxes

    field_scaled_down = mef.smoothen(mef.lower_resolution(field, 2**n_iterations), 2**n_iterations, 2**n_iterations, 2.0)  # Get the smoothen version of the field on the lower scale

    field_extrapolated = mef.fluctuations_distribute(field_scaled_down, fluxes_extrapolated, factor, ratio_bound)  # Distribute the field fluctuations on the lower scale using the 1. fluxes and 2. smoothened field .

    return field_extrapolated
