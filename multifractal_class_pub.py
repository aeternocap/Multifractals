
# Jozef Skakala, PML, 2016.

# Here is the class that has a distribution as an input and returns the complete multifractal information about the distribution. It computes the UM scaling using all the functions defined in the `multifractal_scaling_essential_pub' module. The names of the attributes are self-explanatory, perhaps with the exception of 'K', which is the standard notation for the moment scaling function.



import numpy as np
import multifractal_scaling_essential_pub as mse
import multifractal_parameter_values_pub as pa




class multifractals:

    def __init__(self, field):
        self.field = field
        self.description = "Field with multifractal properties - the subject to the analysis."
        self.author = "JS"

        self.momenta_flux = pa.flux_momenta()
        self.momenta_inc = pa.inc_momenta()
        self.field_inc_scaling = mse.scaling_increments(self.field, self.momenta_inc, pa.max_scale(), pa.min_scale(), pa.scale_coeff(), 'moment')
        self.scales_inc = mse.scaling_increments(self.field, self.momenta_inc, pa.max_scale(), pa.min_scale(), pa.scale_coeff(), 'scale')          
        self.flux = mse.fluxes(self.field, 1.0)
        self.flux_scaling = mse.scaling(self.flux, self.momenta_flux, pa.max_scale(), pa.min_scale(), pa.scale_coeff(), 1.0, 'moment') 
        self.scales_flux = mse.scaling(self.flux, self.momenta_flux, pa.max_scale(), pa.min_scale(), pa.scale_coeff(), 1.0, 'scale') 
        self.K = mse.UM_parameters(self.flux_scaling, self.field_inc_scaling, self.scales_flux, self.scales_inc, self.momenta_flux, self.momenta_inc, 'K')
        self.parameters = mse.UM_parameters(self.flux_scaling, self.field_inc_scaling, self.scales_flux, self.scales_inc, self.momenta_flux, self.momenta_inc, 'parameters')

# The UM_parameters contains: [H, alpha, C_1, outer scale of process, fluctuations proportionality constant, UM fit error]. Please note that the outer_scale calculated through the UM_parameters function is in the units of the regional scale.

        

    def fluxes(self):
        return self.flux

    def increments_scaling(self):
        return self.field_inc_scaling

    def fluxes_scaling(self):
        return self.flux_scaling

    def UM_parameters(self):
        return self.parameters

    def scales_increments(self):    
        return self.scales_inc

    def scales_fluxes(self):
        return self.scales_flux     

    def moment_scaling_function(self):
        return self.K
