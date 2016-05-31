# Jozef Skakala, PML,  2016.

# The module contains the most basic functions of the multifractal scaling analysis. It is as concise as possible...


import numpy as np
import multifractal_basic_functions_pub as mbf




# This function computes the fluxes for a specific distribution given by the `field' variable. The function computes the fluxes at the scale given by the `step_size' variable. The fluxes are computed by a simple method of taking `delta field / mean(delta field)' and averaging this quantity through a circle originating at the point of the flux value. It is assumed that field = 0 means `masked' (land). The details of this computation are just a special case of the function scaling_increments(...).


def fluxes(field, step_size):

    step_size = int(step_size)
    region_height = np.shape(field)[0]
    region_width = np.shape(field)[1]
    flux = np.zeros((region_height, region_width))
    n_steps_lon = int(mbf.round_down(region_width/(step_size+0.0)))
    n_steps_lat = int(mbf.round_down(region_height/(step_size+0.0)))
    center = int((step_size-1)/2.0)
  
    for lon in range(0, n_steps_lon):
        for lat in range(0, n_steps_lat):

            if field[lat*step_size+center, lon*step_size+center] != 0:

                rel_case = 0

                for n in range(max((lon-1)*step_size+center,0),min((lon+1)*step_size+center+1,region_width), step_size):
                    for m in range(max((lat-1)*step_size+center,0),min((lat+1)*step_size+center+1,region_height), step_size): 

                        if field[m,n] != 0:

                            flux[lat*step_size : lat*step_size+2*center+1, lon*step_size : lon*step_size+2*center+1] += np.abs(field[m,n]-field[lat*step_size+center, lon*step_size+center])
                            rel_case+=1
              
                if rel_case > 1: 

                    flux[lat, lon]=flux[lat, lon]/(rel_case-1.0)
          
    return flux/np.mean(flux[flux != 0])



# This function typically calculates scaling of the fluxes (captured by the more general `field' variable!). It calculates the statistical moments scaling for the moments supplied (`momenta') from minimal (`scale_min') to maximal (`scale_max') scale, scales separated by scaling coefficient (`scale_coeff').  The moments are calculated in boxes with the area A = scale**2. The boxes however might be squashed (non-rectangular) by the anisotropy coefficient (`anisotropy'). It is typical to set anisotropy = 1.0, implying that the boxes are squares. The boundaries and the land lead in general to smaller effective box area than A = scale**2. The boxes with smaller effective (not necessarily geometric!) area are included in the analysis with a lower statistical weight (the weight is simply proportional to the box effective area). To resolve the assymetry of the analysis introduced by the regional boundaries, the boxes are defined symmetrically from all 4 corners of the rectangular region.
 
#This function is more general than just for the purpose of calculating statistical moments, it can calculate also mean variance per box, or mean standard deviation per box. It can also return the scale ratio at which the fluxes were computed. What is calculated is determined by the `output' variable with possible four values: output = (moment, variance, st_deviation, scale).


def scaling(field, momenta, scale_max, scale_min, scale_coeff, anisotropy, output):

# Defines main parameters used in the calculation.

    region_height = np.shape(field)[0]
    region_width = np.shape(field)[1]
    n_pixels = region_width*region_height
    mean_field_region = np.mean(field[field != 0])
    mean_het = []
    scale = np.array([])
    total_n_pixels_sea = len(field[field != 0])
    step = 1.0
    length = scale_min
    n_moments = len(momenta)

# This is the main while-loop running through all the scales. 

    while length < scale_max:    

        n_box_pixels = np.round(length**2.0)+0.0
        n_box_pixels_lat = np.round(np.sqrt(n_box_pixels/anisotropy))+0.0
        n_box_pixels_long = np.round(np.sqrt(n_box_pixels*anisotropy))+0.0
        eff_n_box_pixels = n_box_pixels_lat*n_box_pixels_long       # Deals with the discrete grid (lattice) structure where eff_n_box pixels != n_box_pixels

        het = np.zeros((n_moments))

 # If the effective (!) scale is larger than scale_min then it starts the analysis. 
       
        if eff_n_box_pixels >= scale_min**2.0:
    
            n_boxes_long = int(mbf.round_up(region_width/n_box_pixels_long))
            n_boxes_lat =  int(mbf.round_up(region_height/n_box_pixels_lat))
            
            if (n_boxes_long > 1) & (n_boxes_lat > 1):

# Loop through the individual boxes:

                for box_long in range(1, n_boxes_long+1):
                    for box_lat in range(1, n_boxes_lat+1): 
 
# The individual box has 4 different coordinates corresponding to the analysis that originated from 4 different regional vortices.
       
                        for iterate in range(1,5):
                            
                            if iterate == 1:
                                
                                coordinate_long_1 = (box_long-1)*n_box_pixels_long
                                coordinate_long_2 = min(box_long*n_box_pixels_long, region_width) 
                                coordinate_lat_1 = (box_lat-1)*n_box_pixels_lat
                                coordinate_lat_2 = min(box_lat*n_box_pixels_lat, region_height)
             
                            if iterate == 2:
              
                                coordinate_long_2 = region_width - (box_long-1)*n_box_pixels_long
                                coordinate_long_1 = max(region_width - box_long * n_box_pixels_long, 0)  
              
                            if iterate == 3:
            
                                coordinate_lat_2 = region_height - (box_lat-1)*n_box_pixels_lat              
                                coordinate_lat_1 = max(region_height - box_lat*n_box_pixels_lat, 0)
              
                            if iterate == 4: 
            
                                coordinate_long_1 = (box_long-1)*n_box_pixels_long
                                coordinate_long_2 = min(box_long*n_box_pixels_long, region_width)
            
             
                            field_box = field[coordinate_lat_1 : coordinate_lat_2, coordinate_long_1 : coordinate_long_2]
                            box_size_pixels = (coordinate_long_2 + 1 - coordinate_long_1)*(coordinate_lat_2 + 1 - coordinate_lat_1)
                            sea_n_box_pixels = len(field_box[field_box != 0])+0.0
                            box_importance = sea_n_box_pixels/(4*total_n_pixels_sea)


# This records the desired parameter of the analysis.

            
                            if (output == 'moment') & (sea_n_box_pixels > 0):
            
                                het += box_importance * (np.mean(field_box[field_box != 0])/mean_field_region)**momenta
              
                            if (output == 'variance') & (sea_n_box_pixels > 1):
  
                                het += box_importance * np.var(field_box[field_box != 0])
   
                            if (output == 'st_deviation') & (sea_n_box_pixels > 1): 
  
                                het += box_importance * np.sqrt(np.var(field_box[field_box != 0]))/mean_field_region
                            
                mean_het.append(het)
                scale = np.append(scale, np.sqrt(eff_n_box_pixels/n_pixels))

        length = scale_min*scale_coeff**step
        step+=1     
        
    if output == 'scale': 
        return scale 
    else: 
        return np.asarray(mean_het)



# This function computes the field increments scaling. It has the same structure as the previous function (for the details see the function scaling(...)) with the `output' variable having values output = (moment, scale). The increments are computed around the circle with the radius = scale across all the relevant points of the region. It is as always assumed that field = 0 means `masked', or in other words land. The scale is here returned with values in grid pixels, rather than in values of the maximal scale. This is a difference to the previous scaling(...) function.

def scaling_increments(field, momenta, scale_max, scale_min, scale_coeff, output): 

    region_height = np.shape(field)[0]
    region_width = np.shape(field)[1]
    scale = np.array([])
    delta_field = []
    step=1.0
    length = scale_min

    while length < scale_max:
  
        cases=0.0
        delta=np.zeros((len(momenta)))

# This loop goes through individual pixels (/points).
        
        for pixel_long in range(0, region_width):
            for pixel_lat in range(0, region_height): 

# This inner loop goes through the circle with radius = scale and the center at the individual point of the outer loop.
        
                for substep in range(0,int(step)):
                    for sign_1 in range(-1,2,2):
                        for sign_2 in range(-1,2,2):

                            n_x=scale_min*scale_coeff**substep
                            n_y=np.round(np.sqrt(length**2-n_x**2))
                            n_x=np.round(n_x)
                            coordinate_long = pixel_long + n_x*(-1)**sign_1
                            coordinate_lat = pixel_lat + n_y*(-1)**sign_2
          
                            if (region_width > coordinate_long >= 0) & (region_height > coordinate_lat >= 0):
            
                                if (field[pixel_lat, pixel_long] != 0) & (field[coordinate_lat, coordinate_long] != 0):                              
                         
                                    delta += np.abs(field[pixel_lat, pixel_long]-field[coordinate_lat, coordinate_long])**momenta
                                    cases+=1            
          
        if cases > 10: 
            delta_field.append(delta/cases)
            scale=np.append(scale, length)

        length = scale_min*scale_coeff**step
        step+=1    
        
                
    if output == 'moment': 
        return np.asarray(delta_field) 
    else: 
        return scale[0:len(delta_field)]



# This is a function that calculates the universal multifractal (2-parametric) fit of the moments scaling function K. It returns the two parameters (alpha and C_1) and also the error of the fit. The names of the output variables are self-explanatory. The input is the K function to be fitted and also the relevant range of moments of the fit (must be consistent with K!). Since the range of relevant values of the C_1 ('C') parameter is largely case dependent, the fitting function allows to take the estimated relevant range as an argument. Also to keep the speed of the analysis optimal it allows the user to define the C fit resolution through the `C_step' variable. This is set automatically for the alpha parameter. The fit is based on minimizing the relative error.


def UM_fit(K, momenta, C_min, C_max, C_step):

    C_fit = C_min
    alpha_fit = 0.005
    error = np.mean(np.abs(((C_fit/(alpha_fit-1.0))*(momenta[momenta != 1]**alpha_fit - momenta[momenta != 1])-K[momenta != 1])/K[momenta != 1]))
    n_C_elements = mbf.round_up((C_max-C_min)/C_step)
    C_elements = np.linspace(C_min, C_max, num=n_C_elements)
    alpha_elements = np.linspace(0.01, 2.0, num=400)
    alpha_elements = alpha_elements[alpha_elements != 1]

    for C in C_elements:
        for alpha in alpha_elements:  
        
            if np.mean(np.abs(((C/(alpha-1.0))*(momenta[momenta != 1]**alpha - momenta[momenta != 1]) - K[momenta != 1] )/ K[momenta != 1]))  < error:
                error = np.mean(np.abs(((C/(alpha-1.0))*(momenta[momenta != 1]**alpha - momenta[momenta != 1])- K[momenta != 1] ) / K[momenta != 1]))
                alpha_fit = alpha
                C_fit = C
        
    return [alpha_fit, C_fit, error]



# This function calculates the 3 scaling parameters H, C_1, alpha and also the remaining scaling information. The moment scaling function is computed through standard linear interpolation. The alpha, C_1 parameters are computed through the fit provided by UM_fit(...) function. The remaining scaling information is the outer scale of the process `outer_scale' parameter and also the fluctuation characteristic sizes. The output variable has two options: output=(K,parameters). If 'K' is selected then in returns the moment scaling function, otherwise it returns the parameters. The parameter output of the function is: [H, alpha, C_1, outer_scale, fluctuation size, the error of the UM fit]. There are two separate moments: moments for increments (momenta_inc) and for fluxes (momenta_flux). This small complication is due to the fact that due to numerical computational reasons it is sometimes convenient to select the moments for fluxes and increments slightly differently. Note for example module defining the parameters for the analysis.


def UM_parameters(flux_scaling, inc_scaling, scales_flux, scales_inc, momenta_flux, momenta_inc, output):
 
    scales_flux = scales_flux[np.argwhere(flux_scaling[:,10]>1)][:,0]
    flux_scaling = flux_scaling[0:len(scales_flux),:]
    H =  np.cov(np.log(inc_scaling[: , np.argwhere(momenta_inc == 1)[0][0]]), np.log(scales_inc))[0][1]/np.var(np.log(scales_inc))
    a_inc = np.mean(np.log(inc_scaling[: , np.argwhere(momenta_inc == 1)[0][0]])) - H*np.mean(np.log(scales_inc))        
    n_moments = len(momenta_flux)
    K = np.zeros((n_moments))
    a_flux = np.zeros((n_moments))   
    
    for moment in range(0,n_moments):
        K[moment] = np.cov(-np.log(flux_scaling[: , moment]), np.log(scales_flux))[0][1]/np.var(np.log(scales_flux))
        a_flux[moment] = np.mean(np.log(flux_scaling[: , moment])) + K[moment]*np.mean(np.log(scales_flux))

    outer_scale = np.exp(a_flux[9]/K[9])
    sc_flux_parameters = UM_fit(K, momenta_flux, 0, 2.0, 0.001)

    if output == 'K': 
        return K
    else:
        return np.array([H, sc_flux_parameters[0], sc_flux_parameters[1], outer_scale, np.exp(a_inc), sc_flux_parameters[2]])
    






