# Author: Jozef Skakala,   PML, 2016
# This module contains all the necessary functions used in the multifractal data extrapolation. 



import numpy as np
import multifractal_basic_functions_pub as mbf
import random as rn
from copy import deepcopy
import multifractal_parameter_values_pub as mpv



# This function lowers the resolution of a distribution, by a simple identity extrapolation. The code is self-explanatory.

def lower_resolution(field, level):
 
    region_height = np.shape(field)[0]
    region_width = np.shape(field)[1]
    field_low_resolution = np.zeros((level*region_height, level*region_width))
 
    for pixel_lon in range(0, region_width):
        for pixel_lat in range(0, region_height): 
     
            field_low_resolution[pixel_lat*level:(pixel_lat+1)*level, pixel_lon*level:(pixel_lon+1)*level] = field[pixel_lat, pixel_lon]
      
    return field_low_resolution
 


# This function provides a smoothening algorithm. (Haven't found anything equally suitable.) The algirthms is based on obtaining a smoothened field using the `level_1'-scale averaging. The field is also set to recover the correct mean values at the `level_2' scale.

def smoothen(field, level_1, level_2, power):

    step = (level_1-1)/2.0  # this is initial smoothening scale
    region_height = np.shape(field)[0]
    region_width = np.shape(field)[1]
    field_smooth = np.zeros((region_height, region_width))   # this is soothened field
 
    for iterate in range(1, 6):   # there are 5 steps in which the smoothening procedure will be repeated. In each step the averaging and stretching approaches the ideally smoothened distribution more.

        for subpixel_lon in range(0, region_width):
            for subpixel_lat in range(0, region_height): 
     
                if field[subpixel_lat, subpixel_lon] != 0:    # if the sub-pixel is sea then go on
       
                    field_box = field[max(0, subpixel_lat-step):min(region_height,subpixel_lat+step+1), max(0, subpixel_lon-step):min(region_width,subpixel_lon+step+1)]
                    field_smooth[subpixel_lat, subpixel_lon] = np.mean(field_box[field_box != 0])     # This does the smoothening via supplying the averaged value at the level_1 scale.
     
# This next step tries to modify the smoothened values in order to recover the correct mean values at the scale `level_2'. It does so by stretching the differences between field values to obtain the right balance. The level of stretching is provided by the `power' value.

        if power != 0:
            for subpixel_lon in range(0, int(mbf.round_up(region_width/(level_2+0.0)))):
                for subpixel_lat in range(0, int(mbf.round_up(region_height/(level_2+0.0)))): 
       
                    if np.mean(field[subpixel_lat*level_2 : min((subpixel_lat+1)*level_2, region_height), subpixel_lon*level_2 : min((subpixel_lon+1)*level_2, region_width)]) != 0:    # if the pixel is on sea, go on
       
                        # This records the smoothened and the ``raw'' value to compare

                        field_smooth_box = field_smooth[subpixel_lat*level_2 : min((subpixel_lat+1)*level_2, region_height), subpixel_lon*level_2 : min((subpixel_lon+1)*level_2, region_width)]
                        field_box =  field[subpixel_lat*level_2 : min((subpixel_lat+1)*level_2, region_height), subpixel_lon*level_2 : min((subpixel_lon+1)*level_2, region_width)]
                        mean_smooth = np.mean(field_smooth_box[field_smooth_box != 0])
                        mean_field = np.mean(field_box[field_smooth_box != 0])
    
                       # if the values are not equal, this algorithm stretches the smoothened value to get the match. The two cases where the smoothened value is smaller / bigger than the desired value are dealt with separately.

              # Case 1

                        if mean_smooth < mean_field:
         
                            delta = np.mean(np.abs(field_smooth_box[field_smooth_box != 0] - min(field_smooth_box[field_smooth_box != 0]))**power)     
         
                            if delta != 0:
                                C = (mean_field - min(field_smooth_box[field_smooth_box != 0]))/delta
                            else:
                                C = 0
     
                            field_smooth_box[field_smooth_box != 0]= min(field_smooth_box[field_smooth_box != 0]) + C*(field_smooth_box[field_smooth_box != 0]-min(field_smooth_box[field_smooth_box != 0]))**power   # C and power provide the stretching
     
               # Case 2
    
                        if mean_smooth > mean_field:
     
                            delta = np.mean(np.abs(field_smooth_box[field_smooth_box != 0] - max(field_smooth_box[field_smooth_box != 0]))**power)    
    
                            if delta != 0: 
                                C = (max(field_smooth_box[field_smooth_box != 0])-mean_field)/delta 
                            else: 
                                C = 0
      
                            field_smooth_box[field_smooth_box != 0]= max(field_smooth_box[field_smooth_box != 0]) - C*(max(field_smooth_box[field_smooth_box != 0])-field_smooth_box[field_smooth_box != 0])**power
  
     
                        field_smooth[subpixel_lat*level_2: min((subpixel_lat+1)*level_2, region_height), subpixel_lon*level_2: min((subpixel_lon+1)*level_2, region_width)] = field_smooth_box
     
    
        field = field_smooth
     
    return field_smooth
 



# This function redistributes the coordinates on the extrapolated smaller grid. It returns the coordinates of the smaller grid cells. `Info' variable contains information about whether longitudes, or latitudes are being extrapolated.

def extrapolate_coord(coordinates, n_iterations, info):
    region_size = np.shape(coordinates)
    region_height = region_size[0]
    region_width = region_size[1]
    coordinates_extrapolated = np.zeros((region_height*2**n_iterations, region_width*2**n_iterations))
     
    for pix_lon in range(0, int(region_width*2**n_iterations)):
        for pix_lat in range(0, int(region_height*2**n_iterations)):
    
            if info == 'longitudes':
                coord_max = coordinates[min(mbf.round_up(pix_lat/(2.0**n_iterations)), region_height-1), min(round(pix_lon/(2.0**n_iterations)), region_width-1)]
                coord_min = coordinates[min(mbf.round_down(pix_lat/(2.0**n_iterations)), region_height-1)][min(round(pix_lon/(2.0**n_iterations)), region_width-1)]
                coordinates_extrapolated[pix_lat, pix_lon] = coord_min + (pix_lon - mbf.round_down(pix_lon/(2.0**n_iterations)))*(coord_max-coord_min)/(2.0**n_iterations)            
            
            if info == 'latitudes':
                coord_max = coordinates[min(round(pix_lat/(2**n_iterations)), region_height-1)][min(mbf.round_up(pix_lon/(2**n_iterations)), region_width-1)]
                coord_min = coordinates[min(round(pix_lat/(2**n_iterations)), region_height-1)][min(mbf.round_down(pix_lon/(2**n_iterations)), region_width-1)]
                coordinates_extrapolated[pix_lat, pix_lon] = coord_min + (pix_lat - mbf.round_down(pix_lat/(2.0**n_iterations)))*(coord_max-coord_min)/(2.0**n_iterations)
      
    return coordinates_extrapolated  




# This stochastically redistributes fluxes at a lower scale using the universal multifractal model.


def fluxes_extrapolate(flux, n_iterations, UM_parameters):

 # This reads all the multifractal information

    H = UM_parameters[0]   
    C = UM_parameters[2]
    alpha = UM_parameters[1]
    R = UM_parameters[3]

 # Define the R/lowest_scale ratio.

    scale = R*2**n_iterations

 # Obtain the real PDF from the inverse Mellin transform.

    PDF = mbf.inverse_mellin_UM(UM_parameters, scale)

 # The basic regional parameters.

    region_width = int(np.shape(flux)[0]*2**n_iterations)
    region_height = int(np.shape(flux)[1]*2**n_iterations)

  # The outcome of the extrapolation.

    flux_extrapolated = np.zeros((region_height, region_width))

  # The loop that runs through individual pixels on the extrapolated grid.

    for lat in range(0, region_height):

        for lon in range(0, region_width):
    
            value = rn.random()   # Generate random value between 0 and 1 (with uniform distribution).
            factor = mbf.integrate_up_to_value(PDF, mpv.PDF_argument(), value)  # obtain from the randomly generated value the extrapoaltion factor.
            flux_extrapolated[lat, lon] = factor*flux[lat/2**n_iterations, lon/2**n_iterations]  # Determine the flux at the lower scales.

    return flux_extrapolated
    



# This function stochastically redistributes fluctuations corresponding to fluxes at a lower scale. 


def fluctuations_distribute(field, flux, factor, ratio_bound):

# Define the extrapolated field and regional parameters.

    field_extrapolated = deepcopy(field)
    region_height = np.shape(field)[0]
    region_width = np.shape(field)[1]

    case = np.zeros((region_height, region_width))    # This variable stores the information about the pixels that were connected. The pixels that were connected are the ones where the `case' variable has the same integer value. The value grown with the `step' variable.
    case[field == 0] = - 1   # Wherever the pixels are irrelevant (land), the value of case is negative one.
    ratio = 0.0       # `ratio' measures the ratio of pixels at which the fluctuations have been redistributed.
    n_sea_pixels = len(field[field != 0]) + 0.0   # This is number of relevant pixels (sea).
    step = 1  # The value of `step' variable is initially set to 1.
          

# This is the main while loop which runs until the pixels where the fluctuations have been redistributed reach the desired ratio, when compared to all the relevant (sea) pixels. 

    while ratio < ratio_bound:

      # This randomly selects a grid point.
          
        pix_lon = round((region_width-1)*rn.random())
        pix_lat = round((region_height-1)*rn.random())

      # If the grid point is relevant (sea) please go on...
            
        if case[pix_lat, pix_lon] >= 0:
            
      # This randomly selects a neighbouring point: `first random' whether we move by 1 in longitudal direction, or latitudal direction. `Second random' determines whether we move by plus or minus one. 
          
            first_random = rn.random()
            second_random = rn.random()
        
          # To save us computational time, these conditions look on whether the randomly selected couples of points were already connected, and if they were it automatically flips the direction in the opposite.
      
            if (first_random > 0.5) & (pix_lon+1 < region_width) &  (pix_lon-1 >= 0): 
                if (case[pix_lat, pix_lon+1] == case[pix_lat, pix_lon]) & (case[pix_lat, pix_lon-1] == case[pix_lat, pix_lon]) & (case[pix_lat, pix_lon]*case[pix_lat, pix_lon-1]*case[pix_lat, pix_lon+1] != 0): 
                
                    first_random = 1.0 - first_random
                
            if (first_random < 0.5) & (pix_lat+1 < region_height) & (pix_lat-1 >= 0):
                if (case[pix_lat+1, pix_lon] == case[pix_lat, pix_lon]) & (case[pix_lat-1, pix_lon] == case[pix_lat, pix_lon]) & (case[pix_lat, pix_lon]*case[pix_lat-1, pix_lon]*case[pix_lat+1, pix_lon] != 0): 
                  
                    first_random = 1.0 - first_random

          # Again, to save computational time, the same thing is done with `second randomn' value (it automatically flips the motion to the opposite if the randomply selected option was realised).
               
            if (first_random > 0.5) & (second_random > 0.5) & (pix_lon+1 < region_width):
                if  (case[pix_lat, pix_lon+1] == case[pix_lat, pix_lon]) & (case[pix_lat, pix_lon]*case[pix_lat, pix_lon+1] != 0): 
                
                    second_random = 1.0 - second_random
                
            if (first_random > 0.5) & (second_random < 0.5) & (pix_lon-1 >= 0):
                if (case[pix_lat, pix_lon-1] == case[pix_lat, pix_lon]) & (case[pix_lat, pix_lon]*case[pix_lat, pix_lon-1] != 0): 
                
                    second_random = 1.0 - second_random
                
            if (first_random < 0.5) & (second_random > 0.5) & (pix_lat+1 < region_height):
                if (case[pix_lat+1, pix_lon] == case[pix_lat, pix_lon]) & (case[pix_lat+1, pix_lon]*case[pix_lat, pix_lon] != 0): 
  
                    second_random = 1.0 - second_random
                
            if (first_random < 0.5) & (second_random < 0.5) & (pix_lat-1 >= 0): 
                if (case[pix_lat-1, pix_lon] == case[pix_lat, pix_lon]) & (case[pix_lat, pix_lon]*case[pix_lat-1, pix_lon] != 0):
        
                    second_random = 1.0 - second_random
 


           # This is the part where the fluctuations are redistributed using the fluxes and the smoothened field information.

            # These outer conditions run through different combinations of randomly selected options: there are four - moving along longitude or latitude and then moving by plus or minus 1. Each case is dealt separately (the comments from the case 1 can be repeated for each other case).

# Case 1
                      
            if first_random > 0.5:

                if (second_random > 0.5) & (pix_lon+1 < region_width):  

                    # If value is at the sea and the two points were not connected before (they either were not connected with anything, or not mutually), please go on
                
                    if (field_extrapolated[pix_lat, pix_lon+1] != 0) & ((case[pix_lat, pix_lon] != case[pix_lat, pix_lon+1]) | (case[pix_lat, pix_lon]*case[pix_lat, pix_lon+1] == 0)):
                   
                        if field[pix_lat, pix_lon+1] > field[pix_lat, pix_lon]: 
                 
                            delta = flux[pix_lat, pix_lon]*factor

                        else: 

                            delta = - flux[pix_lat, pix_lon]*factor


                   # This distributes the fluctuation and records the connection in the ``case'' array.
                    

                        if case[pix_lat, pix_lon+1] == 0:
                        
                            field_extrapolated[pix_lat, pix_lon+1] = field_extrapolated[pix_lat, pix_lon] + delta
                        
                            if case[pix_lat, pix_lon] != 0:
                       
                                case[pix_lat, pix_lon+1] = case[pix_lat, pix_lon]

                            else:

                                case[pix_lat, pix_lon+1] = step

                                case[pix_lat, pix_lon] = step
                      
                        else:

                            diff = field_extrapolated[pix_lat, pix_lon] + delta - field_extrapolated[pix_lat, pix_lon+1]
                            field_extrapolated[case == case[pix_lat, pix_lon+1]] += diff
                            case[case == case[pix_lat, pix_lon+1]] = step
                            case[case == case[pix_lat, pix_lon]] = step 
                   
        # The same than in the previous step for different random selections...

# Case 2

                if (second_random < 0.5) & (pix_lon-1 >= 0):
                
                    if (field_extrapolated[pix_lat, pix_lon-1] != 0) & ((case[pix_lat, pix_lon] != case[pix_lat, pix_lon-1]) | (case[pix_lat, pix_lon]*case[pix_lat, pix_lon-1] == 0)):
                   
                        if field[pix_lat, pix_lon-1] > field[pix_lat, pix_lon]: 
                 
                            delta = flux[pix_lat, pix_lon]*factor

                        else: 

                            delta = - flux[pix_lat, pix_lon]*factor

                    
                        if case[pix_lat, pix_lon-1] == 0:
                        
                            field_extrapolated[pix_lat, pix_lon-1] = field_extrapolated[pix_lat, pix_lon] + delta
                        
                            if case[pix_lat, pix_lon] != 0:
                       
                                case[pix_lat, pix_lon-1] = case[pix_lat, pix_lon]

                            else:

                                case[pix_lat, pix_lon-1] = step

                                case[pix_lat, pix_lon] = step
                      
                        else:

                            diff = field_extrapolated[pix_lat, pix_lon] + delta - field_extrapolated[pix_lat, pix_lon-1]
                            field_extrapolated[case == case[pix_lat, pix_lon-1]] += diff
                            case[case == case[pix_lat, pix_lon-1]] = step
                            case[case == case[pix_lat, pix_lon]] = step 
                   
                    
# Case 3



            if first_random < 0.5:
              
                if (second_random > 0.5) & (pix_lat+1 < region_height): 
                
                    if (field_extrapolated[pix_lat+1, pix_lon] != 0) & ((case[pix_lat+1, pix_lon] != case[pix_lat, pix_lon]) | (case[pix_lat+1, pix_lon]*case[pix_lat, pix_lon] == 0)):
                   
                        if field[pix_lat+1, pix_lon] > field[pix_lat, pix_lon]: 
                 
                            delta = flux[pix_lat, pix_lon]*factor

                        else: 

                            delta = - flux[pix_lat, pix_lon]*factor

                    
                        if case[pix_lat+1, pix_lon] == 0:
                        
                            field_extrapolated[pix_lat+1, pix_lon] = field_extrapolated[pix_lat, pix_lon] + delta
                        
                            if case[pix_lat, pix_lon] != 0:
                       
                                case[pix_lat+1, pix_lon] = case[pix_lat, pix_lon]

                            else:

                                case[pix_lat+1, pix_lon] = step

                                case[pix_lat, pix_lon] = step
                      
                        else:

                            diff = field_extrapolated[pix_lat, pix_lon] + delta - field_extrapolated[pix_lat+1, pix_lon]
                            field_extrapolated[case == case[pix_lat+1, pix_lon]] += diff
                            case[case == case[pix_lat+1, pix_lon]] = step
                            case[case == case[pix_lat, pix_lon]] = step 
 
# Case 4
                  
                   
                if (second_random < 0.5) & (pix_lat-1 >= 0):

                    if (field_extrapolated[pix_lat-1, pix_lon] != 0) & ((case[pix_lat, pix_lon] != case[pix_lat-1, pix_lon]) | (case[pix_lat-1, pix_lon]*case[pix_lat, pix_lon] == 0)):
                   
                        if field[pix_lat-1, pix_lon] > field[pix_lat, pix_lon]: 
                 
                            delta = flux[pix_lat, pix_lon]*factor

                        else: 

                            delta = - flux[pix_lat, pix_lon]*factor

                    
                        if case[pix_lat-1, pix_lon] == 0:
                        
                            field_extrapolated[pix_lat-1, pix_lon] = field_extrapolated[pix_lat, pix_lon] + delta
                        
                            if case[pix_lat, pix_lon] != 0:
                       
                                case[pix_lat-1, pix_lon] = case[pix_lat, pix_lon]

                            else:

                                case[pix_lat-1, pix_lon] = step

                                case[pix_lat, pix_lon] = step
                      
                        else:

                            diff = field_extrapolated[pix_lat, pix_lon] + delta - field_extrapolated[pix_lat-1, pix_lon]
                            field_extrapolated[case == case[pix_lat-1, pix_lon]] += diff
                            case[case == case[pix_lat-1, pix_lon]] = step
                            case[case == case[pix_lat, pix_lon]] = step 
                      
   # After the step was done we update the ``ratio'' value and increase the ``step'' value.                    
            
            
            nonzero = len(case[np.where(case > 0)])
            ratio = nonzero/n_sea_pixels
            step += 1            
            
    return field_extrapolated








