#! /usr/bin/env python
#
#
def shallow_water_1d_test ( ):   #define the function shallow_water_1d_test
#*****************************************************************************80
#
## SHALLOW_WATER_1D Python Code.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 Oct. 2019
#
#  Author:
#
#    Xin Liu
#
  print ( '' )
  print ( 'Software:  HMX-1D' )
  print ( 'Author:    Xin Liu' )
#  print ( '' )
#  print ( 'Computational Resources Info.:' )
#  print ( '  Python version: %s' % ( platform.python_version ( ) ) )   
#  print ( '  Machine type: %s' % ( platform.machine() ) )
#  print ( '  Processor name: %s' % ( platform.processor() ) )
#  print ( '  OS name: %s %s' % ( platform.system(), platform.release() ))


#call modules and rename them
#  import matplotlib.pyplot as plt  # module for plottings
  import numpy as np               # numpy: numerical python which is a library for scientfic comput.It's in general used to deal with arrays
#  Set parameters.
  nx = 41         # not cell number, but nodes number
  nt = 100        # define number of time step, will modify it later based on CFL
  x_length = 1.0  #domain length
  t_length = 0.2  #total simulation time
  g = 9.81        #g, you know it, if not, delete this program and go to watch 'Peppa Pig'
#
#  Compute h and uh.
#
  h_new, uh_new, x, t = shallow_water_1d ( nx, nt, x_length, t_length, g ) 
#----- output results in arrays  -------------------#
#
##****************  Write data.  ***********************##
  filename_x = 'x.txt'
  filename_h = 'h.txt'
  filename_uh = 'uh.txt'
  r8vec_write ( filename_x, nx, x )      # r8: double precision real, vec:vector
  r8vec_write ( filename_h, nx, h_new )   
  r8vec_write ( filename_uh, nx, uh_new )
#  Terminate.
  print ( '' )
#  print ( 'HMX-1D' )
  print ( 'Execution ends normally.' )
  return




#*****************************************************************************!!
def shallow_water_1d ( nx, nt, x_length, t_length, g ):
#
#      dh/dt + d uh/dx = 0
#      duh/dt + d ( u^2 h + 1/2 g h^2 )/dx = 0
#
  import numpy as np
  import platform    # used to access the underlying platform’s data, such as, hardware, operating system, and interpreter version information

#----- define variales in column array -------#
  h = np.zeros ( nx )    #zero array first
  uh = np.zeros ( nx )
  hm = np.zeros ( nx - 1 )   # hm are cell center values of h
  uhm = np.zeros ( nx - 1 )
  x = np.zeros ( nx )
  t = np.zeros ( nt + 1 )
  h_new = np.zeros ( nx ) #np.zeros ( [ nx, nt + 1 ] )   # matrix with nx rows and nt+1 columns, the first column for initial solution at t=0 
  uh_new = np.zeros ( nx ) #np.zeros ( [ nx, nt + 1 ] )
#
#  Define the locations of the nodes and time steps and the spacing.
#
  x = np.linspace ( 0, x_length, nx )    # np.linspace returns evenly spaced numbers over a specified interval
  t = np.linspace ( 0, t_length, nt + 1 )

 # dx = x_length / float ( nx - 1 )  # float() makes sure dividing by a number, not string
 # dt = t_length / float ( nt )

  dx = x_length /  ( nx - 1 )  # float() makes sure dividing by a number, not string
  dt = t_length /  nt 
  h, uh = initial_conditions ( nx, nt, h, uh, x )        # Apply the initial conditions.
  h, uh = boundary_conditions ( nx, nt, h, uh, t[0] )    # Apply the boundary conditions; t[0] means initial solutions; bc applies to the boundary nodes
  h_new[0:nx] = h[0:nx]       #  Store the initial condition into the arrays.
  uh_new[0:nx] = uh[0:nx]
  
  for it in range ( 1, nt + 1 ):
#
#  Take a half time step, estimating H and UH at the NX-1 spatial midpoints.
#
    hm[0:nx-1] = ( h[0:nx-1] + h[1:nx] ) / 2.0 \
      - ( dt / 2.0 ) * ( uh[1:nx] - uh[0:nx-1] ) / dx

    uhm[0:nx-1] = ( uh[0:nx-1] + uh[1:nx] ) / 2.0 \
      - ( dt / 2.0 ) * ( \
        uh[1:nx] ** 2    / h[1:nx]   + 0.5 * g * h[1:nx] ** 2 \
      - uh[0:nx-1] ** 2  / h[0:nx-1] - 0.5 * g * h[0:nx-1] ** 2 ) / dx
#
#  Take a full time step, evaluating the derivative at the half time step,
#  to estimate the solution at the NX-2 nodes.
#
    h[1:nx-1] = h[1:nx-1] \
      - dt * ( uhm[1:nx-1] - uhm[0:nx-2] ) / dx

    uh[1:nx-1] = uh[1:nx-1] \
      - dt * ( \
        uhm[1:nx-1] ** 2  / hm[1:nx-1] + 0.5 * g * hm[1:nx-1] ** 2 \
      - uhm[0:nx-2] ** 2  / hm[0:nx-2] - 0.5 * g * hm[0:nx-2] ** 2 ) / dx
#
#  Update the boundary conditions.
#
    h, uh = boundary_conditions ( nx, nt, h, uh, t[it] )    # not sure if this is neccessary to be here
#
#  Copy data into the big arrays.
#
    h_new[0:nx] = h[0:nx]   #copy the results of it step into the it column of data array
    uh_new[0:nx] = uh[0:nx]

#  Terminate.
#
#  print ( '' )
#  print ( 'Core Solver:' )
#  print ( '  Normal end of execution.' )

#  return h_array, uh_array, x, t 
  return h_new, uh_new, x, t



##***************************************************************!!
def boundary_conditions ( nx, nt, h, uh, t ):
  bc = 1   # set the type of boundary condition 

  if ( bc == 1 ):      #  Periodic boundary conditions on H and UH.
    h[0] = h[nx-2]     
    h[nx-1] = h[1]
    uh[0] = uh[nx-2]
    uh[nx-1] = uh[1]  

  elif ( bc == 2 ):     #  Free boundary conditions on H and UH.
    h[0] = h[1]
    h[nx-1] = h[nx-2]
    uh[0] = uh[1]
    uh[nx-1] = uh[nx-2]

  elif ( bc == 3 ):      #  Reflective boundary conditions on UH, free boundary conditions on H.
    h[0] = h[1]
    h[nx-1] = h[nx-2]
    uh[0] = - uh[1]
    uh[nx-1] = - uh[nx-2]

  return h, uh

##**************************************************************!!
def initial_conditions ( nx, nt, h, uh, x ):

  import numpy as np

  h = 2.0 + np.sin ( 2.0 * np.pi * x )
  uh = np.zeros ( nx )

  return h, uh


##************ Writes a double presioucs column array to a file **************##
def r8vec_write ( filename, n, a ):
  output = open ( filename, 'w' )
  for i in range ( 0, n ):
    s = '  %g\n' % ( a[i] )
    output.write ( s )
  output.close ( )
  return

##******************** Get timestamp *****************************## 
def timestamp ( ):   #** TIMESTAMP prints the date as a timestamp.
  import time
  t = time.time ( )
  print ( time.ctime ( t ) )
  return None

##******************** Main Program *****************************##
#
if ( __name__ == '__main__' ):
#  print ( 'Start time:'  ) 
  timestamp ( )
  shallow_water_1d_test ( )
  print ( '' )
#  print ( 'End time:' )
  timestamp ( )

