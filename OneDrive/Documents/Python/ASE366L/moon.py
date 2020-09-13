#!/usr/local/bin/python
'''

 Moon-specific functions


 Functions used to get properties of the Moon

'''

__author__ = 'Brandon A. Jones'
__version__ = '$Revision: 652 $'[11:-2]
__date__ = '$Date: 2016-05-03 12:12:02 -0500 (Tue, 03 May 2016) $'[7:26]

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################

import numpy as np

#  Pre-compute the conversion from degrees to radians
D2R_ = np.pi/180.0

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

#-------------------------------------------------------------------------------
def getMoonPos_wrt_Earth( JD_TDB ) :
  '''
  Get the position of the Moon with respect to the Earth based on the 
  2010 Astronomical Almanac algorithms.  Returns position vector in
  J2000, geocentric, equatorial frame.
  '''

  #  Convenience functions to make copy/paste from MATLAB easier
  def sind( arg ) :
    return np.sin( arg * D2R_ )

  def cosd( arg ) :
    return np.cos( arg * D2R_ )

  # Convert the input time into Julian Centuries
  T_TDB = ( JD_TDB - 2451545.0  ) / 36525.0

  long_ecliptic= 218.32  + 481267.8813 *T_TDB \
                         + 6.29 *sind( 134.9 +477198.85 *T_TDB ) \
                         - 1.27 *sind( 259.2 -413335.38 *T_TDB ) \
                         + 0.66 *sind( 235.7 +890534.23 *T_TDB ) \
                         + 0.21 *sind( 269.9 +954397.70 *T_TDB ) \
                         - 0.19 *sind( 357.5 + 35999.05 *T_TDB ) \
                         - 0.11 *sind( 186.6 +966404.05 *T_TDB )      # deg

  lat_ecliptic =   5.13 *sind( 93.3 +483202.03 *T_TDB ) \
                 + 0.28 *sind( 228.2 +960400.87 *T_TDB ) \
                 - 0.28 *sind( 318.3 +  6003.18 *T_TDB ) \
                 - 0.17 *sind( 217.6 -407332.20 *T_TDB )      # deg

  parralax =  0.9508  + 0.0518 *cosd( 134.9 +477198.85 *T_TDB ) \
                      + 0.0095 *cosd( 259.2 -413335.38 *T_TDB ) \
                      + 0.0078 *cosd( 235.7 +890534.23 *T_TDB ) \
                      + 0.0028 *cosd( 269.9 +954397.70 *T_TDB )    # deg


  #  Use the unwrap() function to get an angle in the range [-pi,pi]
  long_ecliptic = np.unwrap( [long_ecliptic*D2R_] )[0]
  lat_ecliptic  = np.unwrap( [lat_ecliptic*D2R_] )[0]
  parralax      = np.unwrap( [parralax*D2R_] )[0]

  obliquity = 23.439291  - 0.0130042 *T_TDB  #deg
  obliquity = obliquity * D2R_

  # Pre-compute cosine and sine values that will be needed at least twice
  cos_lag_ecl = np.cos( lat_ecliptic )
  sin_lon_ecl = np.sin( long_ecliptic )
  sin_lat_ecl = np.sin( lat_ecliptic )
  cos_obliq   = np.cos( obliquity )
  sin_obliq   = np.sin( obliquity )

  # ------------- calculate moon position vector ----------------
  pos_j2000 = np.zeros((3,))
  pos_j2000[0] = cos_lag_ecl * np.cos(long_ecliptic)
  pos_j2000[1] = cos_obliq*cos_lag_ecl*sin_lon_ecl - sin_obliq*sin_lat_ecl
  pos_j2000[2] = sin_obliq*cos_lag_ecl*sin_lon_ecl + cos_obliq*sin_lat_ecl
          
  #  Generate position in km where Re = 6378.1363 km
  pos_j2000 *= 6378.1363/np.sin(parralax)

  #  Return position 
  return pos_j2000


################################################################################
#             U N I T     T E S T     C A S E     F U N C T I O N:
################################################################################

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():

  np.set_printoptions(precision=15)
  reftime = 2458927.07792
  print(getMoonPos_wrt_Earth( reftime ) )

  return

if __name__ == "__main__":
  main()
