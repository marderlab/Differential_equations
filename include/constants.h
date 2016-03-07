#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_



#include <cmath>
#include <limits>



//////////////////////// Special Floating Point Numbers ///////////////////////
#ifndef NAN_DEFINED
#define NAN_DEFINED
const double NaN __attribute__ ((unused)) =
  std::numeric_limits<double>::quiet_NaN();
#endif

#ifndef INF_DEFINED
#define INF_DEFINED
const double Inf __attribute__ ((unused)) =
  std::numeric_limits<double>::infinity();
#endif

#ifndef MACHINE_EPSILON_DEFINED
#define MACHINE_EPSILON_DEFINED
const double machineEpsilon __attribute__ ((unused)) =
  std::numeric_limits<double>::epsilon();
#endif



/////////////// Testing for Special Floating Point Numbers ////////////////////
inline bool
isNaN(const double & x)
{
  // note, this can evidently fail with some gcc compiler options
  return x != x;
}



//////////////////////////// Mathematical constants ///////////////////////////
#ifndef PI_DEFINED
#define PI_DEFINED
const double Pi __attribute__ ((unused)) = 4 * atan( (double) 1);
#endif



////////////////////////////// Physical constants /////////////////////////////
const double Q_Electron __attribute__ ((unused)) = 1.602176487e-19; // Coul
const double N_Avogadro __attribute__ ((unused)) = 6.02214179e23; //mol^-1
const double K_Boltzmann __attribute__ ((unused)) = 1.3806504e-23; //J/K
const double CelsiusToKelvin __attribute__((unused)) = 273.15;  //deg C



#endif
