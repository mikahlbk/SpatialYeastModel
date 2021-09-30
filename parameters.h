// parameters.h

//************************************
// Include Guard
#ifndef _PARAMETERS_H_INCLUDED_
#define _PARAMETERS_H_INCLUDED_

//***********************************
/// Forward Declarations

//***********************************
// Include Dependencies
#include <math.h>

//***********************************
//Simulation Constants

//time in minutes
const int end_time = 1440;
const int NUM_STEPS = 1000000;

//timestep
//minutes
const double dt = (double)end_time/(double)NUM_STEPS;

//frequency of output for visualization
//6250*dt = how many minutes
const int OUTPUT_FREQ = 6250;
const int BIN_UPDATE_INCREMENT = 6250;
//Cell parameters
const double AVERAGE_MAX_RADIUS = 2.58;//microns
const double AVERAGE_G1_MOTHER  = 15;//minutes
const double DAUGHTER_INIT_RADIUS = 0.0;//microns
const double AVERAGE_SEPARATION_RADIUS = 2.35;//microns
const double AVERAGE_G1_DAUGHTER = 45;//minutes
const double AVERAGE_BUD_PHASE = 75;//minutes
const double DIV_SHIFT_RADIANS = .174533;

//Hertz potential parameters
const double SINGLE_BOND_BIND_ENERGY = 25.0;
const double RECEPTOR_SURF_DENSITY = 1E15;
const double TEMPERATURE = 300;
const double KB = 1.38E-23;
const double POISSON = .3;
const double ELASTIC_MOD = 1000;
const double K_ADH = 25;
const double ETA = 2.5;

//Protein conceentration parameters
const double K_LOGISTIC = 100;
const double A_LOGISTIC = 36;
//**************************************************************


#endif
