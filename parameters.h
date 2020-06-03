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
const int end_time = 2880;
const int NUM_STEPS = 2000000;//100000;

//timestep
const double dt = (double)end_time/(double)NUM_STEPS;

//frequency of output for visualization
const int OUTPUT_FREQ = 6250;

//Cell parameters
const double average_G1_mother =  10.0;//minutes
const double average_G2_mother = 80.0;//minutes
const double average_radius =  3.1;//microns
const double extra_G1_daughter = 40.0;//minutes
//average length of G1 for daughter cell is 75 minutes
//45 minutes on mother, 30 minutes growing before
//daughter forms new bud of its own
const double size_at_div = 2.0;
const double DIV_SHIFT_RADIANS = .174533;

//adhesion between mother daughter cells
//const double k_adhesion_mother_daughter = 2*k_adhesion_cell_cell;
const double POISSON = .3;
const double ELASTIC_MOD = 1000;
const double K_ADH = 25;
const double K_BEND = 8;
const double ADHESION_STRENGTH = 12;
const double k_axial_frac = .5;
//chou model params
const double k_r = 1; 
const double k_a = .2;

//signaling params
//const double P_0 = 50;
//const double r_LOGISTIC = 1;
const double K_LOGISTIC = 100;
//const double A_LOGISTIC = 20;
//signaling params DNPM
const double alpha = .0154;
const double beta = 12;
const double Gamma = .00008;
const double mu = .0077;
const int n_0 = 1;
const int array_size = 4;
const double eta =  5.0;
const double THETA = .75;
//const bool ADHESION_ON = true;
const double SIGMA = .5;
const double OMEGA = 500;
//const double SINGLE_BOND_BIND_ENERGY = 25.0;
const double RECEPTOR_SURF_DENSITY = 1E15;
const double TEMPERATURE = 300;
const double KB = 1.38E-23;

//*************not using************************
//These two growth parameters set the timescale
const double k_g1 = .01;
const double k_g2 = .01;

//At initiation, cells are given a random maximum
//radii in the interval (.99*raius_average,1.01*radius_average)
//This is assigned in the cell constructor at the time
//of initiation based on the radius_average given below
const double radius_average = 4.375;

//A cell must reach nearly full size to start producing 
//a bud this is when it enters S
const double k_G1 = .9;

//When the bud reaches a threshold size
//the mother and daughter cell will divide (cleave)
const double k_mitosis = .5;

//the strength of the spring force handling
//the positional dynamics
//must be calibrated in relation to the growth rate 
//to allow the spring to push cells away from each other
const double k_spring = 1;

//calibrated base don length of the adhesion force and also
//used to define which cells are considered neighbors
//to produce adhesion it should be set higher than k^cell-cell_r_frac
const double k_neighbor = 1.1;

//spring force between mother and daughter bud
//updated during force calculation in code
//based on size of mother and bud radius

//since cells are not infinitely hard spheres, some amount
//of overlap is allowed
const double k_repulsion_cell_cell = 1;

//so that bud does not move away from mother
const double k_adhesion_mother_bud = 50;

//adhesion between unrelated cells
const double k_adhesion_cell_cell = .1;
//**************************************************************


#endif
