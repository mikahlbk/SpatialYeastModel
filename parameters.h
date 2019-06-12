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

//timestep
const double dt = .1;

//These two growth parameters set the timescale
const double k_g1 = .01;
const double k_g2 = .01;

//At initiation, cells are given a random maximum
//radii in the interval (.99*raius_average,1.01*radius_average)
//This is assigned in the cell constructor at the time
//of initiation based on the radius_average given below
const double radius_average = .6;

//A cell must reach nearly full size to start producing 
//a bud this is when it enters S
const double k_G1 = .9;

//When the bud reaches a threshold size
//the mother and daughter cell will divide (cleave)
const double k_mitosis = .7;

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
const double k_repulsion_cell_cell = .9;

//so that bud does not move away from mother
const double k_adhesion_mother_bud = 1;

//adhesion between unrelated cells
const double k_adhesion_cell_cell = .5;

//adhesion between mother daughter cells
const double k_adhesion_mother_daughter = 2*k_adhesion_cell_cell;

//chou model params
const double k_r = .9; 
const double k_a = 1.1;
#endif
