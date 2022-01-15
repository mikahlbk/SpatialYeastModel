//cell.cpp

//******************************************
//Include Dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include <chrono>
#include <cstdio>
#include <memory>
#include <random>
#include <omp.h>
#include <boost/array.hpp>
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "externs.h"
#include "colony.h"
#include "mesh_pt.h"
#include "mesh.h"
using namespace std;
//****************************************
//Public Member Functions for Cell.cpp

//Constructor for single founder cell
Cell::Cell(shared_ptr<Colony> my_colony, Coord cell_center, int rank, int Ti){
    //set this first to use random number generator
    this->my_colony    = my_colony;
    //cell location variables
    this->cell_center  = cell_center;
    //set bin id after it is made
    //cell growth variables
    double size_scaler = this->my_colony->roll_my_real_dis(2);
    size_scaler = (size_scaler -1)/10;
    this->max_radius        = (1+size_scaler)*AVERAGE_MAX_RADIUS;
    this->OG_max_radius     = max_radius;
    this->curr_radius       = max_radius;
    //cell cycle variables
    this->curr_phase    = 1;
    size_scaler = this->my_colony->roll_my_real_dis(2);
    size_scaler = (size_scaler -1)/10;
    this->my_G1_length    = (1+size_scaler)*AVERAGE_G1_MOTHER;
    this->my_OG_G1_length = my_G1_length;
    //cell lineage variables
    this->rank             = rank;
    return;
}//MBK
//Constructor for new bud
Cell::Cell(Coord cell_center,double init_radius,int Ti){
    //cell location
    this->cell_center = cell_center;
    //set bin id after made
    //cell growth
    this->curr_radius       = init_radius;
    double size_scaler      = this->my_colony->roll_my_real_dis(2);
    size_scaler             = (size_scaler -1)/10;
    this->max_radius        = (1+size_scaler)*AVERAGE_MAX_RADIUS;
    this->OG_max_radius     = max_radius;
    //cell cycle
    this->curr_phase           = 1;
    size_scaler                = this->my_colony->roll_my_real_dis(2);
    size_scaler                = (size_scaler -1)/10; 
    this->my_G1_length         = (1+size_scaler)*AVERAGE_G1_DAUGHTER;
    this->my_OG_G1_length      = my_G1_length;
   
    //cell lineage
    this->rank             = this->my_colony->get_num_cells();
    return;
}//MBK
//******Functions in order of appearance in colony.cpp***
shared_ptr<Mesh_Pt> Cell::find_my_nearest_bin(){
	shared_ptr<Cell> this_cell = shared_from_this();
	
	//cout << "Step 1" << endl;
	vector<shared_ptr<Mesh_Pt>> mesh_pts;
	//cout << "Step 2" << endl;
	this->get_my_colony()->get_my_mesh()->get_mesh_pts_vec(mesh_pts);
        //cout << "Step 3" << endl;
	double smallest_distance = 1.0e6;
	//cout << "Step 4" << endl;
	shared_ptr<Mesh_Pt> nearest_bin_id;
	//cout << "Step 5" << endl;
	double curr_distance;
	for(unsigned int i = 0; i < mesh_pts.size(); i++){
		curr_distance = (this->get_cell_center()-mesh_pts.at(i)->get_center()).length();
		if(curr_distance < smallest_distance){	
			smallest_distance = curr_distance;
			nearest_bin_id    = mesh_pts.at(i);
			//cout << "In find_bin function here is the current nutrient conc: " << mesh_pts.at(i)->get_my_nutrient_conc() << endl;
		}
	}
	//cout << "Step 6" << endl;
	this->bin_id = nearest_bin_id;		
	return nearest_bin_id;
}//MBK
void Cell::update_my_phase_length(){
	//determine the bin/mesh_pt pointer associated
	//with this cell
	shared_ptr<Mesh_Pt> my_bin = this->bin_id;
	//determine nutrient level in the cell's bin/mesh_pt
	double available_nutrients = my_bin->get_my_nutrient_conc();
	double my_new_G1_length    = 2.0*this->my_OG_G1_length*(1.0-available_nutrients) + this->my_OG_G1_length*available_nutrients;
	this->my_G1_length         = my_new_G1_length;
	return;
}//MBK
bool Cell::at_max_size(){
	bool my_bool;
	if(this->curr_radius >= this->max_radius){
		my_bool = true;
	}else{
		my_bool = false;
	}
	return my_bool;
}//MBK
void Cell::update_my_radius(){
    if(at_max_size()){
	    //do nothing
	}else{
	    double radius_increment = (this->max_radius)/my_G1_length;
		this->curr_radius       = curr_radius + radius_increment*dt;
	}
	return;
}//MBK
void Cell::perform_division(){
	if(this->at_max_size){
        double new_radius = this->get_curr_radius();
        //my new cell center
	    double my_new_center_x = this->get_cell_center().get_X()-.5*new_radius;
	    double my_new_center_y = this->get_cell_center().get_Y();
        this->set_curr_radius(new_radius);
        this->set_curr_phase(1);
	    Coord my_new_center = Coord(my_new_center_x,my_new_center_y);
        //daughter new cell center
	    double d_new_center_x = this->get_cell_center().get_X()+.5*new_radius;
	    double d_new_center_y = this->get_cell_center().get_Y();
	    Coord d_new_center = Coord(d_new_center_x,d_new_center_y);
        auto new_cell = make_shared<Cell>(d_new_center,new_radius);
    	new_cell->find_my_nearest_bin();
}//MBK
void Cell::compute_my_curr_force(){
	Coord force;
    	vector<shared_ptr<Cell>> neighbor_cells;
    	this->bin_id->get_neighbor_cells(neighbor_cells);
    	shared_ptr<Cell> this_cell = shared_from_this();
    	//cout << neighbor_cells.size() << endl;
	#pragma omp parallel
    	{
    		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:force) schedule(static,1)
		for(unsigned int i = 0; i < neighbor_cells.size(); i++){
    			if(neighbor_cells.at(i) != this_cell){
				force += this_cell->calc_forces_Hertz(neighbor_cells.at(i));
				//#pragma omp critical
				//cout << "Cell Rank: " << rank << "And force: " << force << endl;
			}
		}
     	}
	//cout << curr_force << endl;
     	this->curr_force = force;
     	return;
}//MBK
Coord Cell::calc_forces_linear_spring(shared_ptr<Cell> my_neighbor){
	//Describe force calculation in paper
    	Coord my_loc = cell_center;
    	double my_radius = curr_radius;
    	Coord neighbor_loc;
    	double neighbor_radius;
    	Coord force = Coord(0,0);
        //distance between two cell centers
    	double d_ij;
    	//vector for direction of force
    	Coord v_ij;
        shared_ptr<Cell> this_cell = shared_from_this();
    	if(my_neighbor != this_cell){
		    neighbor_loc = my_neighbor->get_cell_center();
		    neighbor_radius = my_neighbor->get_curr_radius();
		    d_ij = (my_loc-neighbor_loc).length();
		    v_ij = (my_loc - neighbor_loc);
            if((d_ij - (my_raduis+neighbor_radius)) < .001){
                force += v_ij*-1*K_ADH*(.001);
            }else{
	            force += v_ij*-1*K_ADH*(d_ij-(my_radius+neighbor_radius));	
			}
		}	
    curr_force = force;
    return curr_force;
}//MBK
void Cell::update_my_location(){
	this->cell_center = this->cell_center + this->curr_force*(1.0/(ETA*(1.0 + this->curr_radius)))*dt;
    	return;
}//MBK
void Cell::print_location_file_format(ofstream& ofs){
	ofs << rank << " " << cell_center.get_X() << " " << cell_center.get_Y() << " " << curr_radius << endl;
	return;
}//MBK

//***END OF FILE***
//**************************************************************************************************************************
