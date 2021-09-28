// colony.cpp

//****************************************
// Inlcude Dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <memory>
#include <random>
#include "externs.h"
#include "parameters.h"
#include "colony.h"
#include "coord.h"
#include "cell.h"
#include "mesh.h"
#include "mesh_pt.h"
using namespace std;
//*****************************************
//Public Member Functions for Colony.cpp

//constructor
Colony::Colony(shared_ptr<Mesh> new_mesh) {
	this->my_mesh = new_mesh;
	int new_seed = time(0);
	this->my_engine.seed(new_seed);		
	return;
}
double Colony::roll_my_real_dis(int upper_bound){
	double random_num;
	if(upper_bound == 1){
		std::uniform_real_distribution<double> my_zero_one_dis(0.0,1.0);
		random_num = my_zero_one_dis(this->my_engine);
	}else{
		std::uniform_real_distribution<double> my_zero_two_dis(0.0,2.0);
		random_num = my_zero_two_dis(this->my_engine);

	}
	return random_num;
}
void Colony::make_founder_cell(shared_ptr<Colony> this_colony,int Ti){
     //shared_ptr<Colony> this_colony = shared_from_this();
     //parameters to be fed into founder cell constructor
     Coord center = Coord(0,0);
     int rank = 0;
     double div_site = 2*M_PI*(this_colony->roll_my_real_dis(1));
     //make the founder cell
     auto new_cell = make_shared<Cell>(this_colony, center, rank, div_site,Ti);
     new_cell->find_my_nearest_bin();
     this->update_colony_cell_vec(new_cell);
     return;
}
void Colony::update_colony_cell_vec(shared_ptr<Cell> new_cell){
	my_cells.push_back(new_cell);
	return;
}
void Colony::update_cell_bin_ids(){
	this->my_mesh->clear_mesh_pt_cells();
	shared_ptr<Mesh_Pt> my_bin;
//	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_bin = my_cells.at(i)->find_my_nearest_bin();
//		#pragma omp critical
		my_bin->update_my_cells(my_cells.at(i));
	}	
	return;
}
void Colony::update_cell_phase_lengths(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i< my_cells.size();i++){
		my_cells.at(i)->update_my_phase_length();
	}
	return;
}
void Colony::update_cell_radii(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_cells.at(i)->update_my_radius();
	}
	return;
}
void Colony::perform_bud_separation(){
    	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i < my_cells.size(); i++){
		my_cells.at(i)->perform_bud_separation();
	}
	return;
}
void Colony::update_cell_cycle_phases(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_cells.at(i)->update_my_phase();
	}
	return;
}
void Colony::add_buds(int Ti){
	//cout << "My guess is here" << endl;
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		shared_ptr<Cell> new_bud = my_cells.at(i)->add_bud(Ti);
		if(new_bud){
			#pragma omp critical
			my_cells.push_back(new_bud);
		}
	}
	//cout << "Am I wrong?" << endl;
	return;
}
void Colony::update_cell_protein_concentrations(){
   #pragma omp parallel for schedule(static,1) 
   for(unsigned int i=0; i< my_cells.size(); i++){
        my_cells.at(i)->update_my_protein_concentration();
    }
	return;
}
void Colony::update_cell_locations(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_cells.at(i)->compute_my_curr_force();
	}
	#pragma omp parallel for schedule(static,1)
        for(unsigned int i = 0; i < my_cells.size(); i++){
        	my_cells.at(i)->update_my_location();
	}
	return;
}
void Colony::update_mesh_nutrient_concentration(){
	vector<shared_ptr<Mesh_Pt>> mesh_pts;
	this->my_mesh->get_mesh_pts_vec(mesh_pts);	
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < mesh_pts.size(); i++){
		//cout << "Mesh: " << i << endl;
		mesh_pts.at(i)->update_my_nutrient_concentration();
	}
	return;
}
void Colony::write_data(ofstream& ofs1, ofstream& ofs2){
    for(unsigned int i = 0; i < my_cells.size();i++){
        my_cells.at(i)->print_location_file_format(ofs1);
	my_cells.at(i)->print_cell_cycle_file_format(ofs2);
    }
    return;
}

//*******************************************************
//End of File
