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
Cell::Cell(shared_ptr<Colony> my_colony, Coord cell_center, int rank, double new_div_site,int Ti){
    //set this first to use random number generator
    this->my_colony    = my_colony;
    //cell location variables
    this->cell_center  = cell_center;
    //set bin id after it is made
    //cell growth variables
    double size_scaler = this->my_colony->roll_my_real_dis(2);
    size_scaler = (size_scaler -1)/10;
    this->max_radius        = (1+size_scaler)*AVERAGE_MAX_RADIUS;
    this->OG_max_radius = max_radius;
    this->curr_radius       = max_radius;
    this->separation_radius = 0;//not used for founder
    this->OG_separation_radius = separation_radius;
    //cell cycle variables
    this->curr_phase    = 1;
    this->my_G1_tracker = 0;
    this->my_bud_phase_length = 0;//not used for founder
    this->my_OG_bud_phase_length = my_bud_phase_length;
    size_scaler = this->my_colony->roll_my_real_dis(2);
    size_scaler = (size_scaler -1)/10;
    this->my_G1_length  = (1+size_scaler)*AVERAGE_G1_MOTHER;
    this->my_OG_G1_length = my_G1_length;
    this->is_mother     = true;
    this->is_daughter   = false;
    this->is_bud        = false;
    this->has_bud       = false;
    //cell lineage variables
    this->rank             = rank;
    this->successor_number = 0;
    this->division_age     = 0;
    int birth_time         = Ti;
    this->time_born        = birth_time;
    shared_ptr<Cell> my_mom(nullptr);
    this->mother           = my_mom;
    this->lineage_vec.push_back(rank);
    this->curr_div_site    = new_div_site;
    this->div_sites_vec.push_back(new_div_site);
    //cell protein variables
    this->lambda       = LAMBDA;
    this->OG_lambda    = lambda;
    this->rho          = RHO;
    this->curr_protein = P_0;
    return;
}//MBK
//Constructor for new bud
Cell::Cell(Coord cell_center,shared_ptr<Cell> mother,double div_site,int Ti){
    //need first
    this->mother    = mother;
    this->my_colony = mother->get_my_colony();
    double size_scaler;
    //cell location
    this->cell_center = cell_center;
    //set bin id after made
    //cell growth
    this->curr_radius       = DAUGHTER_INIT_RADIUS;
    size_scaler             = this->my_colony->roll_my_real_dis(2);
    size_scaler             = (size_scaler -1)/10;
    this->separation_radius = (1+size_scaler)*AVERAGE_SEPARATION_RADIUS;
    this->OG_separation_radius = separation_radius;
    size_scaler             = this->my_colony->roll_my_real_dis(2);
    size_scaler             = (size_scaler -1)/10;
    this->max_radius        = (1+size_scaler)*AVERAGE_MAX_RADIUS;
    this->OG_max_radius     = max_radius;
    //cell cycle
    this->curr_phase           = 1;
    size_scaler                = this->my_colony->roll_my_real_dis(2);
    size_scaler                = (size_scaler -1)/10; 
    this->my_G1_length         = (1+size_scaler)*AVERAGE_G1_DAUGHTER;
    this->my_OG_G1_length      = my_G1_length;
    this->my_G1_tracker        = 0.0;
    size_scaler                = this->my_colony->roll_my_real_dis(2);
    size_scaler                = (size_scaler -1)/10; 
    this->my_bud_phase_length  = (1+size_scaler)*AVERAGE_BUD_PHASE;
    this->my_OG_bud_phase_length = my_bud_phase_length;
    this->is_mother   = false;
    this->is_daughter = false;
    this->is_bud      = true;
    this->has_bud     = false;
    //cell lineage
    this->rank             = this->my_colony->get_num_cells();
    this->successor_number = this->mother->get_division_age();
    this->division_age     = 0;
    int birth_time         = Ti;
    this->time_born        = birth_time;
    vector<int> new_lineage_vec;
    this->mother->get_lineage_vec(new_lineage_vec);
    this->lineage_vec      = new_lineage_vec;
    this->lineage_vec.push_back(mother->get_rank());
    this->curr_div_site   = div_site;
    this->lambda          = LAMBDA;
    this->OG_lambda       = lambda;
    this->rho             = RHO; 
    this->curr_protein    = this->mother->get_curr_protein();
    return;
}//MBK
//***getters that need a whole function***
void Cell::get_lineage_vec(vector<int>& curr_lineage_vec){
	curr_lineage_vec = this->lineage_vec;
	return;
}
void Cell::get_daughter_vec(vector<int>& curr_daughter_vec){
	curr_daughter_vec = this->daughter_vec;
	return;
}
void Cell::get_div_sites_vec(vector<double>& curr_div_sites_vec){
	curr_div_sites_vec = this->div_sites_vec;
	return;
}//ALL 3 GOOD-MBK
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
	if(this->check_is_mother()){//this is a fully grown cell that has had at least one division
		//the only thing the nutrients will affect is G1 length 
		//AND protein replication rate
		//calculate new G1 length
		//cout << "My bins available nutrients: " << available_nutrients << endl;
		double my_new_G1_length = 2.0*this->my_OG_G1_length*(1.0-available_nutrients) + this->my_OG_G1_length*available_nutrients;
		this->my_G1_length      = my_new_G1_length;
		//cout << "This is cell " << rank << " and has G1 " << my_new_G1_length << endl;
		//cout << "This is cell " << rank << "and has G1 " << my_G1_length << endl;
		double my_new_lambda    = .5*this->OG_lambda*(1-available_nutrients) + this->OG_lambda*available_nutrients;
		this->lambda            = my_new_lambda; 
	}else{//this is either a new daughter still in G1 or a bud
		if(this->check_is_bud()){
			//if its a bud we need to adjust its budding phase length and size at separation
			//calculate new budding phase length
			double my_new_budding_length = 2*this->my_OG_bud_phase_length*(1-available_nutrients) + this->my_OG_bud_phase_length*available_nutrients;
		        this->my_bud_phase_length    = my_new_budding_length;
			//calculate new size at separation
			double radius_scaler    = .6*(1-available_nutrients) + available_nutrients;
			this->separation_radius = radius_scaler*OG_separation_radius;
		}else if(this->check_is_daughter()){
			//if its not a bud then its a new daughter that is still in G1
			//we need to adjust its G1 length, max_radius and lambda
			//caculate new G1 length
			double my_new_G1_length = 2*this->my_OG_G1_length*(1-available_nutrients) + this->my_OG_G1_length*available_nutrients;
			this->my_G1_length      = my_new_G1_length;
			//calculate max_radius
			double radius_scaler    = .7*(1-available_nutrients) + available_nutrients;
			this->max_radius        = radius_scaler*OG_max_radius;
			double my_new_lambda    = .5*this->OG_lambda*(1-available_nutrients) + this->OG_lambda*available_nutrients;
			this->lambda            = my_new_lambda;
		}else{
			exit(100);
		}
	}
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
	if(this->check_is_mother()){//this is a fully grown mother cell
		//do nothing
	}else if(this->check_is_daughter()){
		if(at_max_size()){
			//daughter is done growing and bud will be added
		}else{
			//daughter not done growing need to increase radius
			double radius_increment = (this->max_radius-this->separation_radius)/my_G1_length;
			this->curr_radius       = curr_radius + radius_increment*dt;
		}
	}else if(this->check_is_bud()){
		if(this->curr_radius < this->separation_radius){
			double radius_increment = (this->separation_radius - DAUGHTER_INIT_RADIUS)/my_bud_phase_length;
			this->curr_radius       = curr_radius + radius_increment*dt;
		} 
	}else{
		exit(100);
	}
	return;
}//MBK
void Cell::perform_bud_separation(){
	if(this->check_is_bud()){
		if(this->curr_radius >= separation_radius){
			shared_ptr<Cell> my_mom = this->get_mother();	
			//At this point the bud has grown big 
			//enough it will separate from it's mother.
			//The bud becomes a daughter cell, change info. 
			this->curr_phase    = 1;
			//G1 tracker isnt used because G1 for new daughters
			//is tracked by size
			this->is_bud       = false;
			this->is_daughter  = true;
			this->curr_protein = my_mom->get_curr_protein()*this->rho; 
			//The mother loses its bud, change info.
			//when a new daughter cell separates from it's
			//first bud it becomes a mother cell
			if(my_mom->check_is_daughter()){
				my_mom->set_is_daughter(false);
				my_mom->set_is_mother(true);
			}
			my_mom->set_has_bud(false);
			my_mom->set_curr_bud(nullptr);
			my_mom->set_curr_phase(1); 
			my_mom->set_my_G1_tracker(0.0);
			double mother_rho = my_mom->get_rho();
			my_mom->set_curr_protein(my_mom->get_curr_protein()*(1-mother_rho));
		}
	}
	return;
}//MBK
void Cell::update_my_phase(){
	if((this->check_is_mother())){
		if(this->curr_phase ==1){//update G1 tracker accordingly
			double G1_tracker_increment = 1.0/(this->my_G1_length);
			this->my_G1_tracker = my_G1_tracker + G1_tracker_increment*dt;
			if(this->my_G1_tracker >= 1.0){
				this->curr_phase = 2;
			}
		}
	}else if(this->check_is_daughter()){
		if((this->at_max_size())&&(this->curr_phase==1)){
			this->curr_phase = 2;
		}
	}else if(this->check_is_bud()){
		//do nothing
	}else{
		exit(100);
	}
	return;
}//MBK
shared_ptr<Cell>  Cell::add_bud(int Ti){
	shared_ptr<Cell> mother_cell = shared_from_this();
	shared_ptr<Cell> my_bud;
	if(this->check_is_mother()){
		if(this->curr_phase == 2){
			this->curr_div_site = this->determine_new_div_site();
			this->div_sites_vec.push_back(curr_div_site);
			my_bud = this->perform_budding(mother_cell,curr_div_site,Ti);
			//update colony cell vec
    			//this->my_colony->update_colony_cell_vec(my_bud);
			this->has_bud = true;
			this->division_age = division_age+1;
			this->curr_bud = my_bud;
			this->daughter_vec.push_back(my_bud->get_rank());
			int div_time   = Ti;
			this->division_times_vec.push_back(div_time);
			this->curr_phase = 3;
		}
	}else if(this->check_is_daughter()){
		if((this->at_max_size()) && (this->curr_phase == 2)){
			this->curr_div_site = this->determine_new_div_site();
			this->div_sites_vec.push_back(curr_div_site);
			my_bud = this->perform_budding(mother_cell,curr_div_site,Ti);
			//update colony cell vec
    			//this->my_colony->update_colony_cell_vec(my_bud);
			this->has_bud = true;
			this->division_age = division_age+1;
			this->curr_bud = my_bud;
			this->daughter_vec.push_back(my_bud->get_rank());
			int div_time   = Ti;
			this->division_times_vec.push_back(div_time);
			this->curr_phase = 3;
		}
	}else if(this->check_is_bud()){
		//do nothing
	}else{
		exit(100);
	}
	return my_bud;
}//MBK
double Cell::determine_new_div_site(){
	bool used = true;
	double curr_site;
	if(this->my_colony->roll_my_real_dis(1) < .5){//50% axial
		curr_site = this->curr_div_site + DIV_SHIFT_RADIANS;
		if(curr_site > 2*M_PI){
			curr_site = curr_site - 2*M_PI;
		}
		if(std::find(this->div_sites_vec.begin(), this->div_sites_vec.end(), curr_site) != this->div_sites_vec.end()){//need to move
			//move 10 degrees clockwise until new spot found
    			do{
				curr_site = curr_site + DIV_SHIFT_RADIANS;
				if(std::find(this->div_sites_vec.begin(), this->div_sites_vec.end(), curr_site) == this->div_sites_vec.end()){
					used = false;
				}	
     			}while(used);
		}

	}else{//50% distal
		curr_site = this->curr_div_site + M_PI;
		if(curr_site > 2*M_PI){
			curr_site = curr_site - 2*M_PI;
		}
		if(std::find(this->div_sites_vec.begin(), this->div_sites_vec.end(), curr_site) != this->div_sites_vec.end()){//need to move
			//move 10 degrees clockwise until new spot found
    			do{
				curr_site = curr_site + DIV_SHIFT_RADIANS;
				if(std::find(this->div_sites_vec.begin(), this->div_sites_vec.end(), curr_site) == this->div_sites_vec.end()){
					used = false;
				}	
     			}while(used);

		}
	}
	return curr_site;
}//MBK
shared_ptr<Cell> Cell::perform_budding(shared_ptr<Cell> mother_cell,double new_div_site,int Ti){
	//cell center
	double new_center_x = mother_cell->get_cell_center().get_X()+(mother_cell->get_curr_radius() + DAUGHTER_INIT_RADIUS)*cos(mother_cell->get_curr_div_site());
	double new_center_y = mother_cell->get_cell_center().get_Y()+(mother_cell->get_curr_radius() + DAUGHTER_INIT_RADIUS)*sin(mother_cell->get_curr_div_site());
	Coord new_center = Coord(new_center_x,new_center_y);
	auto new_cell = make_shared<Cell>(new_center,mother_cell,new_div_site,Ti);
    	new_cell->find_my_nearest_bin();
	return new_cell;
}//MBK
void Cell::update_my_protein_concentration(){
	double new_protein_conc;
	if(BISTABLE == 1){
		double scalar = this->my_colony->roll_my_real_dis(2);
    		scalar = (scalar -1)/10;
		double A_param = scalar*A_LOGISTIC;
		new_protein_conc   = this->curr_protein + this->lambda*curr_protein*(1-curr_protein/K_LOGISTIC)*(curr_protein/A_param -1)*dt;
    		this->curr_protein = new_protein_conc;
	}else{ 
		new_protein_conc   = this->curr_protein + this->lambda*curr_protein*(1-curr_protein/K_LOGISTIC)*dt;
		this->curr_protein = new_protein_conc;
	}
	if(this->check_has_bud()){
		this->curr_bud->set_curr_protein(new_protein_conc);
	}	
	return;
}//MBK
void Cell::compute_my_curr_force(){
	Coord force;
    	vector<shared_ptr<Cell>> neighbor_cells;
    	this->bin_id->get_neighbor_cells(neighbor_cells);
    	shared_ptr<Cell> this_cell = shared_from_this();
    	#pragma omp parallel
    	{
    		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:force) schedule(static,1)
		for(unsigned int i = 0; i < neighbor_cells.size(); i++){
    			if(neighbor_cells.at(i) != this_cell){
				force += this_cell->calc_forces_Hertz(neighbor_cells.at(i));
			}
		}
     	}
     	this->curr_force = force;
     	return;
}//MBK
Coord Cell::calc_forces_Hertz(shared_ptr<Cell> my_neighbor){
	//Describe force calculation in paper
    	Coord my_loc = cell_center;
    	double my_radius = curr_radius;
    	Coord neighbor_loc;
    	double neighbor_radius;
    	Coord rep_force = Coord(0,0);
    	Coord adh_force = Coord(0,0);
    	Coord adh_force_reg = Coord(0,0);
    	//distance between two cell centers
    	double d_ij;
    	//vector for direction of force
    	Coord v_ij;
    	double E_ij_inverse = 1/((3.0/4.0)*(2*(1-pow(POISSON,2))/ELASTIC_MOD));
    	double sqrt_term;
    	shared_ptr<Cell> this_cell = shared_from_this();
    	if(my_neighbor != this_cell){
		neighbor_loc = my_neighbor->get_cell_center();
		neighbor_radius = my_neighbor->get_curr_radius();
		d_ij = (my_loc-neighbor_loc).length();
		v_ij = (my_loc - neighbor_loc);
		sqrt_term = sqrt((my_radius*neighbor_radius)/(my_radius+neighbor_radius));
		if(Budding_On == 1){
			if((my_neighbor == curr_bud)){
                		if(this->check_has_bud()){
					//mother bud force
                        		adh_force += v_ij*-1*K_ADH*(d_ij-(my_radius+neighbor_radius));
                    		}
                    	}
			else if((my_neighbor == mother)){
                       		if(this->check_is_bud()){
                        		//bud mother force
                        		adh_force += v_ij*-1*K_ADH*(d_ij-(my_radius+neighbor_radius));
					
				}
			}
		}	
		if(my_radius + neighbor_radius - d_ij >= 0){
			rep_force += v_ij*.5*(1/d_ij)*pow(my_radius+neighbor_radius -d_ij,1.5)*E_ij_inverse*sqrt_term;
		}
		if((my_radius + neighbor_radius - d_ij < 1) && (my_radius + neighbor_radius - d_ij > -1)){
			adh_force_reg += v_ij*SINGLE_BOND_BIND_ENERGY*RECEPTOR_SURF_DENSITY*KB*TEMPERATURE*M_PI*(my_radius + neighbor_radius)*.5*-1;
		}		
		
		
	
     }
     curr_force = rep_force + adh_force + adh_force_reg;
     return curr_force;
}//MBK
void Cell::update_my_location(){
	this->cell_center = this->cell_center + this->curr_force*(1.0/(ETA*(1.0 + this->curr_radius)))*dt;
    	return;
}//MBK
void Cell::print_location_file_format(ofstream& ofs){
	ofs << rank << " " << cell_center.get_X() << " " << cell_center.get_Y() << " " << curr_radius << " " << separation_radius << " " << max_radius << " " << curr_protein << " " << my_OG_G1_length << " " << my_OG_bud_phase_length << endl;
	return;
}//MBK
void Cell::print_cell_cycle_file_format(ofstream& ofs){
	int mother_bool;
	int daughter_bool;
	int has_bud_bool;
	int is_bud_bool;
	if(is_mother){
		mother_bool = 1;
	}else{
		mother_bool = 0;
	}
	if(is_daughter){
		daughter_bool = 1;
	}else{
		daughter_bool = 0;
	}
	if(has_bud){
		has_bud_bool = 1;
	}else{
		has_bud_bool = 0;
	}
	if(is_bud){
		is_bud_bool = 1;
	}else{
		is_bud_bool = 0;
	}
	int mother_rank;
	int bud_rank;
	if(mother){
		mother_rank = this->mother->get_rank();
	}else{
		mother_rank = -1;
	}
	if(curr_bud){
		bud_rank = this->curr_bud->get_rank();
	}else{
		bud_rank = -1;
	}
	ofs << rank << " " << curr_phase << " " << mother_bool << " " << daughter_bool << " " << has_bud_bool << " " << is_bud_bool << " " << successor_number << " " << division_age << " " << time_born << " " << mother_rank << " " << bud_rank << " ";
	if(lineage_vec.size() > 0){
		for(unsigned int i = 0; i < this->lineage_vec.size();i++){
			ofs << "/" << lineage_vec.at(i);
		}
	}
	ofs << " ";
	if(this->daughter_vec.size() >0){
		for(unsigned int i = 0; i < this->daughter_vec.size();i++){
			ofs << "/" << daughter_vec.at(i);
		}
	ofs << " ";
	}
	if(div_sites_vec.size() > 0){
	for(unsigned int i = 0; i < this->div_sites_vec.size();i++){
		ofs << "/" << div_sites_vec.at(i);
	}
	}
	ofs << endl;
	return;
}
//***END OF FILE***
//**************************************************************************************************************************
