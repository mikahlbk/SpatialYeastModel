// cell.h

//*********************************************************
// Include Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_

//*********************************************************
// forward declarations
class Colony;
class Mesh;
class Mesh_Pt;
//*********************************************************
// include dependencies
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <memory>
#include "parameters.h"
#include "coord.h"
#include "externs.h"
#include "colony.h"
#include "mesh.h"
#include "mesh_pt.h"
//***********************************************************
// Cell Class Declaration

class Cell: public enable_shared_from_this<Cell>{
	private:
		//variables for cell location
		Coord cell_center;
		Coord curr_force;
		shared_ptr<Mesh_Pt> bin_id;
		
		//variables for cell growth
		double curr_radius;
		double separation_radius;
 		double OG_separation_radius;
		double max_radius;
		double OG_max_radius;

		//variables for cell cycle
		int    curr_phase; //1 is G1, 2 is S/Ready for new Bud, 3 is G2
		double my_OG_G1_length;
		double my_G1_length;
		double my_G1_tracker;
		double my_OG_bud_phase_length;
		double my_bud_phase_length;
		bool   is_mother;
		bool   is_daughter;
		bool   has_bud;
		bool   is_bud;
		

		//variables for cell lineage
		shared_ptr<Colony> my_colony;
		int rank;
		int successor_number;
		int division_age;
		int time_born;
		shared_ptr<Cell> mother;
		shared_ptr<Cell> curr_bud;
		vector<int> lineage_vec;
		vector<int> daughter_vec;
		double curr_div_site;
		vector<double> div_sites_vec;
		vector<int> division_times_vec;
		//variables for protein dynamics
		double lambda; //replication rate
		double OG_lambda;
		double rho;    //transmission bias
		double curr_protein;
        	
	public:
		//Constructor for single founder cell
		Cell(shared_ptr<Colony> my_colony, Coord cell_center, int rank, double new_div_site,int Ti);
        	//Constructor for new bud
        	Cell(Coord cell_center, shared_ptr<Cell> mother, double new_div_site,int Ti);
		
		//***Getters***	
		
		//cell location and forces info
		Coord get_cell_center(){return cell_center;}
		Coord get_curr_force(){return curr_force;}
		shared_ptr<Mesh_Pt> get_bin_id(){return bin_id;}
		
		//cell growth info
 		double get_curr_radius(){return curr_radius;}
		double get_separation_radius(){return separation_radius;}
                double get_max_radius(){return max_radius;}
		
		//cell cycle info
		int    get_curr_phase(){return curr_phase;}
		double get_my_G1_length(){return my_G1_length;}
		double get_my_G1_tracker(){return my_G1_tracker;}
		double get_my_bud_phase_length(){return my_bud_phase_length;}
		bool   check_is_mother(){return is_mother;}
		bool   check_is_daughter(){return is_daughter;}
		bool   check_has_bud(){return has_bud;}
		bool   check_is_bud(){return is_bud;}

		//cell lineage info
		shared_ptr<Colony> get_my_colony(){return my_colony;}
		int get_rank(){return rank;}
		int get_successor_number(){return successor_number;}
		int get_division_age(){return division_age;}
		int get_time_born(){return time_born;}
		shared_ptr<Cell> get_mother(){return mother;}
		shared_ptr<Cell> get_curr_bud(){return curr_bud;}
	        void   get_lineage_vec(vector<int>& curr_lineage_vec);
		void   get_daughter_vec(vector<int>& daughter_vec);
		int    get_num_daughters(){return this->daughter_vec.size();}
		double get_curr_div_site(){return curr_div_site;}
		void   get_div_sites_vec(vector<double>& div_sites_vec);
                	
		//protein dynamics info
		double get_lambda(){return lambda;}
		double get_rho(){return rho;}
		double get_curr_protein(){return curr_protein;}

	        //***Setters***

		//cell location and forces info
		void set_cell_center(Coord new_center){this->cell_center = new_center;}
		void set_curr_force(Coord new_force){this->curr_force = new_force;}
		void set_bin_id(shared_ptr<Mesh_Pt> new_bin_id){this->bin_id = new_bin_id;}

	        //cell growth info
		void set_curr_radius(double new_radius){this->curr_radius = new_radius;}
		void set_separation_radius(double new_radius){this->separation_radius = new_radius;}
		void set_max_radius(double new_radius){this->max_radius = new_radius;}
		 
		//cell cycle info
		void set_curr_phase(int new_phase){this->curr_phase = new_phase;}
		void set_my_G1_length(double new_length){this->my_G1_length = new_length;}
		void set_my_G1_tracker(double new_tracker_val){this->my_G1_tracker = new_tracker_val;}
		void set_my_bud_phase_length(double new_length){this->my_bud_phase_length = new_length;}
		void set_is_mother(bool value){this->is_mother = value;}
		void set_is_daughter(bool value){this->is_daughter = value;}
		void set_has_bud(bool value){this->has_bud = value;}
		void set_is_bud(bool value){this->is_bud = value;}
		
		//cell lineage info
		void set_my_colony(shared_ptr<Colony> new_colony){this->my_colony = new_colony;}
		void set_rank(int new_rank){this->rank = new_rank;}
		void set_successor_number(int new_num){this->successor_number = new_num;}
		void update_division_age(int number){this->division_age = this->division_age + 1;}
		void set_division_age(int new_age){this->division_age = new_age;}
		void set_time_born(int time){this->time_born = time;}
		void set_mother(shared_ptr<Cell> new_mom){this->mother = new_mom;}
		void set_curr_bud(shared_ptr<Cell> new_bud){this->curr_bud = new_bud;}
		void update_lineage_vec(int new_lineage){this->lineage_vec.push_back(new_lineage);}
		void set_lineage_vec(vector<int> new_vec){this->lineage_vec = new_vec;}
		void set_curr_div_site(double new_site){this->curr_div_site = new_site;}
		void update_div_sites_vec(double new_site){this->div_sites_vec.push_back(new_site);}
		void set_div_sites_vec(vector<double> new_vec){this->div_sites_vec = new_vec;}

		//protein dynamics info
		void set_lambda(double new_lambda){this->lambda = new_lambda;}
		void set_rho(double new_rho){this->rho = new_rho;}
		void set_curr_protein(double new_protein){this->curr_protein = new_protein;}

		//***Other Functions in order of appearance in cell.cpp***
		shared_ptr<Mesh_Pt> find_my_nearest_bin();
		void update_my_phase_length();
		bool at_max_size();
		void update_my_radius();
		void perform_bud_separation();
		void update_my_phase();
		shared_ptr<Cell>  add_bud(int Ti);
		double determine_new_div_site();	
		shared_ptr<Cell> perform_budding(shared_ptr<Cell> mother_cell,double new_div_site,int Ti);	
		void update_my_protein_concentration();
		void compute_my_curr_force();
		Coord calc_forces_Hertz(shared_ptr<Cell> my_neighbor);	
		void update_my_location();
		//output
		void print_location_file_format(ofstream& ofs); 
		void print_cell_cycle_file_format(ofstream& ofs);
	
};

//***END OF FILE***
//*******************************************************************************************************
#endif


