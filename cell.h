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
		double max_radius;
		double OG_max_radius;

		//variables for cell cycle
		int    curr_phase; //1 is G1, 0 is Mitosis
		double my_OG_G1_length;
		double my_G1_length;

		//variables for cell lineage
		shared_ptr<Colony> my_colony;
		int rank;

	public:
		//Constructor for single founder cell
		Cell(shared_ptr<Colony> my_colony, Coord cell_center, int rank,int Ti);
        //Constructor for new bud
        Cell(Coord cell_center,double init_radius,int Ti);
		
		//***Getters***	
		
		//cell location and forces info
		Coord get_cell_center(){return cell_center;}
		Coord get_curr_force(){return curr_force;}
		shared_ptr<Mesh_Pt> get_bin_id(){return bin_id;}
		
		//cell growth info
 		double get_curr_radius(){return curr_radius;}
        double get_max_radius(){return max_radius;}
		
		//cell cycle info
		int    get_curr_phase(){return curr_phase;}
		double get_my_G1_length(){return my_G1_length;}
		double get_my_G1_tracker(){return my_G1_tracker;}
	

		//cell lineage info
		shared_ptr<Colony> get_my_colony(){return my_colony;}
		int get_rank(){return rank;}
	
        //***Setters***

		//cell location and forces info
		void set_cell_center(Coord new_center){this->cell_center = new_center;}
		void set_curr_force(Coord new_force){this->curr_force = new_force;}
		void set_bin_id(shared_ptr<Mesh_Pt> new_bin_id){this->bin_id = new_bin_id;}

	    //cell growth info
		void set_curr_radius(double new_radius){this->curr_radius = new_radius;}
		void set_max_radius(double new_radius){this->max_radius = new_radius;}
		 
		//cell cycle info
		void set_curr_phase(int new_phase){this->curr_phase = new_phase;}
		void set_my_G1_length(double new_length){this->my_G1_length = new_length;}
		//void set_my_G1_tracker(double new_tracker_val){this->my_G1_tracker = new_tracker_val;}
	
		//cell lineage info
		void set_my_colony(shared_ptr<Colony> new_colony){this->my_colony = new_colony;}
		void set_rank(int new_rank){this->rank = new_rank;}
		
        //***Other Functions in order of appearance in cell.cpp***
		shared_ptr<Mesh_Pt> find_my_nearest_bin();
		void update_my_phase_length();
		bool at_max_size();
		void update_my_radius();
		void perform_division();
		void compute_my_curr_force();
		Coord calc_forces_linear_spring(shared_ptr<Cell> my_neighbor);	
		void update_my_location();
		//output
		void print_location_file_format(ofstream& ofs); 
	
};

//***END OF FILE***
//*******************************************************************************************************
#endif


