// colony.h

//******************************************
// Include Guards
#ifndef _COLONY_H_INCLUDED_
#define _COLONY_H_INCLUDED_

//******************************************
// forward declarations
class Cell;
class Mesh;
class Mesh_Pt;
//******************************************
// Include Dependencies
#include <string>
#include <vector>
#include <fstream> 
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <random>

#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "mesh.h"
#include "mesh_pt.h"
#include "externs.h"
//******************************************
//COLONY Class Declaration

class Colony: public enable_shared_from_this<Colony>{
	private:
		shared_ptr<Mesh> my_mesh;
		default_random_engine my_engine;
		vector<shared_ptr<Cell>> my_cells;
		
	public:
		//make colony
		Colony(shared_ptr<Mesh> my_mesh);
		
		//make founder cell
		void   make_founder_cell(shared_ptr<Colony> this_colony, int Ti);
		
		//***Getters***
		default_random_engine get_engine(){return my_engine;}
		void get_colony_cell_vec(vector<shared_ptr<Cell>>& curr_cells){curr_cells = my_cells;}
		int  get_num_cells(){return my_cells.size();}
		shared_ptr<Mesh> get_my_mesh(){return my_mesh;}
		
		//***Setters**
		void update_colony_cell_vec(shared_ptr<Cell> new_cell);

		//cell actions
		double roll_my_real_dis(int upper_bound);
		void update_cell_bin_ids();
		void update_cell_phase_lengths();
		void update_cell_radii();
		void perform_bud_separation();
		void update_cell_cycle_phases();
        	void add_buds(int Ti);
		void update_cell_protein_concentrations();
		void update_cell_locations();
        	void update_mesh_nutrient_concentration();
		void write_data(ofstream& ofs1, ofstream& ofs2);
};

//*********************************************
//End of file

#endif

