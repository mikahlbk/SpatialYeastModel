// colony.h

//******************************************
// Include Guards
#ifndef _COLONY_H_INCLUDED_
#define _COLONY_H_INCLUDED_

//******************************************
// forward declarations

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
#include "externs.h"
//******************************************
//COLONY Class Declaration

class Colony: public enable_shared_from_this<Colony>{
	private:
		shared_ptr<Mesh> my_mesh;
		mt19937 dist_generator;
		vector<shared_ptr<Cell>> cells;
		
	public:
		//constructor
		Colony(shared_ptr<Mesh> my_mesh, mt19937 gen);
        	//Colony(shared_ptr<Mesh> my_mesh);
		//make founder cell
        	//void make_founder_cell(string filename);
		void make_founder_cell();
		double uniform_random_real_number(double a, double b);
		//getters and setters
		void get_Cells(vector<shared_ptr<Cell>>& new_cells);
		void update_Colony_Cell_Vec(shared_ptr<Cell> new_cell);
		int get_Num_Cells();
		shared_ptr<Mesh> get_mesh(){return my_mesh;}
		//cell actions
		void find_bin();
		void pull_daughter();
		void grow_cells();
		void update_cell_cycles(int Ti);
        	void perform_budding(int Ti);
		void perform_mitosis(int Ti);
		void match_up();
		void update_locations();
		shared_ptr<Cell> return_cell(int cell_rank);
		void update_protein_concentration();
        	void print_vtk_file(ofstream& ofs);
        	void write_data(ofstream& ofs);
};

//*********************************************
//End of file

#endif

