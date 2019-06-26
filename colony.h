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

#include "parameters.h"
#include "coord.h"
#include "cell.h"

//******************************************
//COLONY Class Declaration

class Colony: public enable_shared_from_this<Colony>{
	private:
		vector<shared_ptr<Cell> > cells;
		
	public:
		//constructor
		Colony();
        //make founder cell
        void make_founder_cell();
		
		//getters and setters
		void get_Cells(vector<shared_ptr<Cell>>& new_cells);
		void update_Colony_Cell_Vec(shared_ptr<Cell> new_cell);
		int get_Num_Cells();

		//cell actions
		void grow_cells();
		void update_cell_cycles();
        void perform_budding();
		void perform_mitosis();
		void update_locations();
		void update_protein_concentration();
        void print_vtk_file(ofstream& ofs);
        void write_data(ofstream& ofs);
};

//*********************************************
//End of file

#endif

