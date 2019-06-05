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
		
		//getters and setters
		void get_Cells(vector<shared_ptr<Cell> >& new_cells);
		void update_cell_vec(shared_ptr<Cell> new_cell);
		void make_cells(int initcells);

		//cell actions
		void grow_cells();
		void update_cell_cycles();
        int get_Num_Cells();
        void perform_budding();
		void perform_mitosis();
		void update_Locations();
		void print_vtk_file(ofstream& ofs);
};

//*********************************************
//End of file

#endif

