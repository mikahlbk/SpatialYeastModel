//mesh_pt.h

//***************************************
//Include Guards
#ifndef _MESH_PT_H_INCLUDED_
#define _MESH_PT_H_INCLUDED_

//**************************************
//forward declarations
class Mesh;
//*************************************
//include dependencies
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <memory>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "parameters.h"
#include "coord.h"
//**************************************************
//mesh_pt class declaration
class Mesh_Pt: public enable_shared_from_this<Mesh_Pt>{
	private:
		vector<shared_ptr<Cell>> cells;
		shared_ptr<Mesh> my_mesh;
		Coord center;
		int index;
		vector<int> neighbors;
	public:
		//constructor
		Mesh_Pt(shared_ptr<Mesh> my_mesh, double x, double y, int index);
		void find_neighbors();
		Coord get_center(){return center;}
		int get_index(){return index;}
		void add_cell(shared_ptr<Cell>& new_cell);
		void clear_cells_vec();
		void get_cells(vector<shared_ptr<Cell>>& neighbors);
		void get_neighboring_bins(vector<int>& neighbors);
		void push_back_cells(vector<shared_ptr<Cell>>& neighbor_cells);
};


//end mesh_pt class
//**********************************************
#endif
