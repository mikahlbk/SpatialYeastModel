//mesh_pt.h

//***************************************
//Include Guards
#ifndef _MESH_PT_H_INCLUDED_
#define _MESH_PT_H_INCLUDED_

//**************************************
//forward declarations
class Colony;
class Mesh;
class Cell;
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
#include "externs.h"
#include "colony.h"
#include "mesh.h"
#include "cell.h"
//**************************************************
//mesh_pt class declaration
class Mesh_Pt: public enable_shared_from_this<Mesh_Pt>{
	private:
		shared_ptr<Mesh> my_mesh;
		Coord center;
		double nutrient_conc;
		vector<shared_ptr<Cell>> my_cells;
		vector<shared_ptr<Mesh_Pt>> my_neighbors;
	public:
		//constructor
		Mesh_Pt(shared_ptr<Mesh> my_mesh, double x, double y);
		Coord get_center(){return center;};
		void update_neighbor_list(shared_ptr<Mesh_Pt> neighbor);
		void get_my_neighbors(vector<shared_ptr<Mesh_Pt>>& neighbor_bins);
		void update_my_nutrient_concentration();
		double get_my_nutrient_conc(){return nutrient_conc;}
		void get_my_cells(vector<shared_ptr<Cell>>& returned_cells);
		void update_my_cells(shared_ptr<Cell> new_cell);
		void clear_my_cells();
		void get_neighbor_cells(vector<shared_ptr<Cell>>& neighbor_cells);
};


//end mesh_pt class
//**********************************************
#endif
