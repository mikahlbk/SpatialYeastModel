//mesh.h

//*************************************
//include guards 
#ifndef _MESH_H_INCLUDED_
#define _MESH_H_INCLUDED_

//**************************************
//forward declarations

//**************************************
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
#include "mesh_pt.h"
//**************************************************
//mesh struct declaration

class Mesh: public enable_shared_from_this<Mesh>{
	private:
		vector<pair<shared_ptr<Mesh_Pt>,int>> mesh_pts;
	public:
		//constructor
		Mesh();
		void make_mesh_pts(double x_start, double y_start, int num_buckets, double increment);
		void update_mesh_pts_vec(shared_ptr<Mesh_Pt>& new_mesh_pt, int index);
		void get_mesh_pts_vec(vector<shared_ptr<Mesh_Pt>>& mesh_points);
		void assign_neighbors();
		void assign_cell_to_bin(int& index, shared_ptr<Cell>& new_cell);
		void get_cells_from_bin(int& index, vector<shared_ptr<Cell>>& neighbors);
		void calculate_nutrient_concentration();
		double get_nutrient_conc(int bin_id);
};

//end mesh class
//*******************************
#endif

