//mesh.h

//*************************************
//include guards 
#ifndef _MESH_H_INCLUDED_
#define _MESH_H_INCLUDED_

//**************************************
//forward declarations
class Colony;
class Cell;
class Mesh_Pt;
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
#include "externs.h"
#include "mesh_pt.h"
#include "cell.h"
#include "colony.h"
//**************************************************
//mesh struct declaration

class Mesh: public enable_shared_from_this<Mesh>{
	private:
		shared_ptr<Mesh_Pt> mesh_pts[156][156];
	public:
		//constructor
		Mesh();
		void make_mesh_pts(double x_start, double y_start, unsigned int num_buckets, double increment);
		void assign_mesh_pt_neighbors();
		void get_mesh_pts_vec(vector<shared_ptr<Mesh_Pt>>& mesh_points);
		void clear_mesh_pt_cells();
		void update_mesh_pt_nutrients();
		//void update_mesh_pt_cells(int i, int j,shared_ptr<Cell> new_cell);
};

//end mesh class
//*******************************
#endif

