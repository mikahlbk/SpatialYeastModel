//mesh_pt.cpp

//****************************************************
//include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <memory>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "colony.h"
#include "mesh.h"
#include "mesh_pt.h"
using namespace std;
//****************************************
//Public member functions for mesh_pts.cpp

//constructor
Mesh_Pt::Mesh_Pt(shared_ptr<Mesh>  my_mesh, double x, double y, int index){
	this->my_mesh = my_mesh;
	this->center = Coord(x, y);
	this->index = index;
	return;
}
void Mesh_Pt::find_neighbors(){
	vector<shared_ptr<Mesh_Pt>> my_mesh_pts;
	this->my_mesh->get_mesh_pts_vec(my_mesh_pts);
	for(unsigned int i = 0; i < my_mesh_pts.size(); i++){
		if((this->center - my_mesh_pts.at(i)->get_center()).length() < 4){
			this->neighbors.push_back(my_mesh_pts.at(i)->get_index());
		}
	}
	return;
}
void Mesh_Pt::add_cell(shared_ptr<Cell>& new_cell){
	this->cells.push_back(new_cell);
	return;
}
void Mesh_Pt::clear_cells_vec(){
	this->cells.clear();
	return;
}
void Mesh_Pt::get_cells(vector<shared_ptr<Cell>>& neighbors){
	neighbors = this->cells;
	return;
}
void Mesh_Pt::get_neighboring_bins(vector<int>& neighbors){
	neighbors = this->neighbors;
	return;
}
void Mesh_Pt::push_back_cells(vector<shared_ptr<Cell>>& neighbor_cells){
	for(unsigned int i = 0; i<cells.size();i++){
		neighbor_cells.push_back(cells.at(i));
	}
	return;
}
