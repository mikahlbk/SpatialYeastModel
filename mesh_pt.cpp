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
#include "externs.h"
#include "coord.h"
#include "cell.h"
#include "colony.h"
#include "mesh.h"
#include "mesh_pt.h"
using namespace std;
//****************************************
//Public member functions for mesh_pts.cpp

//constructor
Mesh_Pt::Mesh_Pt(shared_ptr<Mesh>  my_mesh, double x, double y){
	this->my_mesh = my_mesh;
	this->center = Coord(x, y);
	this->nutrient_conc = 1;
	return;
}
void Mesh_Pt::get_my_neighbors(vector<shared_ptr<Mesh_Pt>>& neighbor_bins){
	neighbor_bins = this->my_neighbors;
	return;
}
void Mesh_Pt::get_my_cells(vector<shared_ptr<Cell>>& returned_cells){
	returned_cells = this->my_cells;
	return;
}
void Mesh_Pt::update_my_cells(shared_ptr<Cell> new_cell){
	this->my_cells.push_back(new_cell);
	return;
}
void Mesh_Pt::update_neighbor_list(shared_ptr<Mesh_Pt> neighbor){
	this->my_neighbors.push_back(neighbor);
	return;
}
void Mesh_Pt::clear_my_cells(){
	this->my_cells.clear();
	return;
}
void Mesh_Pt::update_my_nutrient_concentration(){
	this->nutrient_conc = this->nutrient_conc - Substrate_uptake_rate*this->my_cells.size()*dt;
	//cout << nutrient_conc << endl;
	if(this->nutrient_conc < 0){
		this->nutrient_conc = 0;
	}
	return;
}
void Mesh_Pt::get_neighbor_cells(vector<shared_ptr<Cell>>& neighbor_cells){
	for(unsigned int i = 0; i< my_cells.size(); i++){
		neighbor_cells.push_back(my_cells.at(i));
	}
	vector<shared_ptr<Cell>> other_cells;
	for(unsigned int i = 0; i < my_neighbors.size();i++){
		my_neighbors.at(i)->get_my_cells(other_cells);
		for(unsigned int j = 0; j < other_cells.size();j++){
			neighbor_cells.push_back(other_cells.at(j));
		}
	}
	return;
}
