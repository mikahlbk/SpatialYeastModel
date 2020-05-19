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
	this->nutrient_conc = 1;
	return;
}
void Mesh_Pt::find_neighbor_bins(){
	vector <pair<double,shared_ptr<Mesh_Pt>>> distances;
	vector <shared_ptr<Mesh_Pt>> mesh_pts;
	my_mesh->get_mesh_pts_vec(mesh_pts);
	shared_ptr<Mesh_Pt> this_mesh_pt = shared_from_this();

	for(unsigned int i = 0; i < mesh_pts.size();i++){
		if(mesh_pts.at(i)!= this_mesh_pt){
			distances.push_back(make_pair((mesh_pts.at(i)->get_center()-center).length(),mesh_pts.at(i)));
		}
	}
	sort(distances.begin(),distances.end());
	for(unsigned int i = 0; i<8; i++){
		this->neighbors.push_back(distances.at(i).second);
		//cout << distances.at(i).second->get_index() << endl;
	}
	return;
}

void Mesh_Pt::clear_cells_vec(){
	this->cells.clear();
	return;
}
void Mesh_Pt::get_cells(vector<shared_ptr<Cell>>&  neighbors){
	neighbors = this->cells;
	return;
}
void Mesh_Pt::get_neighbor_bins(vector<shared_ptr<Mesh_Pt>>& neighbor_bins){
	neighbor_bins = this->neighbors;
	return;
}
void Mesh_Pt::add_cells_to_neighbor_vec(vector<shared_ptr<Cell>>& neighbor_cells){
	for(unsigned int i =0; i < cells.size();i++){
		neighbor_cells.push_back(cells.at(i));
	}
	return;
}
void Mesh_Pt::add_cell(shared_ptr<Cell>& new_cell){
	this->cells.push_back(new_cell);
	return;
}
void Mesh_Pt::calculate_nutrient_concentration(){
	double curr_mass;
	double total_mass = 0;
	double multiplier; 
	for(unsigned int i = 0;i<cells.size();i++){
		curr_mass = M_PI*pow(cells.at(i)->get_radius(),2);
		total_mass += curr_mass;
	}
	multiplier = this->nutrient_conc -(1-total_mass/K_MASS)*this->nutrient_conc*dt;
	this->nutrient_conc = multiplier;
	return;
}
