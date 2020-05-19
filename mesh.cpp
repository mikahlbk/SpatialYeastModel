//mesh.cpp

//****************************************
//include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <memory>
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "colony.h"
#include "mesh_pt.h"
#include "mesh.h"
using namespace std;
//****************************************
//Public member functions for mesh.cpp

//constructor
Mesh::Mesh(){
	
	return;
}
void Mesh::make_mesh_pts(double x_start, double y_start, int num_buckets, double increment){
	shared_ptr<Mesh> this_mesh = shared_from_this();
	double x_coord = x_start;
	double y_coord = y_start;
	int index = 1;
	for (int i = 0;i < num_buckets+1 ;i++){
		x_coord = x_start;
		for(int j = 0; j<num_buckets+1; j++){	
			auto new_mesh_pt= make_shared<Mesh_Pt>(this_mesh, x_coord, y_coord, index);
			//should be paired vector
			update_mesh_pts_vec(new_mesh_pt,index);
			x_coord = x_coord + increment;
			index++;
		}
		y_coord = y_coord - increment;
	}
	//cout << "How many mesh points made?" << mesh_pts.size() << endl;
	for(unsigned int i = 0; i< mesh_pts.size(); i++){
		//cout << "index: " << mesh_pts.at(i).second << endl;
		//cout << "center: " << mesh_pts.at(i).first->get_center() << endl;

	}
	return;

}
void Mesh::update_mesh_pts_vec(shared_ptr<Mesh_Pt>& new_mesh_pt, int index){
	mesh_pts.push_back(make_pair(new_mesh_pt,index));
	return;
}
void Mesh::get_mesh_pts_vec(vector<shared_ptr<Mesh_Pt>>& new_mesh_pts){
	for(unsigned int i = 0; i < mesh_pts.size(); i++){
		new_mesh_pts.push_back(this->mesh_pts.at(i).first);
	}
	return;
}
void Mesh::assign_neighbors(){
	for(unsigned int i = 0; i<mesh_pts.size();i++){
		//cout << "Mesh Pt " << i << "neighbors" << endl;
		mesh_pts.at(i).first->find_neighbor_bins();
	}

	return;
}

void Mesh::assign_cell_to_bin(int& index, shared_ptr<Cell>& new_cell){
	this->mesh_pts.at(index).first->add_cell(new_cell);
	return;
}
void Mesh::get_cells_from_bin(int& index, vector<shared_ptr<Cell>>& neighbors){
	//cout << "i think its an index problem" << endl;
	//cout << "num msh pts" << mesh_pts.size()<< endl;
	//cout << "index cell from bin function " << index << endl;
	//cout << "this bin" << endl;
	vector<shared_ptr<Mesh_Pt>> neighbor_bins;
	this->mesh_pts.at(index).first->get_cells(neighbors);
	//cout << "neighbor bins" << endl;
	this->mesh_pts.at(index).first->get_neighbor_bins(neighbor_bins);
	for(unsigned int i = 0; i<neighbor_bins.size();i++){
		neighbor_bins.at(i)->add_cells_to_neighbor_vec(neighbors);
	}
	//cout << "but it could be somthing else" << endl;
	//get neighbors
	//push back cells from neighbors
	return;

}
double Mesh::get_nutrient_conc(int bin_id){
	shared_ptr<Mesh_Pt> curr_mesh_pt = this->mesh_pts.at(bin_id).first;
	double conc = curr_mesh_pt ->get_nutrient_conc();
	return conc;
}
	
	
