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
#include "externs.h"
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
void Mesh::make_mesh_pts(double x_start, double y_start, unsigned int num_buckets, double increment){
	shared_ptr<Mesh> this_mesh = shared_from_this();
	double x_coord = x_start;
	double y_coord = y_start;
	for (unsigned int i = 0; i < num_buckets;i++){
		x_coord = x_start;
		for(unsigned int j = 0; j<num_buckets; j++){	
			auto new_mesh_pt= make_shared<Mesh_Pt>(this_mesh, x_coord, y_coord);
			mesh_pts[i][j] = new_mesh_pt;
			x_coord = x_coord + increment; 
		}
		y_coord = y_coord - increment;
	}
	this->assign_mesh_pt_neighbors();
	return;

}
void Mesh::assign_mesh_pt_neighbors(){
	//cout << "Step 1" << endl;
	for(unsigned int i = 1; i < 155; i++){
		for(unsigned int j = 1; j < 155; j++){
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j-1]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j+1]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j-1]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j+1]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j-1]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j]);
			mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j+1]);
		}
	}
	//cout << "Step 2" << endl;
	int i = 0;
	for(unsigned int j = 1; j< 155; j++){
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j-1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j+1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j-1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j+1]);
	}
	//cout << "Step 3" << endl;
	i = 155;
	for(unsigned int j = 1; j< 155; j++){
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j-1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j+1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j-1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j+1]);
	}	
	//cout << "Step 4" << endl;
	int j = 0;
	for(unsigned int i = 1; i < 155; i++){
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j+1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j+1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j+1]);
	}
	//cout << "Step 5" << endl;
	j = 155;
	for(unsigned int i = 1; i< 155; i++){
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j-1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j-1]);
		mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j-1]);
	}
	//cout << "Step 6" << endl;
	i = 0;
	j = 0;
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j+1]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j+1]);
	//cout << "Step 7" << endl;
	i = 0;
	j = 155;
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j-1]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j-1]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i+1][j]);
	//cout << "Step 8" << endl;
	i = 155;
	j = 0;
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j+1]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j+1]);
	//cout << "Step 9" << endl;
	i = 155;
	j = 155;
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i-1][j-1]);
	mesh_pts[i][j]->update_neighbor_list(mesh_pts[i][j-1]);
	return;
}
void Mesh::get_mesh_pts_vec(vector<shared_ptr<Mesh_Pt>>& new_mesh_pts){
	for(unsigned int i = 0; i < 156; i++){
		for(unsigned int j = 0; j < 156; j++){
			new_mesh_pts.push_back(this->mesh_pts[i][j]);
		}

	}

	return;
}
void Mesh::clear_mesh_pt_cells(){
	for(unsigned int i = 0; i < 156; i++){
		for(unsigned int j = 0; j < 156; j++){
			mesh_pts[i][j]->clear_my_cells();
		}

	}
	return;
}
void Mesh::update_mesh_pt_nutrients(){
	for(unsigned int i = 0; i < 156; i++){
		for(unsigned int j = 0; j < 156; j++){
			//cout << "Mesh Pt " << i << "neighbors" << endl;
			mesh_pts[i][j]->update_my_nutrient_concentration();
		}
	}
	return;
}
//void Mesh::update_mesh_pt_cells(int i, int j,shared_ptr<Cell> new_cell){
//	this->mesh_pts[i][j]->update_my_cells(new_cell);
//	return;
//}	
