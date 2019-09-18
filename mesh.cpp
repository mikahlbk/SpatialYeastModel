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
void Mesh::make_mesh_pts(int start, int end, double increment){
	shared_ptr<Mesh> this_mesh = shared_from_this();
	double x_coord;
	double y_coord;
	int index = 0;
	for (int x = start; x< end+1;x = x+increment){
		x_coord = x;
		for(int y = start; y < end+1; y = y+increment){
			y_coord = y;	
			auto new_mesh_pt= make_shared<Mesh_Pt>(this_mesh, x_coord, y_coord, index);
			update_mesh_pts_vec(new_mesh_pt);
			index++;
		}
	}
	//cout << "How many mesh points made?" << mesh_pts.size() << endl;
	//for(unsigned int i = 0; i< mesh_pts.size(); i++){
	//	cout << mesh_pts.at(i)->get_index() << endl;
	//}
	return;
}
void Mesh::assign_neighbors(){
	for(unsigned int i = 0; i<mesh_pts.size();i++){
		mesh_pts.at(i)->find_neighbors();
	}
	return;
}
void Mesh::get_mesh_pts_vec(vector<shared_ptr<Mesh_Pt>>& curr_mesh_pts){
	curr_mesh_pts = this->mesh_pts;
	return;
}
void Mesh::update_mesh_pts_vec(shared_ptr<Mesh_Pt>& new_mesh_pt){
	mesh_pts.push_back(new_mesh_pt);
	return;
}
void Mesh::give_cell_to_bin(int& index, shared_ptr<Cell>& new_cell){
	this->mesh_pts.at(index)->add_cell(new_cell);
	return;
}
void Mesh::get_cells_from_bin(int& index, vector<shared_ptr<Cell>>& neighbors){
	//cout << "i think its an index problem" << endl;
	//cout << "num msh pts" << mesh_pts.size()<< endl;
	//cout << "index cell from bin function " << index << endl;
	vector<int> neighbor_bins;
	//cout << "this bin" << endl;
	this->mesh_pts.at(index)->get_cells(neighbors);
	//cout << "neighbor bins" << endl;
	this->mesh_pts.at(index)->get_neighboring_bins(neighbor_bins);
	for(unsigned int i = 0; i<neighbor_bins.size();i++){
		this->mesh_pts.at(neighbor_bins.at(i))->push_back_cells(neighbors);
	}
	//cout << "but it could be somthing else" << endl;
	//get neighbors
	//push back cells from neighbors
	return;

}

