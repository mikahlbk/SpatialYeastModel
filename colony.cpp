// colony.cpp

//****************************************
// Inlcude Dependencies
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
#include "mesh.h"
//*****************************************
//Public Member Functions for Colony.cpp

//constructor
Colony::Colony(shared_ptr<Mesh> new_mesh) {
	this->my_mesh = new_mesh;
	return;
}
void Colony::make_founder_cell(){
    //pointer to tell founder cell what colony
    //it belongs to
    shared_ptr<Colony> this_colony = shared_from_this();
    //make founder cell
    //variables needed to 
    //feed to cell constructor
	double new_max_radius;
	double init_radius;
	Coord center;
	int rank = 0;
	new_max_radius =static_cast<double>(rand() % 3 + 99)/(double)(100.0)*radius_average;
	init_radius = .2;
	center = Coord(0,0);
	auto new_cell = make_shared<Cell>(this_colony, rank, center, new_max_radius, init_radius);
	update_Colony_Cell_Vec(new_cell);
    	rank = cells.size();
    	new_cell->set_mother(new_cell);
    /*center = Coord(.5,0);
	auto new_cell2 = make_shared<Cell>(this_colony, rank, center, new_max_radius, init_radius);
	update_Colony_Cell_Vec(new_cell2);
    rank = cells.size();
    new_cell->set_mother(new_cell);
    center = Coord(-.5,0);
	auto new_cell3 = make_shared<Cell>(this_colony, rank, center, new_max_radius, init_radius);
	update_Colony_Cell_Vec(new_cell3);
    rank = cells.size();
    new_cell->set_mother(new_cell);
	center = Coord(0,2);
	auto new_cell4 = make_shared<Cell>(this_colony, rank, center, new_max_radius, init_radius);
	update_Colony_Cell_Vec(new_cell4);
    rank = cells.size();
    new_cell->set_mother(new_cell);
    center = Coord(0,-2);
	auto new_cell5 = make_shared<Cell>(this_colony, rank, center, new_max_radius, init_radius);
	update_Colony_Cell_Vec(new_cell5);
    rank = cells.size();
    new_cell->set_mother(new_cell);*/


    //return;
}

void Colony::get_Cells(vector<shared_ptr<Cell>>& new_cells){
	new_cells = cells;
	return;
}

void Colony::update_Colony_Cell_Vec(shared_ptr<Cell> new_cell){
	cells.push_back(new_cell);
	return;
}

int Colony::get_Num_Cells(){
        return cells.size();
}
void Colony::find_bin(){
	//cout << "error in find bin?" << endl;
	vector<shared_ptr<Mesh_Pt>> mesh_pts;
	//cout << "get mesh" << endl;
	my_mesh->get_mesh_pts_vec(mesh_pts);
	//cout << "mesh pts loop" << endl;
	for(unsigned int i = 0; i< mesh_pts.size();i++){
		mesh_pts.at(i)->clear_cells_vec();
	}
	//cout << "mesh pts end loop" << endl;
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->find_bin();
		//cout << "assigned id" << cells.at(i)->get_bin_id() << endl;
	}
	return;
}

void Colony::grow_cells(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->grow_cell();
	}
	return;
}

void Colony::update_cell_cycles(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->update_cell_cycle();
	}
	return;
}

void Colony::perform_budding(int Ti){
    for(unsigned int i=0; i < cells.size(); i++){
        if(cells.at(i)->get_S()){
		cells.at(i)->perform_budding(Ti);
        }
     }
	return;
}

void Colony::perform_mitosis(int Ti){
	for(unsigned int i = 0; i < cells.size(); i++){
	    if(cells.at(i)->get_M()){
		    cells.at(i)->perform_mitosis(Ti);
		}
    }
	return;
}
void Colony::update_protein_concentration(){
    for(unsigned int i=0; i< cells.size(); i++){
        //cout << "Protein before" << cells.at(i)->get_protein_conc() << endl;
        cells.at(i)->compute_protein_concentration();
		//cout << "Protein after" << cells.at(i)->get_protein_conc() << endl;
	    }
	    return;
	}
void Colony::update_locations(){
	//cout << "in colony" << endl;
	//cout << cells.size() << endl;
	for(unsigned int i = 0; i < cells.size(); i++){
		cout <<"cell: "<< i << endl;
		// cells.at(i)->calc_forces_chou();
		//cout << "update forces" << endl;
		cells.at(i)->calc_forces_Hertz();
		//cells.at(i)->calc_forces_jonsson();
		//cout << "forces updated" << endl;
		//cells.at(i)->calc_forces_exponential();
	        //cells.at(i)->lennard_jones_potential();
        }
        for(unsigned int i = 0; i < cells.size(); i++){
		//cout << "update locations" << endl;
	        cells.at(i)->update_location();
		//cout << "locations updated" << endl;
	}
    return;
}
void Colony::write_data(ofstream& ofs){
    ofs << cells.size() << endl;
    for(unsigned int i = 0; i < cells.size();i++){
        cells.at(i)->print_txt_file_format(ofs);
    }
    return;
}
void Colony::print_vtk_file(ofstream& ofs){
	vector<shared_ptr<Cell>> colony_cells;
	get_Cells(colony_cells);
	
	//every vtk file starts with this
	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Points representing cell centers and sizes for cell-center model of yeast" << endl;
	ofs << "ASCII" << endl << endl;

	//tell paraview the type of data to be read in
	ofs << "DATASET POLYDATA" << endl << endl;
	
	int num_cells = colony_cells.size();
	
	//tell paraview how many points this file is supposed
	//to represent, in our case we have one point per cell
	ofs <<"POINTS " << num_cells << " float64" << endl;
	for(unsigned int i = 0; i< colony_cells.size(); i++){
		//print each coordinate of each cell
		//center with an end line after
		colony_cells.at(i)->print_cell_center(ofs);
	}

	ofs << endl;

	//tell paraview to expect num_cells rows
	//followd by 2*numcells integers
	//the first integer says how many coordinaate points
	//define a vertex and the second integer gives the
	//index of the coordinate points from the above
	//list of all cell centers
	ofs << "Vertices " << num_cells << " " <<  2*num_cells << " " << endl;
	for(unsigned int i = 0; i<colony_cells.size();i++){

		ofs << 1 << " " << i << endl;	
	}
	
	ofs << endl;
	
    ofs << "POINT_DATA " << num_cells << endl;
    ofs << "SCALARS radius float64 " << 1 << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(unsigned int i = 0; i<colony_cells.size(); i++){
            ofs << colony_cells.at(i)->get_radius() << endl;
    }

    ofs << endl;
    return;
}
