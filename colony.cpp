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

//*****************************************
//Public Member Functions for Colony.cpp

//constructor
Colony::Colony() {
	return;
}
void Colony::make_founder_cell(){
	//pointer to tell founder cell what colony
    //it belongs to
    shared_ptr<Colony> this_colony = shared_from_this();
	//make founder cell
    //variables needed to 
    //feed to cell constructor
	double max_radius;
	double init_radius;
	Coord center;
	int rank = 0;
	max_radius = static_cast<double>(rand() % 3 + 99)/(double)(100.0);
	init_radius = .8;
	center = Coord(0,0);
	auto new_cell = make_shared<Cell>(this_colony, rank, center, max_radius, init_radius);
	update_Colony_Cell_Vec(new_cell);
    rank = cells.size();
    //auto new_cell_2 = make_shared<Cell>(this_colony,rank, center+Coord(2,2), max_radius, init_radius);
    //update_Colony_Cell_Vec(new_cell_2);
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

void Colony::perform_budding(){
    for(unsigned int i=0; i < cells.size(); i++){
        if(cells.at(i)->get_S()){
		    cells.at(i)->perform_budding();
        }
	}
	return;
}

void Colony::perform_mitosis(){
	for(unsigned int i = 0; i < cells.size(); i++){
	    if(cells.at(i)->get_M()){
		    cells.at(i)->perform_mitosis();
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
	double indica = 0;
    int counter = 0;
    double largest_value;
    for(unsigned int i = 0; i < cells.size(); i++){
	    cells.at(i)->calc_forces_chou();
        //cells.at(i)->calc_forces_jonsson();
        //cells.at(i)->calc_forces_exponential();
	    // cells.at(i)->lennard_jones_potential();
    }
    for(unsigned int i = 0; i < cells.size(); i++){
	    cells.at(i)->update_location();
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
