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
void Colony::make_cells(int init_cells){
	shared_ptr<Colony> this_colony = shared_from_this();
	//make new cell(s)
	double max_radius;
	double init_radius;
	Coord center;
	int rank;
	//for(int i = 0; i< init_cells; i++){
		rank = cells.size();
		max_radius = static_cast<double>(rand() % 3 + 99)/(double)(100.0);
		init_radius = .8;
		center = Coord(0,0);
		//cout << "Assigned max radius" << max_radius << endl;
		
		auto new_cell = make_shared<Cell>(this_colony, rank, center, max_radius, init_radius);
		update_cell_vec(new_cell);
    	/*rank = cells.size();
		max_radius = static_cast<double>(rand() % 3 + 99)/(double)(100.0);
		init_radius = .2;
		center = Coord(.3,0);
		//cout << "Assigned max radius" << max_radius << endl;
		
		new_cell = make_shared<Cell>(this_colony, rank, center, max_radius, init_radius);
		update_cell_vec(new_cell);*/
	
	return;
}

void Colony::get_Cells(vector<shared_ptr<Cell> >& new_cells){
	new_cells = cells;
	return;
}

void Colony::update_cell_vec(shared_ptr<Cell> new_cell){
	cells.push_back(new_cell);
	return;
}

void Colony::grow_cells(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->grow_cell();
	}
	return;
}
int Colony::get_Num_Cells(){
    return cells.size();
}
void Colony::update_cell_cycles(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->update_cell_cycle();
	}
	return;
}

void Colony::perform_budding(){
	cout << "budding" << endl;
    cout << cells.size() << endl;
    //for(unsigned int i=0; i < 1; i++){
		cout << "calling truth" << endl;
        if(cells.at(0)->get_S()){
			cout << "entering budding" << endl;
            perform_budding();
            cout << "check two" << endl;
		}
	//}
	return;
}

void Colony::perform_mitosis(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cout << "mitosis check" << endl;
        if(cells.at(i)->get_M()){
			perform_mitosis();
		}
        cout << "tissue" << endl;
    }
	return;
}

void Colony::update_Locations(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->calc_forces();
	}

	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->update_location();
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
