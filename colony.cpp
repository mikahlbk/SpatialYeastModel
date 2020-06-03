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
#include <random>

#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "colony.h"
#include "mesh.h"
#include "mesh_pt.h"
//*****************************************
//Public Member Functions for Colony.cpp

//constructor
Colony::Colony(shared_ptr<Mesh> new_mesh, mt19937  gen) {
	this->my_mesh = new_mesh;
	this->dist_generator = gen;
	return;
}
void Colony::make_founder_cell(){
     shared_ptr<Colony> this_colony = shared_from_this();
     int num_cells = 0;
     //make founder cell
     //variables needed to 
     //feed to cell constructor
     double init_radius;
     double div_site = 2*M_PI*(this->uniform_random_real_number(0.0,1.0));
     Coord center;
     int rank = 0;
     init_radius = 0;
     center = Coord(0,0);
     auto new_cell = make_shared<Cell>(this_colony, rank, center, init_radius, div_site);
     update_colony_cell_vec(new_cell);
     new_cell->set_mother(new_cell);
     new_cell->update_lineage_vec(rank); 
    return;
}
/*void Colony:: make_founder_cell(string filename){
	shared_ptr<Colony> this_colony = shared_from_this();
	int num_cells = 0;
	ifstream ifs(filename.c_str());
	if (!ifs){
		cout << "No File" << endl;
		return;
	}
	stringstream ss;
	string line;
	string item;
	int rank;
	double x_coord;
	double y_coord;
	double radius;
	int bud_status;
	int phase;
	double cell_prog;
	int mother;
	while(getline(ifs,line)){
		ss.str(line);
		getline(ss,item,' ');
		rank = stoi(item);
		getline(ss,item,' ');
		x_coord = stod(item);
		getline(ss,item,' ');
		y_coord = stod(item);
		getline(ss,item,' ');
		radius = stod(item);
		getline(ss,item,' ');
		bud_status = stoi(item);
		getline(ss,item,' ');
		phase = stoi(item);
		getline(ss,item,' ' );
		cell_prog = stod(item);
		getline(ss,item,' ');
		mother = stoi(item);
		double new_max_radius;
     		double init_radius;
     		uniform_real_distribution<> div_dist(0.0,1.0);
     		double color_roll = div_dist(this->dist_generator);
		int color;
		//int START_CELL = rand() % 11;
		if(rank == 1){//START_CELL){
			color = 0;
		}else{
			color = 1;
		}
		//uniform_real_distribution<> div_dist(0.0,1.0);
     		double div_site = 2*pi*div_dist(this->dist_generator);
     		Coord center;
     		new_max_radius = average_radius;// + static_cast<double>(rand() % 10 + 99)/(double)(100.0);
     		init_radius = radius;
     		center = Coord(x_coord,y_coord);
     		cout << rank << bud_status << cell_prog << endl;
		auto new_cell = make_shared<Cell>(this_colony, rank, center, new_max_radius, init_radius, div_site,bud_status,phase,cell_prog,mother,color);
     		update_Colony_Cell_Vec(new_cell);
		new_cell->set_mother(new_cell);
     		new_cell->update_lineage_vec(new_cell);
   
		ss.clear();
	}
	return;
}*/
double Colony::uniform_random_real_number(double a, double b){
	uniform_real_distribution<> dis(a,b);
	return dis(this->dist_generator);
}
/*void Colony::match_up(){
	for(unsigned int i = 0; i < cells.size(); i++){
		//cout << "rank: " << i << "pointer: " << cells.at(i) << endl;
		cells.at(i)->mother_rank_to_ptr();
		cells.at(i)->mother_bud_check();

	}
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->return_bud_status();
	}
	return;
}*/

void Colony::get_colony_cell_vec(vector<shared_ptr<Cell>>& curr_cells){
	curr_cells = my_cells;
	return;
}
void Colony::update_colony_cell_vec(shared_ptr<Cell> new_cell){
	my_cells.push_back(new_cell);
	return;
}
int Colony::get_num_cells(){
        return my_cells.size();
}
void Colony::find_bin(){
	//cout << "error in find bin?" << endl;
	vector<shared_ptr<Mesh_Pt>> mesh_pts;
	//cout << "get mesh" << endl;
	my_mesh->get_mesh_pts_vec(mesh_pts);
	//cout << "mesh pts loop" << endl;
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i< mesh_pts.size();i++){
		mesh_pts.at(i)->clear_cells_vec();
	}
	//cout << "mesh pts end loop" << endl;
	//#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_cells.at(i)->find_bin();
		//cout << "assigned id" << cells.at(i)->get_bin_id() <<" rank: " << cells.at(i)->get_rank() << endl;
	}
	return;
}
void Colony::update_growth_rates(){
	//cout << "error in find bin?" << endl;
	vector<shared_ptr<Mesh_Pt>> mesh_pts;
	//cout << "get mesh" << endl;
	my_mesh->get_mesh_pts_vec(mesh_pts);
	//cout << "mesh pts loop" << endl;
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i< mesh_pts.size();i++){
		mesh_pts.at(i)->calculate_nutrient_concentration();
	}
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i< my_cells.size();i++){
		my_cells.at(i)->update_growth_rate();
	}
	return;
}
void Colony::grow_cells(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_cells.at(i)->grow_cell();
	}
	return;
}

void Colony::update_cell_cycles(int Ti){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
		my_cells.at(i)->update_cell_cycle();
	}
	return;
}

void Colony::perform_budding(int Ti){
    //#pragma omp parallel for schedule(static,1)	
    for(unsigned int i=0; i < my_cells.size(); i++){
        if(my_cells.at(i)->is_S()){
		my_cells.at(i)->perform_budding(Ti);
        }
     }
	return;
}
/*void Colony::pull_daughter(){
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->pull_daughter();
	}
	return;
}*/
void Colony::perform_mitosis(int Ti){
	//#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < my_cells.size(); i++){
	    if(my_cells.at(i)->is_M()){
		    my_cells.at(i)->perform_mitosis(Ti);
		}
    }
	return;
}
void Colony::update_protein_concentration(){
   #pragma omp parallel for schedule(static,1) 
   for(unsigned int i=0; i< my_cells.size(); i++){
        //cout << "Protein before" << cells.at(i)->get_protein_conc() << endl;
        my_cells.at(i)->compute_protein_concentration();
	//cout << "Protein after" << cells.at(i)->get_protein_conc() << endl;
    }
	return;
}
/*shared_ptr<Cell> Colony::return_cell(int cell_rank){
	return this->cells.at(cell_rank);
}*/
void Colony::update_locations(){
	//cout << "in colony" << endl;
	//cout << cells.size() << endl;
	//double force_check = 0;
	//int counter = 0;
	//do {
	//force_check = 0;
//	#pragma omp parallel 
//	{
		//#pragma omp for reduction(+:force_check) schedule(static,1)
		#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < my_cells.size(); i++){
			//cout <<"cell: "<< i << endl;
			//cells.at(i)->calc_forces_chou();
			//cout << "update forces" << i <<  endl;
			my_cells.at(i)->get_cell_force();
			//force_check+= cells.at(i)->get_curr_force().length();
			//counter++;
			//cells.at(i)->pull_daughter();
			//cout << "forces updated" << i << endl;
			//cells.at(i)->calc_forces_jonsson();
			//cout << "forces updated" << endl;
			//cells.at(i)->calc_forces_exponential();
	        	//cells.at(i)->lennard_jones_potential();
        	}
//	}
	#pragma omp parallel for schedule(static,1)
        for(unsigned int i = 0; i < my_cells.size(); i++){
		//cout << "update locations" << endl;
	        my_cells.at(i)->update_location();
		//cout << "locations updated" << endl;
	}
	//} while((force_check > 100));
    return;
}
void Colony::write_data(ofstream& ofs){
    ofs << my_cells.size() << endl;
    for(unsigned int i = 0; i < my_cells.size();i++){
        my_cells.at(i)->print_txt_file_format(ofs);
    }
    return;
}
/*void Colony::print_vtk_file(ofstream& ofs){
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
}*/
