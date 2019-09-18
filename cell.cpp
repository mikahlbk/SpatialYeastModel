//cell.cpp

//******************************************
//Include Dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <memory>
#include <boost/array.hpp>
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "colony.h"
using namespace std;
//****************************************
//Public Member Functions for Cell.cpp

//constructor
Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center, double max_radius, double init_radius){
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->max_radius = max_radius;
    this->curr_radius = init_radius;
    this->age = 0;
    this->G1 = true;
    this->G2 = false;
    this->S = false;
    this->M = false;
    this->type = 0;
    this->mother = nullptr;
    this->curr_bud = nullptr;
    this->is_bud = false;
    this->has_bud = false;
    this->curr_protein = P_0;
    this->div_site1 = 0;
    this->div_site2 = div_site1;
    this->sector = 0;
    //this->lineage is vector
    griesemer_lineage.push_back(0);
    return;
}
Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center, double max_radius, double init_radius, shared_ptr<Cell> mother, double protein, double div_site, vector<shared_ptr<Cell>> lineage, vector<int>g_lineage,int sector, int bin_id){
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->max_radius = max_radius;
    this->curr_radius = init_radius;
    this->age = 0;
    this->type = 0;
    this->G1 = true;
    this->G2 = false;
    this->S = false;
    this->M = false;
    this->mother = mother;;
    this->curr_bud = nullptr;
    this->is_bud = true;
    this->has_bud = false;
    this->curr_protein = protein;
    this->div_site1 = div_site;
    this->div_site2 = div_site1;
    this->lineage = lineage;
    this->griesemer_lineage = g_lineage;
    this->sector = sector;
    this->bin_id = bin_id;
    return;
}
void Cell::get_daughters(vector<shared_ptr<Cell>>& daughter_cells){
	shared_ptr<Cell> this_cell = shared_from_this();
	daughter_cells = this_cell->daughters;
	return;
}
void Cell:: get_lineage_vec(vector<shared_ptr<Cell>>& new_lineages){
    new_lineages = this->lineage;
    return;
}
void Cell:: get_lineage_g_vec(vector<int>& new_g_lineages){
    new_g_lineages = this->griesemer_lineage;
    return;
}
void Cell::add_daughter(shared_ptr<Cell> daughter){
        this->daughters.push_back(daughter);
        return;
}
void Cell::find_bin(){
	shared_ptr<Cell> this_cell = shared_from_this();	
	vector<shared_ptr<Mesh_Pt>> mesh_pts;
	this->get_Colony()->get_mesh()->get_mesh_pts_vec(mesh_pts);
	double smallest_dist=100000;
	int smallest_index;
	double curr_dist;
	//cout << "Smallest index before starting" << smallest_index << endl;
	for(unsigned int i = 0; i< mesh_pts.size();i++){
		curr_dist = (this->cell_center - mesh_pts.at(i)->get_center()).length();
		if(curr_dist < smallest_dist){
			smallest_dist = curr_dist;
			smallest_index = mesh_pts.at(i)->get_index();
			//cout << "smallest index " << smallest_index << endl;
		}
	}
	//push cell back onto cell vector for that index
	//cout << "smallest index " << smallest_index << endl;
	my_colony->get_mesh()->give_cell_to_bin(smallest_index, this_cell);
	this->bin_id = smallest_index;
	return;
}
double Cell::calc_closest_center(){
    vector<shared_ptr<Cell>> neighbor_cells;
    my_colony->get_Cells(neighbor_cells);
    Coord my_loc = cell_center;
    shared_ptr<Cell> this_cell = shared_from_this();
    Coord neighbor_loc;
    double smallest_diff = 100;
    double diff_len;
    for(unsigned int i = 0; i<neighbor_cells.size(); i++){
        if(neighbor_cells.at(i) != this_cell){
            neighbor_loc = neighbor_cells.at(i)->get_Cell_Center();
            diff_len = (my_loc - neighbor_loc).length();
            if(diff_len < smallest_diff){
                smallest_diff = diff_len;
            }
        }
    }
    return smallest_diff;
}

void Cell::set_bud(shared_ptr<Cell> bud){
    this->curr_bud = bud;
    return;
}
void Cell::set_mother(shared_ptr<Cell> mother){
        this->mother = mother;
        return;
}
void Cell::reset_is_bud(){
        this->is_bud = false;
        return;
}
void Cell::reset_has_bud(){
        this->has_bud = false;
}
void Cell::grow_cell(){
	if(G1){
		double f_radius = k_g1*curr_radius + k_g2;
		double f_radius_max = k_g1*max_radius + k_g2;
		double slowing_factor = 1- (f_radius/f_radius_max);
		curr_radius = curr_radius + f_radius*.8*slowing_factor*dt;
    }
	return;
}
void Cell::update_cell_cycle(){
	shared_ptr<Cell> this_cell = shared_from_this();
	//cell cycle phase of a bud is started at G1
	//when the bud grows big enough
	//its cell cycle phase gets set to M
	//in addition the cell cycle phase of the mother
	//is set to M so that in the perform_mitosis 
	//function the mother cell and bud are both
	//given the status of regular cell 
	//without bud so the G1 phase is set to 
	//true for both cells
	if(is_bud){
		if(curr_radius > k_mitosis*max_radius){
			this_cell->get_mother()->enter_mitosis();
		}
	}	
	//G1 is set to true when cell object is made
	//if cell gets big enough and there
	//is not already a bud present it enters 
	//S phase so that a bud is created in the
	// perform_budding function. In this function, once 
	//the bud is made the cell cycle phase of the 
	//mother is set to G2 and the cell cycle of the daughter
	//is set to G1.
	else{
		if((curr_radius > k_G1*max_radius) && (!G2)){
			S = true;
		} 
	}

	return;
}

void Cell::enter_mitosis(){
    G1 = false;
	S = false;
	G2 = false;
	M = true;
	return;
}

void Cell::perform_budding(int Ti){
    //make bud
    //cout << this->get_rank() << "bud formed" << Ti <<endl;
    this-> age = age+1;
    shared_ptr<Cell> this_cell = shared_from_this();
    shared_ptr<Colony> this_colony = this->get_Colony();
    int new_rank = this_colony->get_Num_Cells();
    double init_radius = .2;//this->curr_radius*.2;
    double division_site;
    if(rand() % 100 < k_axial_frac*100){
        //axial same
        division_site = this->div_site2;
    }
    else{
            //polar opposite
            division_site = this->div_site2 + pi;
    }
    do{
        //move a little
        division_site = division_site + .01;
    }while((division_site == div_site2)&&(division_site == div_site1));
    this->div_site1 = this->div_site2;
    this->div_site2 = division_site;
    double new_max_radius = static_cast<double>(rand()%3 + 99)/(double)(100.0)*radius_average; 
    double new_center_x = this->cell_center.get_X()+curr_radius*cos(division_site);
    double new_center_y = this->cell_center.get_Y()+curr_radius*sin(division_site);
    Coord new_center = Coord(new_center_x, new_center_y);
    double protein = this->get_protein_conc()*.3;
    this->curr_protein = curr_protein*.7;
    vector<shared_ptr<Cell>> new_lineage;
    this->get_lineage_vec(new_lineage);
    new_lineage.push_back(this_cell);
    int sector;
    vector<int> new_lineage_g;
    this->get_lineage_g_vec(new_lineage_g);
    new_lineage_g.push_back(this->age);
    if(this->rank == 0){
    	sector = this-> age;
    }
    else{
    	sector = this->sector;
    }
	//cout << "New cell rank: " << new_rank << endl;
    auto new_cell = make_shared<Cell>(this_colony, new_rank, new_center,new_max_radius, init_radius,this_cell,protein, division_site+pi, new_lineage, new_lineage_g, sector, this->bin_id);
    my_colony->get_mesh()->give_cell_to_bin(this->bin_id,new_cell);
    this->add_daughter(new_cell);
    this->curr_bud = new_cell;
    this->has_bud = true;
    //this->curr_radius = curr_radius*.8;
    G1 = false;
    S = false;
    G2 = true;
    M = false;
    this_colony->update_Colony_Cell_Vec(new_cell);
	//new cell is G1 and old cell is G2
	return;
}

void Cell::perform_mitosis(int Ti){
	//separate mother and daughter
    //cout << "reset is bud" << endl;
    //cout << "mitosis time " << Ti << " " << this->get_rank() << endl;
    //cout << this->curr_bud << endl;
    this->curr_bud->reset_is_bud();
    //cout << "reset has bud" << endl;
    this->reset_has_bud(); 
    M = false;
    G1 = true;
    //daughter remains G1 and mother gets set to G1
	return;
}
void Cell::calc_forces_chou(){
    double d_ij; 
    double delta_d;
    Coord v_ij;
    vector<shared_ptr<Cell>> neighbor_cells;
    this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id,neighbor_cells);
	
    //my_colony->get_Cells(neighbor_cells);
    Coord my_loc = cell_center;
    double my_radius = curr_radius;
    Coord neighbor_loc;
    double neighbor_radius;
    Coord rep_force;
    Coord adh_force;
    shared_ptr<Cell> this_cell = shared_from_this();
    for(unsigned int i=0; i<neighbor_cells.size(); i++){
            if(neighbor_cells.at(i) != this_cell){
                neighbor_loc = neighbor_cells.at(i)->get_Cell_Center();
                neighbor_radius = neighbor_cells.at(i)->get_radius();
                d_ij = (neighbor_loc-my_loc).length();
                delta_d = my_radius + neighbor_radius - d_ij;
                v_ij = (neighbor_loc-my_loc)/d_ij;
                if(d_ij < (my_radius+neighbor_radius)){
                    rep_force += v_ij*delta_d*k_r*-1;
                    adh_force += v_ij*delta_d*k_a;
                }
            }       
    }
    this->curr_force = rep_force +  adh_force;
    return;
}
void Cell::calc_forces_jonsson(){
	vector<shared_ptr<Cell>> neighbor_cells;
	//put cells from that bin and neighbors  into the neighbor cells vector
	//cout << "get cells from bin" << endl;
	//cout << "this bin id " << this->bin_id << endl;
	this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id,neighbor_cells);
	//cout << "bad" << endl;
	//my_colony->get_Cells(neighbor_cells);
	Coord my_loc = cell_center;
	double my_radius = curr_radius;
	Coord neighbor_loc;
	double neighbor_radius;
	Coord rep_force;
	Coord adh_force;
	shared_ptr<Cell> this_cell = shared_from_this();
	Coord diff_vec;
    double diff_len;
    //cout << "for loop" << endl;
    for(unsigned int i = 0; i < neighbor_cells.size(); i++){
		if(neighbor_cells.at(i)!=this_cell){
            neighbor_loc = neighbor_cells.at(i)->get_Cell_Center();
		    neighbor_radius = neighbor_cells.at(i)->get_radius();
		    diff_len = (neighbor_loc - my_loc).length();
            diff_vec = (neighbor_loc - my_loc)/diff_len;
            if(diff_len < k_neighbor*(my_radius + neighbor_radius)){
            if((neighbor_cells.at(i) == curr_bud)){
                    if(this->has_bud){
                        //cout << "Cell rank " << this->rank << " has daughter " << this->curr_bud->get_rank() << endl;
                        rep_force += diff_vec*(diff_len-k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring;
                	    adh_force += diff_vec*(diff_len-k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring*k_adhesion_mother_bud;
                    }
                    else{
                         //cout << "regular cell" << rank << endl;
                        rep_force += diff_vec*(diff_len -k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring;
                        adh_force += diff_vec*(diff_len -k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring*k_adhesion_cell_cell;
                    }
                }
                else if((neighbor_cells.at(i) == mother)){
                       if(this->is_bud){
                            //cout << "Cell rank " << rank << " has mother " << mother->get_rank() << endl;
                	        rep_force += diff_vec*(diff_len-k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring;
	            	        adh_force += diff_vec*(diff_len-k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring*k_adhesion_mother_bud;  
                       }
                       else {
                             //cout << "regular cell" << rank << endl;
                            rep_force += diff_vec*(diff_len -k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring;
                            adh_force += diff_vec*(diff_len -k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring*k_adhesion_mother_daughter;
               
                       }
                }
                else {
		            //cout << "regular cell" << rank << endl;
                    rep_force += diff_vec*(diff_len -k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring;
                    adh_force += diff_vec*(diff_len -k_repulsion_cell_cell*(my_radius+neighbor_radius))*k_spring*k_adhesion_cell_cell;
               }
           }

	    }	
    }	
	this->curr_force = rep_force + adh_force;
	//cout << curr_force << endl;
	return;
}

void Cell::calc_forces_exponential(){
    double dd;
    double aaa = .5;
    double ppp = .4;
	vector<shared_ptr<Cell>> neighbor_cells;
	my_colony->get_Cells(neighbor_cells);
	Coord my_loc = cell_center;
	double my_radius = curr_radius;
	Coord neighbor_loc;
	double neighbor_radius;
    double potential_double;
	Coord potential = Coord(0,0);
	Coord adh_force = Coord(0,0);
	shared_ptr<Cell> this_cell = shared_from_this();
	Coord diff_vect;
    double diff_len;
    for(unsigned int i = 0; i < neighbor_cells.size();i++){
            if(neighbor_cells.at(i)!= this_cell){
                neighbor_loc = neighbor_cells.at(i)->get_Cell_Center();
                neighbor_radius = neighbor_cells.at(i)->get_radius();
                diff_vect = (my_loc-neighbor_loc);
                diff_len = diff_vect.length();
                dd = diff_len/((my_radius + neighbor_radius)*.5);
                if(dd < 2){
                    potential_double = exp(-dd*ppp)*dd;
                    potential = potential + (diff_vect/dd)*potential_double/ppp;
                }
            }
    }
    this->curr_force = potential;
    return;
}
void Cell::lennard_jones_potential(){
    vector<shared_ptr<Cell>> neighbor_cells;
    my_colony->get_Cells(neighbor_cells);
    Coord my_loc = cell_center;
    double my_radius = curr_radius;
    Coord neighbor_loc;
    double neighbor_radius;
    Coord rep_force;
    Coord adh_force;
    shared_ptr<Cell> this_cell = shared_from_this();
    Coord diff_vec;
    double diff_len;
    Coord derivative;
    double inside;
    for(unsigned int i = 0; i < neighbor_cells.size(); i++){
        if(neighbor_cells.at(i) != this_cell){
                neighbor_loc = neighbor_cells.at(i)->get_Cell_Center();
                neighbor_radius = neighbor_cells.at(i)->get_radius();
                diff_vec = neighbor_loc - my_loc;
                diff_len = diff_vec.length();
                inside = (my_radius + neighbor_radius)/diff_len;
                derivative = (diff_vec*(my_radius + neighbor_radius))/pow(diff_len,3);
                rep_force += (derivative*12*pow(inside,11) -derivative*12*pow(inside,5))*LJ_EPS;
        }
    }
    this->curr_force = rep_force;
    return;

}
void Cell::update_location(){
    //cout << this->curr_bud << endl;
    cell_center = cell_center + curr_force*dt;
	return;
}
/*void Cell::dnpm(const state_type &x, state_type& dxdt, double t){
	int fsum = 0;
	int size_x = x.size()-1;
	for(unsigned int i = 1; i<size_x; i++){
		fsum+= x[i];
	}
	int fsum_greater = 0;

	dxdt[0] = alpha -mu*x[0]-2*beta*x[0]*fsum + Gamma*n_0*(n_0-1)*fsum;
	for(unsigned int i = 1; i<size_x;i++){
		fsum_greater = 0;
		for(int j = i+1; j< size_x; j++){
			fsum_greater += x[j];
		}
		dxdt[i] = -2*beta*x[0]*(x[i]-x[i-1])-mu*x[i]-Gamma*(i-1)*x[i]+2*Gamma*fsum_greater;
	}
	return;

}*/
/*void Cell::compute_protein_conc_DNPM(){
	state_type x = {10, 1, 1, 1};
	integrate(ref(*this), x, 0,10,.1);
	return;
}*/
void Cell::compute_protein_concentration(){
    this->curr_protein = curr_protein + r_LOGISTIC*curr_protein*(1-curr_protein/K_LOGISTIC)*(curr_protein/A_LOGISTIC -1);
    //this->curr_protein = curr_protein + r_LOGISTIC*curr_protein*(1-curr_protein/K_LOGISTIC);
    return;
} 
void Cell::print_txt_file_format(ofstream& ofs){
    ofs << rank << " " << cell_center.get_X() << " " << cell_center.get_Y() << " " << curr_radius << " " << curr_protein << " ";
    for(unsigned int i = 0; i < this->griesemer_lineage.size();i++){
    	ofs << "/" << griesemer_lineage.at(i);
    }
    ofs << " " << this->get_sector() << endl;

    return;
}
void Cell::print_cell_center(ofstream& ofs){
    ofs << cell_center.get_X() << " " << cell_center.get_Y() << " " << 0 << endl;
    return;
}
	
