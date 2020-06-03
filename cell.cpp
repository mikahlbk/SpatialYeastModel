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
#include <random>
#include <omp.h>
#include <boost/array.hpp>
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "colony.h"
#include "mesh_pt.h"
#include "mesh.h"
using namespace std;
//****************************************
//Public Member Functions for Cell.cpp

//Constructor for simulations starting with a single founder
Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center,double init_radius, double div_site){
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->max_radius = average_radius + average_radius*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0);
    this->curr_radius = max_radius;
    this->at_max_size = false;
    //curr_force set in function
    //bin_id set in function
    this->age = 0;
    this->T_age = 0;
    this->my_G1_length = average_G1_mother + average_G1_mother*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    this->my_G2_length = average_G2_mother + average_G2_mother*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0);
    //reset this after every mitosis because this determines the
    // length of time the next daughterwill spend on the mother cell
    this->G1 = true;
    this->G2 = false;
    this->S = false;
    this->M = false;
    this->growth_rate = max_radius/(my_G1_length);
    this->cell_cycle_increment = 1.0/(my_G1_length + my_G2_length);
    this->theoretical_cci = 1.0/(my_G1_length + my_G2_length);
    this->CP = .1;
    this->is_mother = false;
    this->has_bud = false;
    //curr_bud assigned in budding function
    //each new daughter is added to daughters vec in budding function
    this->curr_div_site = div_site;
    this->div_site_vec.push_back(curr_div_site);
    this->is_bud = false;
    //set in make founder function
    this->mother_rank = rank;  
    //lineage vec filled out in make founder function
    griesemer_lineage.push_back(0);
    this->sector = 0;
    
    //**** get rid of these ASAP****
    this->curr_protein = P_0;
    this->color = 0;
    //*************************
 
    return;
}
//Constructor for new daugher after division
Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center, double init_radius, shared_ptr<Cell> mother,int mother_rank, double div_site, vector<int> lineage, vector<int>g_lineage,int sector, int my_col, double g_two_from_mother){
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->curr_radius = init_radius;
    this->max_radius = average_radius + average_radius*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0);
 ;
    this->at_max_size = false;
    //curr_force set in function
    //bin id assigned in function
    this->age = 0;
    this->T_age = 0;
    this->my_G1_length = g_two_from_mother + extra_G1_daughter + extra_G1_daughter*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0);
    this->my_G2_length = average_G2_mother + average_G2_mother*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    this->G1 = true;
    this->G2 = false;
    this->S = false;
    this->M = false;
    this->growth_rate = max_radius/(my_G1_length + my_G2_length);
    this->cell_cycle_increment = 1.0/(my_G1_length + my_G2_length);
    this->theoretical_cci = 1.0/(my_G1_length + my_G2_length);
    this->CP = 0;
    this->is_mother = false;
    this->has_bud = false;
    //curr_bud assigned in function
    //daughters vector updated in function
    this->curr_div_site = div_site;
    this->div_site_vec.push_back(curr_div_site);
    this->is_bud = true;
    this->mother = mother;
    this->mother_rank = mother_rank;
    this->lineage = lineage;
    this->griesemer_lineage = g_lineage;
    this->sector = sector;
   
    //******get rid of these ASAP**********
    this->equi_point = cell_center;
    this->curr_protein = 0;
    this->at_max_size = false;
    this->color = my_col;
    //**************************************

    return;
}
//constructor for starting with multiple cells
/*Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center, double max_radius, double init_radius, double div_site, int bud_status,int phase, double cell_prog,int Mother, int my_color){
    //cout << "im in" << endl;
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->max_radius = max_radius;
    this->curr_radius = init_radius;
    this->age = 0;
    this->T_age = 0;
    this-> CP = cell_prog;
    if(phase == 1){
    	this->G1 = true;
	this->G2 = false;
	this->S = false;
	this->M = false;
    }
    else if(phase == 2){
    	this->G1 = false;
	this->G2 = true;
	this->S = false;
	this->M = false;
    }
    else if(phase == 3){

    	this->G1 = false;
	this->G2 = false;
	this->S = true;
	this->M = false;
    }
    else{

    	this->G1 = false;
	this->G2 = false;
	this->S = false;
	this->M = true;
    }
    this->G_one_length = G_one + G1*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    this->is_mother = false;
    this->mother_rank = Mother;
    this->mother = nullptr;
    this->curr_bud = nullptr;
    if(bud_status == 1){
	this->is_bud = true;
    }else{
	this->is_bud = false;
    }
    this->has_bud = false;
    this->equi_point = cell_center;
    this->curr_protein = P_0;
    this->curr_div_site = div_site;
    this->div_site_vec.push_back(curr_div_site);
    this->sector = 0;
    this->at_max_size = false;
    this->slow_grow = false;
    //this->lineage.push_back(1);
    //is vector
    this->color = my_color;
    griesemer_lineage.push_back(0);
    this->growth_rate = max_radius/((average_G1_daughter) + average_G1_daughter*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100));
    return;
}*/
//**********functions only to be used when starting colony with multiple cells*****************
/*void Cell::mother_bud_check(){
     shared_ptr<Cell> this_cell = shared_from_this();
     if(this->is_bud){
     	change_mother_vars(this->mother,this_cell);
	get_bud_status_mom(this->mother);
	this->color = mother->get_color();
	//cout << "bud rank: " <<  this_cell << endl;
     }	
     return;
}
void Cell::get_bud_status_mom(shared_ptr<Cell> mother){
	if(mother->return_has_bud()){
		cout << "true mother bud" << endl;
	}else {
		cout << "false mother bud" << endl;
	}
}
void Cell::mother_rank_to_ptr(){
	this->mother = this->my_colony->return_cell(mother_rank);
	return;
}
void Cell::change_mother_vars(shared_ptr<Cell> mother_cell, shared_ptr<Cell> bud){
     mother_cell->set_is_mother();
     mother_cell->set_has_bud();
     //cout << "bud: rank" << bud << endl;
     mother_cell->set_bud(bud);
     return;
}

void Cell::return_bud_status(){
	int my_rank = this->rank;
	int my_bud_rank;
	if(has_bud){
		my_bud_rank = this->curr_bud->get_rank();
	}else {
		my_bud_rank = 1000;
	}
	if(this->is_bud){
		cout << "bud"<< endl;
	}
	int phase = this->get_phase();
	cout << "rank: " << this->rank << "bud: " << my_bud_rank   << "phase: " << phase << endl;
}*/
//*************************************************************************************************

//***getters&setters that need a function,in order of cell.h***
void Cell::set_mother(shared_ptr<Cell> mother){
     this->mother = mother;
     return;
}
void Cell::get_daughters_vec(vector<shared_ptr<Cell>>& curr_daughters){
     shared_ptr<Cell> this_cell = shared_from_this();
     curr_daughters = this_cell->daughters;
     return;
}
void Cell::get_div_site_vec(vector<double>& previous_div_sites){
     shared_ptr<Cell> this_cell = shared_from_this();
     previous_div_sites = this_cell->div_site_vec;
     return;
}
void Cell::get_lineage_vec(vector<int>& curr_lineage_vec){
     shared_ptr<Cell> this_cell = shared_from_this();
     curr_lineage_vec = this_cell->lineage;
     return;
}
void Cell::update_lineage_vec(int cell_rank){
     this->lineage.push_back(cell_rank);
     return;
}
void Cell::get_griesemer_lineage_vec(vector<int>& curr_griesemer_lineage){
     shared_ptr<Cell> this_cell = shared_from_this();
     curr_griesemer_lineage = this_cell->griesemer_lineage;
     return;
}
int Cell::get_phase(){
	if(this->G1){
		return 1;
	}else if(this->G2){
		return 2;
	}else if(this->S){
		return 3;
	}else if(this->M){
		return 4;
	}
}
//****functions in order of cell.h***
void Cell::find_bin(){
     shared_ptr<Cell> this_cell = shared_from_this();	
     vector<shared_ptr<Mesh_Pt>> mesh_pts;
     this->get_colony()->get_mesh()->get_mesh_pts_vec(mesh_pts);
     double smallest_dist=1000000;
     int smallest_index;
     double curr_dist;
     //cout << "smallest index before starting" << smallest_index << endl;
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
     this->my_colony->get_mesh()->assign_cell_to_bin(smallest_index,this_cell);
     this->bin_id = smallest_index;
     return;
}
void Cell::update_growth_rate(){
     double multiplier = this->get_nutrient_conc(bin_id);
     if(multiplier > .1){
	this->cell_cycle_increment = this->theoretical_cci;
     }else{
	this->cell_cycle_increment = this->theoretical_cci*multiplier;
     }
     //cout << "Nutrient Conc: " << multiplier << "Cell cycle: " << cell_cycle_increment << endl; 
     return;
}
double Cell::get_nutrient_conc(int curr_bin_id){
     shared_ptr<Colony> my_colony= this->get_colony();
     shared_ptr<Mesh> my_mesh = my_colony->get_mesh();
     double nutrient_conc = my_mesh->get_nutrient_conc(curr_bin_id);
     return nutrient_conc;
}
void Cell::grow_cell(){
    if(!at_max_size){
     	//cout << "Rank " << rank << " size " << curr_radius << "set max size " << max_radius << endl;
	curr_radius = curr_radius + growth_rate*dt;
	//cout << "Rank: " << rank << "Curr radius: " << curr_radius << endl;
    }
    if(curr_radius >= max_radius){
	//cout << "made it" << endl;
	this->at_max_size = true;
    }
    return;
}
void Cell::update_cell_cycle(){
     this->T_age++;
     shared_ptr<Cell> this_cell = shared_from_this();
     if(this->G1 && this->CP >= .1 && this->age >0){
     	//cout << "cell cycle" << cell_cycle_increment << "actual " << CP << endl;
	this->G1 = false;
	this->S = true;
     }
     else if(this->G1 && this->CP >= 1 && this->age < 1){
	this->G1 = false;
	this->S = true;
     }
     if(is_bud && (curr_radius >= size_at_div)){
	this_cell->get_mother()->enter_mitosis();
     }
     //cout << "cell cycle" << cell_cycle_increment << "actual " << CP << endl;
     CP = CP + cell_cycle_increment*dt;	
     //cout << "cell cycle" << cell_cycle_increment << "actual " << CP << endl;
	
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
    //cout << "Rank: " << this->get_rank() << "Bud Formed: " << Ti << "Curr Radius: " << this->curr_radius << endl;
    //increment budding age
    this-> age = age+1;
    shared_ptr<Cell> this_cell = shared_from_this();
    shared_ptr<Colony> this_colony = this->get_colony();
    int daughter_cell_rank = this_colony->get_num_cells();
    double daughter_init_radius = 0;
    double mother_division_site;
    bool used = false;
    
    if(Division_Pattern==0){//axial
    	mother_division_site = this->curr_div_site;//axial same side
    }else if(Division_Pattern == 1){//bipolar
	mother_division_site = this->curr_div_site + M_PI;//bipolar opposite side
    }else if(Division_Pattern == 2){//mixed
	if(mother_status()){
		//50% chance axial or bipolar
		if(this->get_colony()->uniform_random_real_number(0.0,1.0)<=.5){
			mother_division_site = this->curr_div_site;
		}else{
			mother_division_site = this->curr_div_site + M_PI;
		}
	}else{
		mother_division_site = this->curr_div_site + M_PI;
	}
    }    
    //50% chance move up or down a little
    if(this->get_colony()->uniform_random_real_number(0.0,1.0)<=.5){
    	//move 10 degrees counter clockwise until new spot found
    	do{	
		used = false;
		mother_division_site = mother_division_site + DIV_SHIFT_RADIANS;
		for(unsigned int i = 0;i<this->div_site_vec.size();i++){
			if(mother_division_site == div_site_vec.at(i)){
				used = true;
			}	
		}	
     	}while(used);
     }else{
	//move 10 degress clockwise until new spot found
	do{
		used = false;
		mother_division_site = mother_division_site - DIV_SHIFT_RADIANS;
		for(unsigned int i = 0;i<this->div_site_vec.size();i++){
			if(mother_division_site == div_site_vec.at(i)){
				used = true;
			}
		}
	}while(used);
    }
    //update mother division site vector
    this->div_site_vec.push_back(mother_division_site);
    this->curr_div_site = mother_division_site;
    //cout << "rank " << this-> rank << " divsite " << division_site << endl;
    double new_center_x = this->cell_center.get_X()+(curr_radius+daughter_init_radius)*cos(curr_div_site);
    double new_center_y = this->cell_center.get_Y()+(curr_radius+daughter_init_radius)*sin(curr_div_site);
    Coord new_center = Coord(new_center_x,new_center_y);
    vector<int> new_lineage;
    this->get_lineage_vec(new_lineage);
    if(this->rank != 0){
    	new_lineage.push_back(this->rank);
    }
    vector<int> new_g_lineage;
    this->get_griesemer_lineage_vec(new_g_lineage);
    new_g_lineage.push_back(this->age);
    int sector;
    if(this->rank == 0){
    	sector = this-> age;
    }
    else{
    	sector = this->sector;
    }
    //cout << "New cell rank: " << new_rank << endl;
    //****new cell stuff***
    auto new_cell = make_shared<Cell>(this_colony, daughter_cell_rank, new_center, daughter_init_radius,this_cell,this->rank,mother_division_site+M_PI,new_lineage, new_g_lineage, sector,this->color,this->my_G2_length);
    new_cell->find_bin();
    //***mother cell stuff***
    this->G1 = false;
    this->S = false;
    this->G2 = true;
    this->M = false;
    this->is_mother = true;
    this->has_bud = true;
    this->curr_bud = new_cell;
    this->daughters.push_back(new_cell);
    this->div_site_vec.push_back(mother_division_site);
    //**colony stuff***
    this_colony->update_colony_cell_vec(new_cell);
    this->equi_point = Coord(curr_radius*cos(mother_division_site + M_PI/2),curr_radius*sin(mother_division_site + M_PI/2));
    return;
}
void Cell::perform_mitosis(int Ti){
    //separate mother and daughter
    //cout << "reset is bud" << endl;
    //cout << "mitosis time " << Ti << " " << this->get_rank() << endl;
    //cout << this->curr_bud << endl;
    double bud_prot = this->curr_protein*.4;
    double mother_prot = this->curr_protein*.6;
    //cout << "Mother before: " << mother_prot << " " << bud_prot << endl;
    this->set_protein_conc(this->curr_protein*.6);
    this->curr_bud->set_protein_conc(bud_prot);
    this->set_has_bud_to_false();
    this->curr_bud->set_is_bud_to_false(); 
    this->G1 = true;
    this->S = false;
    this->G2 = false;
    this->M = false;
    this->my_G1_length = average_G1_mother + average_G1_mother*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    this->my_G2_length = average_G2_mother + average_G2_mother*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    this->theoretical_cci = 1.0/(my_G1_length + my_G2_length);
    this->cell_cycle_increment = theoretical_cci;
    this->CP = 0;
    //daughter remains G1 and mother gets set to G1
    return;
}
void Cell::set_has_bud_to_false(){
    this->has_bud = false;
    return;
}
void Cell::set_is_bud_to_false(){
    this->is_bud = false;
    return;
}
void Cell::compute_protein_concentration(){
    if(!is_bud){
	this->curr_protein = curr_protein + r_LOGISTIC*curr_protein*(1-curr_protein/K_LOGISTIC)*(curr_protein/A_LOGISTIC -1)*dt;
    //this->curr_protein = curr_protein + r_LOGISTIC*curr_protein*(1-curr_protein/K_LOGISTIC);
    }
    return;
}
void Cell::set_protein_conc(double protein){
	this->curr_protein = protein;
	return;
}
void Cell::get_cell_force(){
    Coord force;
    vector<shared_ptr<Cell>> neighbor_cells;
    this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id,neighbor_cells);
    shared_ptr<Cell> this_cell = shared_from_this();
    #pragma omp parallel
    {
    	#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
	#pragma omp for reduction(+:force) schedule(static,1)
	for(unsigned int i = 0; i< neighbor_cells.size();i++){
    		if(neighbor_cells.at(i) != this_cell){
			//cout << "force on cell: " << i << endl;
			force += this_cell->calc_forces_Hertz(neighbor_cells.at(i));
		}
	}
     }
     this->curr_force = force;
     return;
}
Coord Cell::calc_forces_Hertz(shared_ptr<Cell> my_neighbor){
    //Describe force calculation in README!!!
    //vector<shared_ptr<Cell>> neighbor_cells;
    Coord my_loc = cell_center;
    double my_radius = curr_radius;
    Coord neighbor_loc;
    double neighbor_radius;
    Coord rep_force = Coord(0,0);
    Coord adh_force = Coord(0,0);
    Coord adh_force_reg = Coord(0,0);
    Coord bending_force = Coord(0,0);
    //distance between two cell centers
    double d_ij;
    //vector for direction of force
    Coord v_ij;
    double E_ij_inverse = 1/((3.0/4.0)*(2*(1-pow(POISSON,2))/ELASTIC_MOD));
    double sqrt_term;
    //bending force necessities
    //the angle made by my cell center, mom cell center and equilibrium point on mom membrane
    double curr_angle;
    //equilibrium angle is 90
    double eps = 0.0001;
    Coord mom_center = this->mother->get_cell_center();
    //Coord equi_point = this->mother->get_equi_point();
    shared_ptr<Cell> this_cell = shared_from_this();
    //for(unsigned int i = 0; i< neighbor_cells.size();i++){
    	if(my_neighbor != this_cell){
		neighbor_loc = my_neighbor->get_cell_center();
		neighbor_radius = my_neighbor->get_curr_radius();
		d_ij = (my_loc-neighbor_loc).length();
		v_ij = (my_loc - neighbor_loc);
		sqrt_term = sqrt((my_radius*neighbor_radius)/(my_radius+neighbor_radius));
		if(Budding_On == 1){
			//cout << "ADH_VARIABLE == " << ADH_ON << endl;
			//ADHESION_ON IS A BOOLEAN IN PARAMETERS FILE
			//turns in adhesion between mother and bud
			if((my_neighbor == curr_bud)){
                		if(this->has_bud){
                        		//cout << "Cell rank " << this->rank << " has daughter " << this->curr_bud->get_rank() << endl;
					//if(my_radius+neighbor_radius - d_ij < 0){
                        			adh_force += v_ij*-1*K_ADH*(d_ij-(my_radius+neighbor_radius));
					//}
                    		}
                    	}
			else if((my_neighbor == mother)){
                       		if(this->is_bud){
                        		//if(my_radius+neighbor_radius - d_ij < 0){
                        			adh_force += v_ij*-1*K_ADH*(d_ij-(my_radius+neighbor_radius));
						//curr_angle = this->compute_angle();
						//if(abs(curr_angle - pi/2) < eps){
							//do nothing theyre close
						//}
						/*else{
							double bend_constant = K_BEND*(curr_angle-pi/2)/(sqrt(1-pow(cos(curr_angle),2)));
							Coord me_to_mom = mom_center- cell_center;
							double me_to_mom_length = me_to_mom.length();
							Coord mom_to_equi_point = equi_point-mom_center;
							double mom_to_equi_point_length = mom_to_equi_point.length();
							Coord term1 = mom_to_equi_point/(mom_to_equi_point_length*me_to_mom_length);
							Coord term2 = me_to_mom*cos(curr_angle)/pow(me_to_mom_length,2);
							bending_force = (term1 + term2)*bend_constant;
						}*/
					//}
				}
			}
		}	
		if(my_radius + neighbor_radius - d_ij >= 0){
			rep_force += v_ij*.5*(1/d_ij)*pow(my_radius+neighbor_radius -d_ij,1.5)*E_ij_inverse*sqrt_term;
			//cout << "rep_force" << rep_force << endl;
		}
		if((my_radius + neighbor_radius - d_ij < 1) && (my_radius + neighbor_radius - d_ij > -1)){
			adh_force_reg += v_ij*SINGLE_BOND_BIND_ENERGY*RECEPTOR_SURF_DENSITY*KB*TEMPERATURE*M_PI*(my_radius + neighbor_radius)*.5*-1;
		}		
		
		
	
     }
     //cout << "Rep_force" << rep_force << endl;
     //Coord unit_vec = (rep_force + adh_force + adh_force_reg);
     //double length_vec =  (rep_force + adh_force + adh_force_reg).length();
     //cout  <<"unit vec" << unit_vec<< "length" << length_vec <<  endl;
     curr_force = rep_force + adh_force + adh_force_reg;
     //cout <<"curr force"  << curr_force << endl;// + adh_force;// + bending_force;
     return curr_force;
}
void Cell::update_location(){
    cell_center = cell_center + curr_force*1.0/(1+(curr_radius/2))*dt;
    return;
}
void Cell::print_txt_file_format(ofstream& ofs){
    ofs << rank << " " << cell_center.get_X() << " " << cell_center.get_Y() << " " << curr_radius << " ";
    for(unsigned int i = 0; i < this->griesemer_lineage.size();i++){
    	ofs << "/" << griesemer_lineage.at(i);
    }
    ofs << " ";
    //cout << "getting vec" << endl;
    vector<int> lineages;
    get_lineage_vec(lineages);
    //cout << "lineages loop" << endl;
    for(unsigned int i = 0; i < lineages.size();i++){
    ofs << "/" << lineages.at(i);
    }
    ofs << " " << this->get_sector() << " " << this->get_age() << " " << this->get_T_age() << " " << this->bud_status() << " " << this->get_phase() << " " << this->get_CP() << " " << this->mother->get_rank() << " " << this->curr_protein << " " << bin_id << endl;
    //ofs << " " << this->get_color() << endl;
    return;
}
//*******OLD******************
/*void Cell::slow_grow_on(){
	this->slow_grow = true;
}
void  Cell::return_lineage_vec(vector<int> & ranks){
	for(unsigned int i = 0; i< lineage.size(); i++){
		ranks.push_back(lineage.at(i)->get_rank());
	}
	return;
}
void Cell::change_gr(){
	this->G_one_length = 120;
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
int Cell::get_bud_status(){
	if(this->is_bud){
		return 1;
	}
	else{
		return 0;
	}
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
void Cell::set_has_bud(){
	this->has_bud = true;
}
void Cell::set_is_mother(){
	this->is_mother = true;	
	return;
}
void Cell::set_at_max_size(){
	this->at_max_size = true;
}

 
void Cell::grow_cell(){
	curr_radius = curr_radius + growth_rate*dt;
	return;
}

void Cell::update_cell_cycle(int Ti){
	if(curr_radius >= max_radius){
		S = true;
	}
	return;
}

void Cell::pull_daughter(){
	if(is_bud){
		//pulling force
		this->curr_force += Coord(-0,0);
	}else{
		this->curr_force += Coord(0,0);
	}
	return;
}
double Cell::compute_angle(){
	double angle;
	Coord left_vect = mother->get_equi_point() - mother->get_Cell_Center();
	Coord right_vect = cell_center - mother->get_Cell_Center();
	double left_len = left_vect.length();
	double right_len = right_vect.length();

	double costheta = left_vect.dot(right_vect)/(left_len*right_len);
	double theta = acos(min(max(costheta,-1.0),1.0));
	double crossProd = left_vect.cross(right_vect);
	if(crossProd < 0.0){
		theta = 2*pi - theta;
	}
	angle = theta;
	return angle;
}
bool Cell::far_enough_from_neighbors(){
    bool grow = false;
    double distance = 0;
    vector<shared_ptr<Cell>> neighbor_cells;
    this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id,neighbor_cells);
    shared_ptr<Cell> this_cell = shared_from_this();
    #pragma omp parallel
    {
    	#pragma omp for reduction(+:distance) schedule(static,1)
	for(unsigned int i = 0; i< neighbor_cells.size();i++){
    		if(neighbor_cells.at(i) != this_cell){
			//cout << "force on cell: " << i << endl;
			distance += this_cell->compute_distance(neighbor_cells.at(i));
		}
	}
     }
     if(distance >= .2){
     	grow = false;
     }
     else{
     	grow = true;
     }
     return grow;
}
double Cell::compute_distance(shared_ptr<Cell> neighbor_cell){
	double distance = 0;
	double d_ij = (cell_center - neighbor_cell->get_Cell_Center()).length();
	distance = curr_radius + neighbor_cell->get_radius() - d_ij;
	return distance;
}
void Cell::calc_forces_chou(){
    double d_ij;
    double delta_d;
    Coord v_ij;
    vector<shared_ptr<Cell>> neighbor_cells;
    //this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id,neighbor_cells);
	
    my_colony->get_Cells(neighbor_cells);
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
	//this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id,neighbor_cells);
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
    for(unsigned int i = 0; i < neighbor_cells.size(); i++)
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
 	
	this->curr_force = rep_force + adh_force;
	//cout << curr_force << endl;
	return;
}

void Cell::calc_forces_exponential(){
    double dd;
    //double aaa = .5;
    double ppp = .4;
	vector<shared_ptr<Cell>> neighbor_cells;
	my_colony->get_Cells(neighbor_cells);
	Coord my_loc = cell_center;
	double my_radius = curr_radius;
	Coord neighbor_loc;
	double neighbor_radius;
    double potential_double;
	Coord potential = Coord(0,0);
	//Coord adh_force = Coord(0,0);
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
                rep_force += (derivative*12*pow(inside,11) -derivative*12*pow(inside,5))*.01;//LJ_EPS;
        }
    }
    this->curr_force = rep_force;
    return;

}

void Cell::dnpm(const state_type &x, state_type& dxdt, double t){
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

}
void Cell::compute_protein_conc_DNPM(){
	state_type x = {10, 1, 1, 1};
	integrate(ref(*this), x, 0,10,.1);
	return;
}*/
	
