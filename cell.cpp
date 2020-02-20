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
using namespace std;
//****************************************
//Public Member Functions for Cell.cpp

//constructor
Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center, double max_radius, double init_radius, double div_site){
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->max_radius = max_radius;
    this->curr_radius = init_radius;
    this->age = 0;
    this-> CP = 0;
    this->G1 = true;
    this->G_one_length = G_one + G1*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    //cout << "G1" << G_one << " " << G_one_length << endl;
    this->G2 = false;
    this->S = false;
    this->M = false;
    this->is_mother = true;
    this->mother = nullptr;
    this->curr_bud = nullptr;
    this->is_bud = false;
    this->has_bud = false;
    this->equi_point = cell_center;
    this->curr_protein = P_0;
    this->curr_div_site = div_site;
    this->div_site_vec.push_back(curr_div_site);
    this->sector = 0;
    this->at_max_size = true;
    this->slow_grow = false;
    //this->lineage is vector
    griesemer_lineage.push_back(0);
    this->growth_rate = max_radius/((average_G1_daughter) + average_G1_daughter*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100));
    return;
}
Cell::Cell(shared_ptr<Colony> my_colony, int rank, Coord cell_center, double max_radius, double init_radius, shared_ptr<Cell> mother, double protein, double div_site, vector<shared_ptr<Cell>> lineage, vector<int>g_lineage,int sector, int bin_id){
    this->my_colony = my_colony;
    this->rank = rank;
    this->cell_center = cell_center;
    this->max_radius = max_radius;
    this->curr_radius = init_radius;
    this->age = 0;
    this-> CP = 0;
    this->G1 = true;
    this->G_one_length = G_one + G1*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100.0); 
    this->G2 = false;
    this->S = false;
    this->M = false;
    this->is_mother = false;
    this->mother = mother;
    this->curr_bud = nullptr;
    this->is_bud = true;
    this->has_bud = false;
    this->equi_point = cell_center;
    this->curr_protein = protein;
    this->curr_div_site = div_site;
    this->div_site_vec.push_back(curr_div_site);
    this->sector = sector;
    this->at_max_size = false;
    this->lineage = lineage;
    this->griesemer_lineage = g_lineage;
    this->bin_id = bin_id;
    this->slow_grow = false;
    this->growth_rate = max_radius/((average_G1_daughter) + average_G1_daughter*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100));
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
     //#pragma omp parallel for reduction(min:smallest_index)
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
     this->my_colony->get_mesh()->assign_cell_to_bin(smallest_index, this_cell);
     this->bin_id = smallest_index;
     return;
}
void Cell::slow_grow_on(){
	this->slow_grow = true;
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
void Cell::set_is_mother(){
	this->is_mother = true;	
	return;
}
void Cell::set_at_max_size(){
	this->at_max_size = true;
}
void Cell::grow_cell(){
     if(!is_mother){
		if(!at_max_size && curr_radius >= max_radius){
			//cout << "Rank " << rank << " size " << curr_radius << "set max size " << endl;
			this->at_max_size = true;
			this->CP = 0;
			//cout << "progress " << CP << endl;
		}
		else if (!at_max_size){
			if(this->far_enough_from_neighbors()){
				curr_radius = curr_radius + growth_rate*dt;
				//cout << "Rank: " << rank << "Curr radius: " << curr_radius << endl;
			}
    		}	
     }
     return;
}
/*void Cell::grow_cell(){
	curr_radius = curr_radius + growth_rate*dt;
	return;
}

void Cell::update_cell_cycle(int Ti){
	if(curr_radius >= max_radius){
		S = true;
	}
	return;
}*/
void Cell::update_cell_cycle(int Ti){
	shared_ptr<Cell> this_cell = shared_from_this();
	if(at_max_size && G1){
		//cout << 2.1/G_one_length << " " << G_one_length << endl;
		this->CP = CP + (1/G_one_length)*dt;
		//cout << "Rank" << rank<<  "Cell Prog" << CP << endl;
	}
	if((CP>=1) && G1){
		S = true;
		//cout << "Rank" << rank << "gets bud" << endl;
		//this->perform_budding(Ti);
	}
	if(is_bud){
		if(curr_radius >= size_at_div){
			//cout << "Rank" << rank << "is bud off" << endl;
			this_cell->get_mother()->enter_mitosis(Ti);
			//this_cell->perform_mitosis(Ti);
		}
	}	
	return;
}

void Cell::enter_mitosis(int Ti){
     G1 = false;
     S = false;
     G2 = false;
     M = true;
     return;
}

void Cell::perform_budding(int Ti){
    //make bud
    //cout << "Rank: " << this->get_rank() << "Bud Formed: " << Ti << "Curr Radius: " << this->curr_radius << endl;
    this-> age = age+1;
    shared_ptr<Cell> this_cell = shared_from_this();
    shared_ptr<Colony> this_colony = this->get_Colony();
    int new_rank = this_colony->get_Num_Cells();
    double init_radius = 0;
    double new_division_site;
    bool used = false;
    uniform_real_distribution<> dis(0.0,1.0);
		
    if(HAPLOID){
    	//division will happen next to last division site
	//first one assigned is randomly placed
	//curr div site is most recent
	//div site vector stores them all
	//division_site = this->curr_div_site;
	//50% chance to go above or below
	//how to implement this?
    }
    else{
    	if(is_mother){//need a mother check here
        	//mother has 50% chance opposite or adjacent
		if(this->get_Colony()->uniform_random_real_number(0.0,1.0) <= .5){
			//opposite
			//cout << "Rank: " << rank << endl;
			//cout << "current: " << curr_div_site << endl;
			new_division_site = curr_div_site + pi;
			if(new_division_site > 2*pi){
				new_division_site = new_division_site - 2*pi;
			}
			//cout << "new site: " << new_division_site << endl;
			//scan through current vector to determine if this site was used before
			for(unsigned int i=0;i<this->div_site_vec.size();i++){
				if(new_division_site == div_site_vec.at(i)){
					used = true;
					//cout << "used true" << endl;
				}
			}	
			//if not used new_division_site is unchanged
			//if used decide to move up or down
			if(used){
				//cout << "got in here " << endl;
				if(this->get_Colony()->uniform_random_real_number(0.0,1.0) <= .5){
					//move up
					do{
						used = false;
						//cout << "Stuck in while???" << endl;
						new_division_site = new_division_site - DIV_SHIFT_RADIANS;
						for(unsigned int i=0;i<this->div_site_vec.size();i++){
							if(new_division_site == div_site_vec.at(i)){
								used = true;
							}
						}	
					}while(used);	
			
				}
				else{
					//move down
					do{
						used = false;
						//cout << "stuck in while false???" << endl;
						new_division_site = new_division_site + DIV_SHIFT_RADIANS;
						for(unsigned int i=0;i<this->div_site_vec.size();i++){
							if(new_division_site == div_site_vec.at(i)){
								used = true;
							}
						}	
					}while(used);	
				}
			}	
		}
		else{
			//adjacent or same side
			//cout << "Rank: " << rank << endl;
			//cout << "current: " << curr_div_site << endl;
			new_division_site = curr_div_site;
			//cout << "new site: " << new_division_site << endl;
			
			if(this->get_Colony()->uniform_random_real_number(0.0,1.0) <= .5){
				//move 10 degrees up until new
				do{
					used = false;
					//cout << "stuck in while" << endl;
					new_division_site = new_division_site - DIV_SHIFT_RADIANS;
					for(unsigned int i=0;i<this->div_site_vec.size();i++){
						if(new_division_site == div_site_vec.at(i)){
							used = true;
						}
					}	
				}while(used);	
			}	
			else{
				//move 10 down until new
				do{
					used = false;
					//cout << "stuck in while" << endl;
					new_division_site = new_division_site + DIV_SHIFT_RADIANS;
					for(unsigned int i=0;i<this->div_site_vec.size();i++){
						if(new_division_site == div_site_vec.at(i)){
							used = true;
						}
					}	
				}while(used);	
			}

		}
        }
	else{ 
		//is daughter
		new_division_site = curr_div_site + pi;
		set_is_mother();
	}
	
    }
    
    //new_division_site = 2*pi*this->get_Colony()->uniform_random_real_number(0.0,1.0);
    //uniform_real_distribution<> div_dist(0.0,1.0);
     //new_division_site = pi*(this->get_Colony()->uniform_random_real_number(0.0,1.0));
     
    this->div_site_vec.push_back(new_division_site);
    this->curr_div_site = new_division_site;
    //cout << "rank " << this-> rank << " divsite " << division_site << endl;
    double new_max_radius = average_radius + average_radius*(this->my_colony->uniform_random_real_number(-10.0,10.0)/100);
    /*double new_center_x_D = this->cell_center.get_X()+(init_radius)*cos(curr_div_site);
    double new_center_y_D = this->cell_center.get_Y()+(init_radius)*sin(curr_div_site);
    double new_center_x_M = this->cell_center.get_X()+(init_radius)*cos(curr_div_site+pi);
    double new_center_y_M = this->cell_center.get_Y()+(init_radius)*sin(curr_div_site+pi);*/
    //double new_center_x_M = this->cell_center.get_Y();
    double new_center_x = this->cell_center.get_X()+(curr_radius+init_radius)*cos(curr_div_site);
    double new_center_y = this->cell_center.get_Y()+(curr_radius+init_radius)*sin(curr_div_site);
    Coord new_center = Coord(new_center_x,new_center_y);
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
    auto new_cell = make_shared<Cell>(this_colony, new_rank, new_center,new_max_radius, init_radius,this_cell,protein, new_division_site+pi, new_lineage, new_lineage_g, sector, this->bin_id);
    my_colony->get_mesh()->assign_cell_to_bin(this->bin_id,new_cell);
    this->add_daughter(new_cell);
    this->curr_bud = new_cell;
    this->has_bud = true;
    this->CP = 0;
    //this->curr_radius = init_radius;
    //this->cell_center = Coord(new_center_x_M,new_center_y_M); 
    G1 = false;
    S = false;
    G2 = true;
    M = false;
    /*if(this->slow_grow){
    	new_cell->slow_grow_on();
    }
    else if(this->get_Colony()->uniform_random_real_number(0.0,1.0) <= .5){
    	new_cell->slow_grow_on();
	new_cell->change_gr();
    }*/
    this_colony->update_Colony_Cell_Vec(new_cell);
    //new cell is G1 and old cell is G2
    this->equi_point = Coord(curr_radius*cos(new_division_site + pi/2),curr_radius*sin(new_division_site + pi/2));
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
void Cell::get_cell_force(){
    Coord force;
    vector<shared_ptr<Cell>> neighbor_cells;
    this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id, neighbor_cells);
    //this->my_colony->get_Cells(neighbor_cells);
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
    //this->my_colony->get_mesh()->get_cells_from_bin(this->bin_id, neighbor_cells);
    //this->my_colony->get_Cells(neighbor_cells);
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
    Coord mom_center = this->mother->get_Cell_Center();
    Coord equi_point = this->mother->get_equi_point();
    shared_ptr<Cell> this_cell = shared_from_this();
    //for(unsigned int i = 0; i< neighbor_cells.size();i++){
    	if(my_neighbor != this_cell){
		neighbor_loc = my_neighbor->get_Cell_Center();
		neighbor_radius = my_neighbor->get_radius();
		d_ij = (my_loc-neighbor_loc).length();
		v_ij = (my_loc - neighbor_loc);
		sqrt_term = sqrt((my_radius*neighbor_radius)/(my_radius+neighbor_radius));
		if(ADHESION_ON){
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
			//adh_force_reg += v_ij*SINGLE_BOND_BIND_ENERGY*RECEPTOR_SURF_DENSITY*KB*TEMPERATURE*M_PI*(my_radius + neighbor_radius)*.5*-1;
			//cout << "rep_force" << rep_force << endl;
		}
		
		
	
     }
     //cout << "Rep_force" << rep_force << endl;
     curr_force = rep_force + adh_force; //adh_force_reg + adh_force;// + bending_force;
     return curr_force;
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
    ofs << " " << this->get_sector() << " " << this->get_age() << " " << this->get_bud_status() << endl;

    return;
}
void Cell::print_cell_center(ofstream& ofs){
    ofs << cell_center.get_X() << " " << cell_center.get_Y() << " " << 0 << endl;
    return;
}
	
