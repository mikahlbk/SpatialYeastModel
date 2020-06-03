// cell.h

//*********************************************************
// Include Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_

//*********************************************************
// forward declarations
class Colony;

//*********************************************************
// include dependencies
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <memory>
#include "parameters.h"
#include "coord.h"
#include "externs.h"
//***********************************************************
// Cell Class Declaration

class Cell: public enable_shared_from_this<Cell>{
	private:
		shared_ptr<Colony> my_colony;
		int rank;
		Coord cell_center;
		double curr_radius;
		double max_radius;
		bool at_max_size;
		Coord curr_force;
		int bin_id;
		int age;
		int T_age;
		double my_G1_length;
		double my_G2_length;
		bool G1;
		bool S;
		bool G2;
		bool M;
		double growth_rate;
		double cell_cycle_increment;
		double theoretical_cci;
		double CP;
		bool is_mother;
		bool has_bud;
		shared_ptr<Cell> curr_bud;
		vector<shared_ptr<Cell>> daughters;
		vector<double> div_site_vec;
		double curr_div_site;
		bool is_bud;
		shared_ptr<Cell> mother;
		int mother_rank;
		vector<int> lineage;
		vector<int> griesemer_lineage;
		int sector;
		double curr_protein;
        	int color;
		Coord equi_point;
	public:
		//Constructor for single founder
		Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double init_radius, double div_site);
        	//Constructor for new daughter after division
        	Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double init_radius, shared_ptr<Cell> mmother,int mother_rank,double div_site, vector<int> lineage, vector<int> g_lineage, int sector, int my_col,double g2_from_mother);	
		/*Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double max_radius, double init_radius, double div_site, int bud_status, int phase, double CP, int Mother, int my_col);*/	
		//***Getters***	
		shared_ptr<Colony> get_colony(){return my_colony;}
		int get_rank(){return rank;}
		Coord get_cell_center(){return cell_center;}
 		double get_curr_radius(){return curr_radius;}
                double get_max_radius(){return max_radius;}
		Coord get_curr_force(){return curr_force;}
		int get_bin_id(){return bin_id;}
		int get_age(){return age;}
		int get_T_age(){return T_age;}
		double get_G1_length(){return my_G1_length;}
		double get_G2_length(){return my_G2_length;}
		bool is_G1(){return G1;}
		bool is_G2(){return G2;}
		bool is_S(){return S;}
		bool is_M(){return M;}
		double get_growth_rate(){return growth_rate;}
		double get_cell_cycle_increment(){return cell_cycle_increment;}
		double get_CP(){return CP;}
		bool mother_status(){return is_mother;}
		bool currently_has_bud(){return has_bud;}
		shared_ptr<Cell> get_curr_bud(){return curr_bud;}
		void get_daughters_vec(vector<shared_ptr<Cell>>& curr_daughters);
		void get_div_site_vec(vector<double>& previous_div_sites);
		double get_curr_div_site(){return curr_div_site;}
		bool bud_status(){return is_bud;}
		shared_ptr<Cell>get_mother(){return mother;}
		void set_mother(shared_ptr<Cell> mother);
		int get_mother_rank(){return mother_rank;}
		void get_lineage_vec(vector<int>& curr_lineage_vec);
		void update_lineage_vec(int mother_rank);
		void get_griesemer_lineage_vec(vector<int>& curr_griesemer_lineage);
		int get_sector(){return sector;}
		bool grown_to_full_size(){return at_max_size;}
		double get_curr_protein(){return curr_protein;}
		int get_color(){return color;}
		int get_phase();	
		//***functions used when starting with > 1 cell***
		void mother_bud_check();
		void get_bud_status_mom(shared_ptr<Cell> mother);
		void mother_rank_to_ptr();
		void change_mother_vars(shared_ptr<Cell> mother_cell,shared_ptr<Cell> bud);
		void return_bud_status();
		//***************************************************

		//functions used to put cell in correct bin
		void find_bin();
		//functions used to adjust growth rate of cells based
		//on nutrient concentration of their current bin
		void update_growth_rate();	 		
   		double get_nutrient_conc(int bin_id);
		void grow_cell();
		void update_cell_cycle();
		void enter_mitosis();
		void perform_budding(int Ti);
		void perform_mitosis(int Ti);
		void set_has_bud_to_false();
		void set_is_bud_to_false();
		void compute_protein_concentration();
		void set_protein_conc(double protein);
		void get_cell_force();
 		Coord calc_forces_Hertz(shared_ptr<Cell> my_neighbor);
		void update_location();
		void print_txt_file_format(ofstream& ofs); 
		/*void mother_rank_to_ptr();
		int get_phase();
		void get_bud_status_mom(shared_ptr<Cell> mother);
		void set_is_mother();
		void set_at_max_size();
		void get_daughters(vector<shared_ptr<Cell>>& daughter_cells);
	        void add_daughter(shared_ptr<Cell> daughter);
        	shared_ptr<Cell> get_bud() {return curr_bud;}
		void set_bud(shared_ptr<Cell> bud);
		shared_ptr<Cell> get_mother(){return mother;}
		void set_mother(shared_ptr<Cell> mother);
        	void reset_is_bud();
        	void reset_has_bud();
        	int get_bud_status();	
		void update_growth_rate();	
		double get_nutrient_conc(int bin_id);
		void mother_bud_check();
		void set_has_bud();
		void change_mother_vars(shared_ptr<Cell> mother_cell,shared_ptr<Cell> bud_cell);
		//cell cycle
		double compute_distance(shared_ptr<Cell> neighbor_cell);
		bool far_enough_from_neighbors(); 
		void grow_cell();
		void update_cell_cycle(int Ti);
		void perform_budding(int Ti);
		void enter_mitosis(int Ti);
		void perform_mitosis(int Ti);
		void slow_grow_on();
		void change_gr();
		//force calculations
		void get_cell_force();
		Coord calc_forces_Hertz(shared_ptr<Cell> my_neighbor);
		void calc_forces_jonsson();
        	void calc_forces_chou();
		void calc_forces_exponential();
        	void lennard_jones_potential();
        	void update_location();
		Coord get_curr_force() {return curr_force;}	
		//angle of bud
		Coord get_equi_point(){return equi_point;}
		double compute_angle();

		//prion functions
		double get_protein_conc(){return curr_protein;}
		void compute_protein_conc_DNPM();
        	void compute_protein_concentration();
        	
		//lineage output
		void update_lineage_vec(shared_ptr<Cell> me); 
		void get_lineage_vec(vector<shared_ptr<Cell>>& lineage_cells);
		void return_lineage_vec(vector<int>& new_vec);
		void get_lineage_g_vec(vector<int>& lineage_g_cells);
		int get_sector(){return sector;}
      		
		//binning functions
		double calc_closest_center();
        	void find_bin();
		int get_bin_id() {return bin_id;}
		
		//visualization
		void print_txt_file_format(ofstream& ofs);
        	void print_cell_center(ofstream& ofs);	
		void pull_daughter();*/
		
};

//End Cell Class
//**************************************************************
#endif


