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
		int age;
		int T_age;
		double CP;
		bool G1;
		double G_one_length;
		bool S;
		bool G2;
		bool M;
		int mother_rank;
		double growth_rate;
		bool is_mother;
		vector<shared_ptr<Cell>> daughters;
		shared_ptr<Cell> curr_bud;
		shared_ptr<Cell> mother;
		bool is_bud;
		bool has_bud;
        	bool at_max_size;
		Coord curr_force;
		double curr_protein;
        	double curr_div_site;
		vector<double> div_site_vec;
        	vector<shared_ptr<Cell>> lineage;//this holds the pointers
		vector<int> griesemer_lineage;
		int sector;
		int bin_id;
		int color;
		bool slow_grow;
		//int mother;
		Coord equi_point;
	public:
		//Constructors
		Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double max_radius, double init_radius, double div_site, int phase, double CP, int bud_status, int Mother, int my_col);
		Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double max_radius, double init_radius, double div_site);
        	Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double max_radius, double init_radius, shared_ptr<Cell> mother, double protein, double div_site, vector<shared_ptr<Cell>> lineage, vector<int> g_lineage, int sector, int bin_id, int my_col);		
		shared_ptr<Colony> get_Colony(){return my_colony;}
		int get_rank() {return rank;}
		Coord get_Cell_Center(){return cell_center;}
		double get_radius() {return curr_radius;}
		double get_max_radius() {return max_radius;}
		int get_age() {return age;}
		bool get_G1() {return G1;}
		int get_G1_length() {return G_one_length;}
		bool get_G2() {return G2;}
		bool get_S() {return S;}
		bool get_M() {return M;}
		int get_color() {return color;}
		double get_growth_rate() {return growth_rate;}
		void mother_rank_to_ptr();
		bool get_mother_status() {return is_mother;}
		int get_phase();
		int get_T_age(){return T_age;}
		double get_CP(){return CP;}
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
		void pull_daughter();
		
};

//End Cell Class
//**************************************************************
#endif


