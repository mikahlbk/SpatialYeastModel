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
		int type;
		bool G1;
		bool S;
		bool G2;
		bool M;
		vector<shared_ptr<Cell>> daughters;
		shared_ptr<Cell> curr_bud;
		shared_ptr<Cell> mother;
		bool is_bud;
		bool has_bud;
        Coord curr_force;
		double curr_protein;
        double div_site1;
		double div_site2;
        vector<shared_ptr<Cell>> lineage;
	
	public:
		//Constructor
		Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double max_radius, double init_radius);
        Cell(shared_ptr<Colony> colony, int rank, Coord cell_center, double max_radius, double init_radius, shared_ptr<Cell> mother, double protein, double div_site, vector<shared_ptr<Cell>> lineage);		
		Coord get_Cell_Center(){return cell_center;}
		shared_ptr<Cell> get_mother(){return mother;}
		shared_ptr<Cell> get_bud() {return curr_bud;}
		void get_daughters(vector<shared_ptr<Cell>>& daughter_cells);
		bool get_G1() {return G1;}
		bool get_G2() {return G2;}
		bool get_S() {return S;}
		bool get_M() {return M;}
        void get_lineage_vec(vector<shared_ptr<Cell>>& lineage_cells);
		double get_protein_conc(){return curr_protein;}
        double get_radius() {return curr_radius;}
		double get_max_radius() {return max_radius;}
		shared_ptr<Colony> get_Colony() {return my_colony;}
		int get_rank() {return rank;}
        void add_daughter(shared_ptr<Cell> daughter);
        int get_generation();
        double calc_closest_center();
        void set_bud(shared_ptr<Cell> bud);
		void set_mother(shared_ptr<Cell> mother);
        void reset_is_bud();
        void reset_has_bud();
        void grow_cell();
		void update_cell_cycle();
		void enter_mitosis();
		void perform_budding();
		void perform_mitosis();
		void calc_forces_jonsson();
        void calc_forces_chou();
		void calc_forces_exponential();
        void lennard_jones_potential();
        void update_location(); 
        void compute_protein_concentration();
        void print_txt_file_format(ofstream& ofs);
        void print_cell_center(ofstream& ofs);	
};

//End Cell Class
//**************************************************************
#endif


