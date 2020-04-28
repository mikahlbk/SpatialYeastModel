// Main.cpp

//***********************************
// Include Dependencies
#include <iostream>
#include <math.h>
#include <cstring>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <memory>
#include <random>
#include <stdio.h>
#include "colony.h"
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "mesh.h"
#include "mesh_pt.h"
//****************************************

using namespace std;

//*****************************************

int ADH_ON = 1;
double SINGLE_BOND_BIND_ENERGY = 25;
double P_0 = 50;
double r_LOGISTIC = 1;
double A_LOGISTIC = 35;
int main(int argc, char* argv[]) {
    //cout << "Starting" << endl;
    //reads in name of folder to store output for visualization
    string anim_folder = argv[1];
    for(int i = 1; i < argc; i++){
    	if(!strcmp(argv[i], "-ADH")){
		ADH_ON = stod(argv[i+1]);
	}else if(!strcmp(argv[i],"-HERTZ_ADH")){
		SINGLE_BOND_BIND_ENERGY = stod(argv[i+1]);
	}else if(!strcmp(argv[i],"-Initial_Protein")){
		P_0 = stod(argv[i+1]);
	}else if(!strcmp(argv[i],"-growth_rate")){
		r_LOGISTIC = stod(argv[i+1]);
	}else if(!strcmp(argv[i],"-competition_term")){
		A_LOGISTIC = stod(argv[i+1]);
	}



    }
    //keeps track of simulation time
    int start = clock();
    //cout << "clock" << endl;
    string init_colony = "mixed_initial.csv";
    //cout << "read in text file" << endl;
    //initialize seed for generating random numbers
    //is fed to colony constructor so that colony holds
    //the same seed and can give to other classes
    
    std::random_device seed;
    std::mt19937 gen(seed());
    //std::mt19937 gen(6);
   //make mesh for bucketing
    //cout << "make mesh" << endl;
    auto mesh_for_bins = make_shared<Mesh>();
    //leftmost point for mesh
    int start_1 = -300;
    //rightmost point for mesh
    int start_2 = 300;
    //each square unit on mesh will be this many units
    double increment = 20.0;

    int num_buckets = 2*ceil(start_2/increment);

    //cout << " make bins " << endl;
    mesh_for_bins->make_mesh_pts(start_1,start_2,num_buckets,increment);
    
    //cout << "make neighbors" << endl;
    mesh_for_bins->assign_neighbors();
    
    //make colony object
    //gen is the seed for random numbers
    auto growing_Colony = make_shared<Colony>(mesh_for_bins,gen);
    //cout << "Made Colony" << endl;
    
    //make founder cell
    //growing_Colony->make_founder_cell(init_colony);
    growing_Colony->make_founder_cell();
    //cout << "Made Founder Cell" << endl;
    //growing_Colony->match_up();
    //cout << "Match up" << endl;
    //variables for writing output files
    string format = ".txt";
    string initial = "/locations";
    int out = 1;
    ofstream myfile;
    string Number;
    string Filename;

    //not in use********************************
    //some variables for writing vtk files
	/*int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Spatial_Model_Yeast_";
	string Filename;
	ofstream ofs_anim;
	int out = 0;*/
   //*****************************************
   
   //loop for time steps
   for (int Ti = 0; Ti < NUM_STEPS; Ti++) {
  	 //write data to txt file
	//change OUTPUT_FREQ to smaller number in parameters.h
	//if want to see more timesteps 
	//cout << "In time loop" << endl;
	if(Ti%OUTPUT_FREQ == 0){
        	//open txt file for writing cell data
            	Number = to_string(out);
           	Filename = anim_folder + initial + Number + format;
            	myfile.open(Filename.c_str());
            	growing_Colony->write_data(myfile);
            	myfile.close();
            	out++;
        } 
   	
	//cout << "Time: " << Ti << endl;
        
	//assign each cell to closest bin
	//for computing forces
	if(Ti%1000 == 0){
		growing_Colony->find_bin();
		//cout << "bins" << endl;
	}
	//growth
	//cout << "grow" << endl;
	growing_Colony->grow_cells();
		
	//cell_cycle
	//cout << "cell cycle" << endl;
        growing_Colony->update_cell_cycles(Ti);

    	//budding
        //cout<< "budding" << endl;
	growing_Colony->perform_budding(Ti);

	//remove buds that are big enough
        //cout << "mitosis" << endl;
	growing_Colony->perform_mitosis(Ti);
		
	//pulling test for mother-bud adhesion
	//growing_Colony->pull_daughter();

	//spatial rearrangment
	//cout << "rearrange" << endl;
	growing_Colony->update_locations();
	//cout << "rearranged" << endl;
       
        //compute protein concentration
        //cout << "Protein Conc" << endl;
        growing_Colony->update_protein_concentration();
	//cout << "protein end" << endl;
      	//not in use********************************
	//make vtk files
	/*if(Ti%100 == 0){digits = ceil(log10(out +1));
		if(digits == 1 || digits == 0){
			Number = "0000" + to_string(out);
		}
		else if(digits == 2){
			Number = "000" + to_string(out);
	    	}
		else if(digits == 3){
			Number = "00" + to_string(out);
		}
		else if(digits == 4){
			Number = "0" + to_string(out);
		}
		
		Filename = anim_folder + initial + Number + format;
		ofs_anim.open(Filename.c_str());
		growing_Colony->print_vtk_file(ofs_anim);
		ofs_anim.close();
		out++;
	}*/
     }
     //open txt file for writing cell data
     Number = to_string(out);
     Filename = anim_folder + initial + Number + format;
     myfile.open(Filename.c_str());
     growing_Colony->write_data(myfile);
     myfile.close();
  
     int stop = clock();
     //cout << "Time: " << (stop-start) / double(CLOCKS_PER_SEC)*1000 << endl;
     //Need to add way to store data over multiple runs
    
    return 0;
}

