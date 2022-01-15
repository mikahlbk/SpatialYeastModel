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
//Nutrient Rich (0) vs. Nutrient Limited (1)
int Nutrient_On = 0;
//Rate of substrate consumption by cells
double K_spring = .05;
int main(int argc, char* argv[]) {
    //cout << "Program Starting" << endl;
    //reads in name of folder to store output for visualization
    string anim_folder = argv[1];
    for(int i = 1; i < argc; i++){
    	if(!strcmp(argv[i], "-Nutrient_limitation")){
		    Nutrient_On = stod(argv[i+1]);
	    }else if(!strcmp(argv[i],"-Adhesion_strength")){
		    K_spring = stod(argv[i+1]);
	    }	
    }
    cout << "Nutrient Condition = " << Nutrient_On << "Adhesion Strength = " << K_spring <<  endl;
    //keeps track of simulation time
    int start = clock(); 
    //****Make mesh for bins*****
    auto Main_mesh = make_shared<Mesh>();
    //cout << "Mesh Created" << endl;
    //Topmost point for mesh
    int start_1 = -400;
    //Rightmost point for mesh
    int start_2 = 400;
    //Each mesh poiint will be this many units
    double increment = 2*AVERAGE_MAX_RADIUS;//~1 cell per bin
    //Number of buckets across a row of mesh
    unsigned int num_buckets = 2*ceil(start_2/increment);
    //****Make mesh point****
    Main_mesh->make_mesh_pts(start_1,start_2,num_buckets,increment);
    //cout << "Mesh Points Created" << endl;
    //****Make Colony****
    auto growing_Colony = make_shared<Colony>(Main_mesh);
    //cout << "Colony Created" << endl;
    //****Make Founder Cell****
    growing_Colony->make_founder_cell(growing_Colony,0);
    //cout << "Founder Cell Created" << endl;
    //****Variables for writing output files*****
    string format          = ".txt";
    string locations_init  = "/locations";
    string cell_cycle_init = "/CellCycle";
    int out               = 1;
    ofstream myfile1;
    ofstream myfile2;
    string Number;
    string Filename1;
    string Filename2;

   //****Loop for time steps****
   for (int Ti = 0; Ti < NUM_STEPS; Ti++) {
   //cout << "Top of Time Loop" << endl;
   //****Write data to txt file
   //Change OUTPUT_FREQ to smaller number in parameters.h
   //if you want to see more timesteps 
   //cout << "In Time Loop" << endl;
   if(Ti%OUTPUT_FREQ == 0){
        //open txt file for writing cell data
        Number = to_string(out);
        Filename1 = anim_folder + locations_init + Number + format;
	    Filename2 = anim_folder + cell_cycle_init + Number + format;
        myfile1.open(Filename1.c_str());
	    myfile2.open(Filename2.c_str());
        growing_Colony->write_data(myfile1,myfile2);
        myfile1.close();
	    myfile2.close();
        out++;
    } 
   //cout << "Time: " << Ti << endl;
        
	if(Ti%BIN_UPDATE_INCREMENT == 0){
		growing_Colony->update_cell_bin_ids();
		//cout << "All Cell Bins Assigned" << endl;
		//Assigns each cell to closest bin
	}
	if(Nutrient_On){
		growing_Colony->update_cell_phase_lengths();
		//cout << "Nutrient Dependent Cell Parameters Updated" << endl;
		//Cell phase lengths, target sizes (max_size and size of bud at separation), and growth rates
		//are updated according to nutrient concentration in each bin
	}
	growing_Colony->update_cell_radii();
	//cout << "Cell Sizes Updated" << endl;
	//Updates (i.e. grows) the radius of buds
	//and unbudded daughters that are
	//still growing in G1
		
    growing_Colony->perform_division();
	//cout << "Buds Separated" << endl;
	//Separates buds who have reached their assigned
	//size of separation from their mothers
	growing_Colony->update_cell_cycle_phases();
	//cout << "Cell Cycles Updated" << endl;
	//Updates the current cell cycle phase of cells
	//growing_Colony->add_buds(Ti);
	//cout << "New Buds Added" << endl;
	//Adds a bud to mother and new daughter cells
	//that have finished G1

	//growing_Colony->update_cell_protein_concentrations();
	//cout << "Protein Concentrations Updated"<< endl;
	
	growing_Colony->update_cell_locations();
	//cout << "Cell Locations Updated" << endl; 
	
	if(Nutrient_On){
		growing_Colony->update_mesh_nutrient_concentration();
		//cout << "Mesh Nutrient Concentrations Updated" << endl;
		//For each mesh point in the mesh this function
		//updates the nutrient concentration based on the
		//substrate uptake rate of cells and the number
		//of cells in that bin.
	}

	//cout << "End of Time Loop" << endl;
     }
     //open txt file for writing cell data
	Number = to_string(out);
        Filename1 = anim_folder + locations_init + Number + format;
        Filename2 = anim_folder + cell_cycle_init + Number + format;
        myfile1.open(Filename1.c_str());
	myfile2.open(Filename2.c_str());
        growing_Colony->write_data(myfile1,myfile2);
        myfile1.close();
	myfile2.close();
        out++;     
 
      	int stop = clock();
     	cout << "Time: " << (stop-start) / double(CLOCKS_PER_SEC)*1000 << endl;
     	//Need to add way to store data over multiple runs
    
    return 0;
}

