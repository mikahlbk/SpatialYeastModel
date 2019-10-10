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

#include "colony.h"
#include "parameters.h"
#include "coord.h"
#include "cell.h"
#include "mesh.h"
#include "mesh_pt.h"
//****************************************

using namespace std;

//*****************************************

int main(int argc, char* argv[]) {
    
    //reads in name of folder to store output
    string anim_folder = argv[1];

    //keeps track of simulation time
    int start = clock();
    //cout << "clock" << endl;
    srand(1); 
    //make mesh
    //cout << " make mesh " << endl;
    auto mesh_for_bins = make_shared<Mesh>();
    int start_1 = -100;
    int end = 100;
    double increment = 10.0;
    //cout << " make bins " << endl;
    mesh_for_bins->make_mesh_pts(start_1,end,increment);
    //cout << "make neighbors" << endl;
    mesh_for_bins->assign_neighbors();
    //make colony object
    auto growing_Colony = make_shared<Colony>(mesh_for_bins);
    //cout << "Made Colony" << endl;
    growing_Colony->make_founder_cell();
    //cout << "Made founder cell" << endl;
	
    //variables for writing output files
    string format = ".txt";
    string initial = "/locations";
    int out = 1;
    ofstream myfile;
    string Number;
    string Filename;

    //not in use
    //some variables for writing vtk files
	/*int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Spatial_Model_Yeast_";
	string Filename;
	ofstream ofs_anim;
	int out = 0;*/

	//variable for main loop
	int numSteps = 5000;

	//loop for time steps
	for (int Ti = 0; Ti*dt < numSteps; Ti++) {
		
		//here I will call the make sector function
		
		cout << "Time: " << Ti << endl;
        	growing_Colony->find_bin();
		//growth
		cout << "grow" << endl;
		growing_Colony->grow_cells();
		
		//cell_cycle
	    	cout << "cell cycle" << endl;
        	growing_Colony->update_cell_cycles();

    		//budding
        	cout<< "budding" << endl;
		growing_Colony->perform_budding(Ti);

		//remove buds that are big enough
        	cout << "mitosis" << endl;
		growing_Colony->perform_mitosis(Ti);

		//spatial rearrangment
		cout << "rearrange" << endl;
		growing_Colony->update_locations();
	    	//cout << "rearranged" << endl;
       
        	//compute protein concentration
        	cout << "Protein Conc" << endl;
        	growing_Colony->update_protein_concentration();

        	//write data to txt file
		//change to smaller number if need to see stuff more frequently 
	    	if(Ti%100 == 0){
            		//open txt file for writing cell data
            		Number = to_string(out);
           		Filename = anim_folder + initial + Number + format;
            		myfile.open(Filename.c_str());
            		growing_Colony->write_data(myfile);
            		myfile.close();
            		out++;
        	}
        
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
     int stop = clock();
     cout << "Time: " << (stop-start) / double(CLOCKS_PER_SEC)*1000 << endl;
     //Need to add way to store data over multiple runs
    
    return 0;
}

