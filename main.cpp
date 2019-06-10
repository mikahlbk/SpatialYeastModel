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

//****************************************

using namespace std;

//*****************************************

int main(int argc, char* argv[]) {
    
    string anim_folder = argv[1];

    int start = clock();
	//cout << "clock" << endl;
	//create initial colony with specified number of cells
	int init_cells = 1;
	cout << init_cells << endl;
	auto growing_Colony = make_shared<Colony>();
	cout << "Made colony" << endl;
	growing_Colony->make_cells(init_cells);
	cout << "Made new colony with initial number of cells" << endl;
	
	//some variables for writing data
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Spatial_Model_Yeast_";
	string Filename;
	ofstream ofs_anim;
	int out = 0;

	//variable for main loop
	int numSteps = 500;

	//loop for time steps
	for (int Ti = 0; Ti*dt < numSteps; Ti++) {
		
		//growth
		growing_Colony->grow_cells();
		
		//cell_cycle
	    cout << "growth" << endl;
        growing_Colony->update_cell_cycles();

		//budding
        cout<< "budding" << endl;
		growing_Colony->perform_budding();

		//remove buds that are big enough
        cout << "mitosis" << endl;
		growing_Colony->perform_mitosis();

		//spatial rearrangment
		cout << "rearrange" << endl;
		growing_Colony->update_Locations();
	
		//make vtk files
		if(Ti%100 == 0){digits = ceil(log10(out +1));
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
	}
    }
	int stop = clock();
	cout << "Time: " << (stop-start) / double(CLOCKS_PER_SEC)*1000 << endl;
	//Need some kind of component to save
	//data over each run????
	
}

