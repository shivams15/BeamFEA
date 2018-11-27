#include "Grid.h"
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;

Grid::Grid(char* fname){
	FILE *f1 = fopen(fname,"r");
	if(ParseDataFile(f1)) SetupGrid();
	else cout<<"Grid setup failed!\n";
	fclose(f1);
}

//Generates the elements and nodes
void Grid::SetupGrid(){
	pts = new double[2*n+1];
	eList = new Element[n];
	h = l/n;

	for(int i=0; i<2*n+1; i++){
		pts[i] = i*h/2;
	}

	for(int i=0; i<n; i++){
		eList[i].nodes[0] = 2*i;
		eList[i].nodes[1] = 2*i + 1;
		eList[i].nodes[2] = 2*i + 2;
		if(i==0){
			eList[i].boundary[0] = BOUNDARY_LEFT;
		}
		else if(i==n-1){
			eList[i].boundary[1] = BOUNDARY_RIGHT;
		}
	}
	setup = true;
	cout<<"Grid generated successfully\n";
}


//Reads user-provided parameters from the grid data file
int Grid::ParseDataFile(FILE* f1){
	char inp[255];
	int err = 0;

	if(!f1){
		cout << "Grid data file not found!\n";
		err++;
	}
	else{
		cout<<"Reading grid data file\n";
		while(!feof(f1)){
			if(fscanf(f1, "%s", inp) > 0){
				if(!strcmp(inp, "l")){
					if(fscanf(f1, "%lf", &l) < 1) err++;
					else if(l <= 0){
						cout << "Beam length must be positive!\n";
						err++;
					}
				}
				else if(!strcmp(inp, "n")){
					if(fscanf(f1, "%d", &n) < 1) err++;
					else if(n < 1){
						cout << "Number of elements must be positive!\n";
						err++;
					}
				}
				else err++;
			}
		}		
		if(err > 0) cout << "Invalid data file format!\n";
	}
	if(err > 0) return 0;
	else return 1;
}