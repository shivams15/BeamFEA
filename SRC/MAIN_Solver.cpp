#include "Grid.h"
#include "FEMSolver.h"
#include "petsc.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
	PetscErrorCode ierr;
	ierr = PetscInitialize(&argc, &argv, (char *)0, "Initializing Program");

	if(argc < 3)
		cout<<"Either grid data file or simulation data file not provided\n";
	else{
		Grid *grid = new Grid(argv[1]);
		if(grid->setup){
			FEMSolver *fSolver = new FEMSolver(grid, argv[2]);
			if(fSolver->setup) fSolver->Solve();
		}
	}
	return 0;
}

