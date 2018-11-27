#ifndef GRID_H
#define GRID_H

//types of boundaries
#define INTERNAL 0
#define BOUNDARY_LEFT 1
#define BOUNDARY_RIGHT 2

#include <fstream>

struct Element{
	int nodes[3];	//list of nodes in the element
	int boundary[2] = {INTERNAL, INTERNAL};
};

class Grid{
private:
	double l;	//length of beam
	void SetupGrid();
	int ParseDataFile(FILE* f1);
public:
	int n;	//number of elements
	double h;	//element length
	double *pts;
	Element *eList;	//list of elements
	bool setup = false;
	Grid(char* fname);
};

#endif