#include "FEMSolver.h"
#include "Grid.h"
#include "petsc.h"
#include <string.h>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

FEMSolver::FEMSolver(Grid* grid, char* fname){
	this->grid = grid;
	J = 0.5*grid->h;
	beta = 0.25; 
	gamma = 0.5;

	FILE *f1 = fopen(fname,"r");

	if(ParseDataFile(f1) && SolverInitialize())
		SolverSetup();
	else cout<<"Solver setup failed!\n";

	fclose(f1);
}

//Sets up the matrices and configures the solver
void FEMSolver::SolverSetup(){
	AssembleMatrices();
	MatAXPY(K0,1.0,M,DIFFERENT_NONZERO_PATTERN);
	MatAXPY(K0,beta*pow(dt,2),K,DIFFERENT_NONZERO_PATTERN);

	if(solverType == STEADY)
		ConfigureKSPSolver(&wSolver, &K);
	else
		ConfigureKSPSolver(&wSolver, &K0);
	
	ConfigureKSPSolver(&fSolver, &M);

	setup = true;
	cout<<"Solver setup successful\n";
}

void FEMSolver::ConfigureKSPSolver(KSP* solver, Mat* A){
	PC pc;
	KSPCreate(PETSC_COMM_SELF, solver);
	KSPSetOperators(*solver, *A, *A);
	KSPSetType(*solver, KSPGMRES);
	KSPGetPC(*solver, &pc);
	PCSetType(pc, PCILU);
	KSPSetTolerances(*solver, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetUp(*solver);
}

//Assembles the mass and stiffness matrices
void FEMSolver::AssembleMatrices(){
	cout<<"Assembling mass and stiffness matrices\n";
	double d;
	double **N0 = new double*[6];
	double **N_XX = new double*[6];
	int l1, l2;
	EvaluateAtGaussPoints(N0, &FEMSolver::N);
	EvaluateAtGaussPoints(N_XX, &FEMSolver::N_xx);

	for(int i=0; i<grid->n; i++){
		for(int j=0; j<6; j++){
			for(int k=0; k<6; k++){
				d = E*I*IntegralAB(N_XX, N_XX, j, k);
				l1 = 2*grid->eList[i].nodes[j/2] + j%2;
				l2 = 2*grid->eList[i].nodes[k/2] + k%2;
				MatSetValues(K,1,&l1,1,&l2,&d,ADD_VALUES);
				d = mu*IntegralAB(N0, N0, j, k);
				MatSetValues(M,1,&l1,1,&l2,&d,ADD_VALUES);
			}
		}
	}

	MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

	for(int i=0; i<grid->n; i++){
		if(grid->eList[i].boundary[0] == BOUNDARY_LEFT){
			l1 = 2*grid->eList[i].nodes[0];
			if(bcType_w[0] == DISPLACEMENT){
				MatZeroRows(K,1,&l1,1.0,NULL,NULL);
				MatZeroRows(M,1,&l1,1.0,NULL,NULL);
			}
			l1 = l1+1;
			if(bcType_phi[0] == DISPLACEMENT){
				MatZeroRows(K,1,&l1,1.0,NULL,NULL);
				MatZeroRows(M,1,&l1,1.0,NULL,NULL);
			}
		}
		else if(grid->eList[i].boundary[1] == BOUNDARY_RIGHT){
			l1 = 2*grid->eList[i].nodes[2];
			if(bcType_w[1] == DISPLACEMENT){
				MatZeroRows(K,1,&l1,1.0,NULL,NULL);
				MatZeroRows(M,1,&l1,1.0,NULL,NULL);
			}
			l1 = l1+1;
			if(bcType_phi[1] == DISPLACEMENT){
				MatZeroRows(K,1,&l1,1.0,NULL,NULL);
				MatZeroRows(M,1,&l1,1.0,NULL,NULL);
			}
		}
	}
}

//Evaluates the force vector at time t
void FEMSolver::ComputeForceVector(Vec *f0, double t){
	int l1;
	double d;
	double **N0 = new double*[6];
	double *f = new double[6];
	EvaluateAtGaussPoints(N0, &FEMSolver::N);
	for(int i=0; i<grid->n; i++){
		EvaluateForceAtGaussPoints(f,t,i);
		for(int j=0; j<6; j++){
			l1 = 2*grid->eList[i].nodes[j/2] + j%2;
			if(grid->eList[i].boundary[0] == BOUNDARY_LEFT){
				if(j == 0){
					VecSetValues(*f0,1,&l1,&bc_CTX[0][0],ADD_VALUES);
					if(bcType_w[0] == DISPLACEMENT)
						continue;
				}
				else if(j == 1){
					d = -bc_CTX[1][0];
					VecSetValues(*f0,1,&l1,&d,ADD_VALUES);
					if(bcType_phi[0] == DISPLACEMENT)
						continue;
				}
			}
			else if(grid->eList[i].boundary[1] == BOUNDARY_RIGHT){
				if(j == 4){
					d = -bc_CTX[0][1];
					VecSetValues(*f0,1,&l1,&d,ADD_VALUES);
					if(bcType_w[1] == DISPLACEMENT)
						continue;
				}
				else if(j == 5){
					VecSetValues(*f0,1,&l1,&bc_CTX[1][1],ADD_VALUES);
					if(bcType_phi[1] == DISPLACEMENT)
						continue;
				}
			}
			if((int)load_CTX[0] == POINT && j%2 == 0){
				d = F(J*(j/2-1) + (i+0.5)*grid->h,t);
				VecSetValues(*f0,1,&l1,&d,ADD_VALUES);
			}
			else{
				d = IntegralNq(N0,f,j);
				VecSetValues(*f0,1,&l1,&d,ADD_VALUES);
			}
		}
	}
}

/*
Runs the solver to find displacement, velocity, and acceleration fields at each time step
Newmark's method is employed for time integration
*/
void FEMSolver::Solve(){

	cout << "Initiating solver\n";
	int iter = 0;
	double t = 0.0;

	Vec f;
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &f);

	if(solverType == STEADY){
		ComputeForceVector(&f,t);
		KSPSolve(wSolver, f, nextSol->w);
		CalculateShearForce_BendingMoment();
		ExportData(0);
	}
	else{
		Vec f0;
		VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &f0);
		do{
			t += dt;
			iter++;
			VecSet(f0, 0.0);
			ComputeForceVector(&f0,t); 
			VecAXPY(prevSol->w, dt, prevSol->v);
			VecAXPY(prevSol->w, (0.5-beta)*pow(dt,2), prevSol->a);
			MatMult(K,prevSol->w,f);
			VecAYPX(f,-1.0,f0);	//this computes the RHS vector for acceleration at the next time step
			KSPSolve(wSolver, f, nextSol->a);	//solving for acclerations
			VecAXPY(prevSol->w, beta*pow(dt,2), nextSol->a);	//computing the displacements at the next time step
			VecCopy(prevSol->w, nextSol->w);
			VecAXPY(prevSol->v,(1-gamma)*dt,prevSol->a); 
			VecAXPY(prevSol->v,gamma*dt,nextSol->a);	//computing velocities at the next time step
			VecCopy(prevSol->v, nextSol->v);
			VecCopy(nextSol->a, prevSol->a);
			CalculateShearForce_BendingMoment();
			if(iter%saveIter == 0)
				ExportData(iter);
		}while(abs(t-finalTime) > 1e-10 && t < finalTime);
		VecDestroy(&f0);
	}
	cout<<"Solve complete\n";
	VecDestroy(&f);
}

void FEMSolver::CalculateShearForce_BendingMoment(){
	double d,w;
	int l1;
	Vec f;
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &f);
	double **N0 = new double*[6];
	double **N_XX = new double*[6];
	double **N_XXX = new double*[6];
	PetscScalar *u0 = new PetscScalar[4*grid->n+2];
	PetscScalar *u1 = new PetscScalar[4*grid->n+2];
	VecGetArray(nextSol->w, &u0);
	VecGetArray(nextSol->a, &u1);
	EvaluateAtGaussPoints(N0, &FEMSolver::N);
	EvaluateAtGaussPoints(N_XX, &FEMSolver::N_xx);
	EvaluateAtGaussPoints(N_XXX, &FEMSolver::N_xxx);

	for(int i=0; i<grid->n; i++){
		for(int j=0; j<6; j++){
			w = 0;
			for(int k=0; k<6; k++){
				d = mu*E*I*IntegralAB(N0, N_XX, j, k);
				w += d*u0[2*grid->eList[i].nodes[k/2] + k%2];
			}
			l1 = 2*grid->eList[i].nodes[j/2] + j%2;
			VecSetValues(f, 1, &l1, &w, ADD_VALUES);
		}
	}

	KSPSolve(fSolver,f,nextSol->M);
	VecSet(f, 0.0);

	for(int i=0; i<grid->n; i++){
		for(int j=0; j<6; j++){
			w = 0;
			for(int k=0; k<6; k++){
				d = mu*E*I*IntegralAB(N0, N_XXX, j, k);
				w += d*u0[2*grid->eList[i].nodes[k/2] + k%2];
				d = -mu*rho*I*IntegralAB(N0, N0, j, k);
				w += d*u1[2*grid->eList[i].nodes[k/2] + k%2];
			}
			l1 = 2*grid->eList[i].nodes[j/2] + j%2;
			VecSetValues(f, 1, &l1, &w, ADD_VALUES);
		}
	}
	KSPSolve(fSolver,f,nextSol->V);

	VecDestroy(&f);
}

//Evaluates the applied force at a given location at time t
double FEMSolver::F(double x, double t){
	if((int)load_CTX[0] == POINT){
		if(abs(x - load_CTX[1]) > 1e-10)
			return 0;
		if((int)load_CTX[3] == STATIC)
			return load_CTX[5];
		if((int)load_CTX[3] == DYNAMIC_LINEAR)
			return min(load_CTX[5], load_CTX[4]*t);
		if((int)load_CTX[3] == DYNAMIC_SINUSOIDAL)
			return load_CTX[5]*sin(load_CTX[4]*t);
		if(t < load_CTX[4] + 1e-10)
			return load_CTX[5];
		return 0;
	}
	else{
		if((load_CTX[1] - x) > 1e-10 || (x - load_CTX[2]) > 1e-10)
			return 0;
		if((int)load_CTX[3] == STATIC)
			return load_CTX[5];
		if((int)load_CTX[3] == DYNAMIC_LINEAR)
			return min(load_CTX[5], load_CTX[4]*t);
		if((int)load_CTX[3] == DYNAMIC_SINUSOIDAL)
			return load_CTX[5]*sin(load_CTX[4]*t);
		if(t < load_CTX[4] + 1e-10)
			return load_CTX[5];
		return 0;
	}
}

//Evaluates the shape functions at a given location
double* FEMSolver::N(double xi){
	double *N0 = new double[6];
	N0[0] = 0.75*pow(xi,5) - 0.5*pow(xi,4) - 1.25*pow(xi,3) + pow(xi,2);
	N0[1] = J*(0.25*pow(xi,5) - 0.25*pow(xi,4) - 0.25*pow(xi,3) + 0.25*pow(xi,2));
	N0[2] = pow(xi,4) -2*pow(xi,2) + 1;
	N0[3] = J*(pow(xi,5) -2*pow(xi,3) + xi);
	N0[4] = -0.75*pow(xi,5) - 0.5*pow(xi,4) + 1.25*pow(xi,3) + pow(xi,2);
	N0[5] = J*(0.25*pow(xi,5) + 0.25*pow(xi,4) - 0.25*pow(xi,3) - 0.25*pow(xi,2));
	return N0;
}

//Evaluates the first derivative shape functions at a given location
double* FEMSolver::N_x(double xi){
	double *N0_x = new double[9];
	N0_x[0] = (3.75*pow(xi,4) - 2*pow(xi,3) - 3.75*pow(xi,2) + 2*xi)/J;
	N0_x[1] = 1.25*pow(xi,4) - pow(xi,3) - 0.75*pow(xi,2) + 0.5*xi;
	N0_x[2] = (4*pow(xi,3) - 4*xi)/J;
	N0_x[3] = 5*pow(xi,4) - 6*pow(xi,2) + 1;
	N0_x[4] = (-3.75*pow(xi,4) -2*pow(xi,3) + 3.75*pow(xi,2) + 2*xi)/J;
	N0_x[5] = 1.25*pow(xi,4) + pow(xi,3) - 0.75*pow(xi,2) -0.5*xi;
	return N0_x;
}

//Evaluates the second derivative shape functions at a given location
double* FEMSolver::N_xx(double xi){
	double *N0_xx = new double[9];
	N0_xx[0] = (15*pow(xi,3) - 6*pow(xi,2) - 7.5*xi + 2)/pow(J,2);
	N0_xx[1] = (5*pow(xi,3) - 3*pow(xi,2) - 1.5*xi + 0.5)/J;
	N0_xx[2] = (12*pow(xi,2) - 4)/pow(J,2);
	N0_xx[3] = (20*pow(xi,3) - 12*xi)/J;
	N0_xx[4] = (-15*pow(xi,3) - 6*pow(xi,2) + 7.5*xi + 2)/pow(J,2);
	N0_xx[5] = (5*pow(xi,3) + 3*pow(xi,2) - 1.5*xi - 0.5)/J;
	return N0_xx;
}

//Evaluates the third derivative shape functions at a given location
double* FEMSolver::N_xxx(double xi){
	double *N0_xxx = new double[9];
	N0_xxx[0] = (45*pow(xi,2) - 12*xi - 7.5)/pow(J,3);
	N0_xxx[1] = (15*pow(xi,2) - 6*xi - 1.5)/pow(J,2);
	N0_xxx[2] = 24*xi/pow(J,3);
	N0_xxx[3] = (60*pow(xi,2) - 12)/pow(J,2);
	N0_xxx[4] = (-45*pow(xi,2) - 12*xi + 7.5)/pow(J,3);
	N0_xxx[5] = (15*pow(xi,2) + 6*xi - 1.5)/pow(J,2);
	return N0_xxx;
}

//Evaluates the forces at the Gauss quadrature points of a given element at time t
void FEMSolver::EvaluateForceAtGaussPoints(double* y, double t, int j){
	double x;
	double xi[6] = {0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521};
	for(int i=0; i<6; i++){
		x = J*xi[i] + (j+0.5)*grid->h;
		y[i] = F(x,t);
	}
}

//Evaluates the given function at the Gauss quadrature points of an element
void FEMSolver::EvaluateAtGaussPoints(double** x, double* (FEMSolver::*f)(double)){
	double xi[6] = {0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521};
	for(int i=0; i<6; i++){
		x[i] = (this->*f)(xi[i]);
	}
}

//Evaluates the integral of a product of two shape functions(or their derivatives)
double FEMSolver::IntegralAB(double** A, double** B, int j, int k){
	double w[6] = {0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704};
	double I = 0;
	for(int i=0; i<6; i++){
		I += w[i]*J*A[i][j]*B[i][k];
	}
	return I;
}

//Evaluates the integral of the product of a shape function and the force q over a given element
double FEMSolver::IntegralNq(double** A, double* q, int j){
	double w[6] = {0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704};
	double I = 0;
	for(int i=0; i<6; i++){
		I += w[i]*J*A[i][j]*q[i];
	}
	return I;
}

//Creates a new Solution structure
void FEMSolver::CreateSolution(Solution** sol){
	*sol = new Solution;
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &((*sol)->w));
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &((*sol)->v));
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &(*sol)->a);
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &(*sol)->M);
	VecCreateSeq(PETSC_COMM_SELF, 4*grid->n+2, &(*sol)->V);
	VecSet((*sol)->w,0.0);
	VecSet((*sol)->v,0.0);
	VecSet((*sol)->a,0.0);
	VecSet((*sol)->M,0.0);
	VecSet((*sol)->V,0.0);
}

//Creates a new matrix of size nxn
void FEMSolver::CreateMatrix(Mat* A, int n){
	MatCreate(PETSC_COMM_SELF, A);
	MatSetType(*A, MATAIJ);
	MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, n, n);
	MatSetFromOptions(*A);
	MatSetUp(*A);
	MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
}

//Exports  data
void FEMSolver::ExportData(int iter){
	char f1[255], f2[255], f3[255], f4[255];
	sprintf(f1, "w_%d.csv", iter);
	sprintf(f2, "phi_%d.csv", iter);
	sprintf(f3, "V_%d.csv", iter);
	sprintf(f4, "M_%d.csv", iter);

	FILE *f = fopen(f1,"w");
	PetscScalar *u0 = new PetscScalar[4*grid->n+2];
	VecGetArray(nextSol->w, &u0);

	for(int i = 0; i<2*grid->n+1; i++){
		fprintf(f, "%lf,%lf,%lf\n", grid->pts[i], 0.0, u0[2*i]);
	}
	fclose(f);

	f = fopen(f2,"w");
	for(int i = 0; i<2*grid->n+1; i++){
		fprintf(f, "%lf,%lf,%lf\n", grid->pts[i], 0.0, u0[2*i+1]);
	}
	fclose(f);

	VecGetArray(nextSol->V, &u0);
	f = fopen(f3,"w");
	for(int i = 0; i<2*grid->n+1; i++){
		fprintf(f, "%lf,%lf,%lf\n", grid->pts[i], 0.0, u0[2*i]);
	}
	fclose(f);

	VecGetArray(nextSol->M, &u0);
	f = fopen(f4,"w");
	for(int i = 0; i<2*grid->n+1; i++){
		fprintf(f, "%lf,%lf,%lf\n", grid->pts[i], 0.0, u0[2*i]);
	}
	fclose(f);
}

int FEMSolver::SolverInitialize(){
	if(solverType == TRANSIENT && dt <= 0){
		cout << "Time step must be positive!\n";
		return 0;
	}
	if(E <= 0){
		cout << "Young's modulus must be positive!\n";
		return 0;
	}
	if(rho <= 0){
		cout << "Density must be positive!\n";
		return 0;
	}
	if(A <= 0){
		cout << "Area must be positive!\n";
		return 0;
	}
	if(I <= 0){
		cout << "Moment of inertia must be positive!\n";
		return 0;
	}
	if(solverType == TRANSIENT && (finalTime <= 0 || finalTime < dt)){
		cout << "Final time must be positive and greater than dt!\n";
		return 0;
	}
	if(bcType_w[0] < 0 || bcType_w[1] < 0 || bcType_phi[0] < 0 || bcType_phi[1] < 0){
		cout << "Specify boundary conditions at all boundaries!\n";
		return 0;
	}
	if(saveIter <= 0) 
		saveIter = 50;

	if(solverType == STEADY){
		dt = 0.0;
		finalTime = 0.0;
	}

	mu = rho*A;

	CreateMatrix(&K, 4*grid->n+2);
	CreateMatrix(&K0, 4*grid->n+2);
	CreateMatrix(&M, 4*grid->n+2);
	
	CreateSolution(&prevSol);
	CreateSolution(&nextSol);

	return 1;
}

//Reads user-provided parameters from the data file
int FEMSolver::ParseDataFile(FILE* f1){
	char inp[255];
	int err = 0;

	if(!f1){
		cout << "Simulation data file not found!\n";
		err++;
	}
	else{
		cout << "Reading simulation data file\n";
		while(!feof(f1)){
			if(fscanf(f1, "%s", inp) > 0){
				if(!strcmp(inp, "analysis")){
					if(fscanf(f1, "%d", &solverType) < 1) err++;
					else if(solverType < 0 || solverType > 1){
						cout<<"Invalid analysis type!\n";
						err++;
					}
				}
				else if(!strcmp(inp, "dt")){
					if(fscanf(f1, "%lf", &dt) < 1) err++;
				}
				else if(!strcmp(inp, "E")){
					if(fscanf(f1, "%lf", &E) < 1) err++;
				}
				else if(!strcmp(inp, "rho")){
					if(fscanf(f1, "%lf", &rho) < 1) err++;
				}
				else if(!strcmp(inp, "A")){
					if(fscanf(f1, "%lf", &A) < 1) err++;
				}
				else if(!strcmp(inp, "leftbc_w")){
					if(fscanf(f1, "%d", &bcType_w[0]) < 1) err++;
					else if(bcType_w[0] > 1 || bcType_w[0] < 0){
						cout << "Invalid boundary condition type!\n";
						err++;
					}
					else if(fscanf(f1, "%lf", &bc_CTX[0][0]) < 1)
						err++;
				}
				else if(!strcmp(inp, "leftbc_phi")){
					if(fscanf(f1, "%d", &bcType_phi[0]) < 1) err++;
					else if(bcType_phi[0] > 1 || bcType_phi[0] < 0){
						cout << "Invalid boundary condition type!\n";
						err++;
					}
					else if(fscanf(f1, "%lf", &bc_CTX[1][0]) < 1)
						err++;
				}
				else if(!strcmp(inp, "rightbc_w")){
					if(fscanf(f1, "%d", &bcType_w[1]) < 1) err++;
					else if(bcType_w[1] > 1 || bcType_w[1] < 0){
						cout << "Invalid boundary condition type!\n";
						err++;
					}
					else if(fscanf(f1, "%lf", &bc_CTX[0][1]) < 1)
						err++;
				}
				else if(!strcmp(inp, "rightbc_phi")){
					if(fscanf(f1, "%d", &bcType_phi[1]) < 1) err++;
					else if(bcType_phi[1] > 1 || bcType_phi[1] < 0){
						cout << "Invalid boundary condition type!\n";
						err++;
					}
					else if(fscanf(f1, "%lf", &bc_CTX[1][1]) < 1)
						err++;
				}
				else if(!strcmp(inp, "final_time")){
					if(fscanf(f1, "%lf", &finalTime) < 1) 
						err++;
				}
				else if(!strcmp(inp, "saveIter")){
					if(fscanf(f1, "%d", &saveIter) < 1) 
						err++;
				}
				else if(!strcmp(inp, "I")){
					if(fscanf(f1, "%lf", &I) < 1) 
						err++;
				}
				else if(!strcmp(inp, "load")){
					if(fscanf(f1, "%lf", &load_CTX[0]) < 1) 
						err++;
					else if(load_CTX[0] < 0 || load_CTX[0] > 1)
						err++;
					else if(load_CTX[0] == POINT){
						if(fscanf(f1, "%lf", &load_CTX[1]) < 1) err++;
						else if(fscanf(f1, "%lf", &load_CTX[3]) < 1)
							err++;
						else if(load_CTX[3] < 0 || load_CTX[3] > 3)
							err++;
						else if((int)load_CTX[3] == STATIC){
							if(fscanf(f1, "%lf", &load_CTX[5]) < 1)
								err++;
						}
						else{
							if(fscanf(f1, "%lf", &load_CTX[4]) < 1)
								err++;
							else if(fscanf(f1, "%lf", &load_CTX[5]) < 1)
								err++;
						}
					}
					else{
						if(fscanf(f1, "%lf", &load_CTX[1]) < 1) err++;
						else if(fscanf(f1, "%lf", &load_CTX[2]) < 1) err++;
						else if(load_CTX[1] >= load_CTX[2])
							err++;
						else if(fscanf(f1, "%lf", &load_CTX[3]) < 1)
							err++;
						else if(load_CTX[3] < 0 || load_CTX[3] > 3)
							err++;
						else if((int)load_CTX[3] == STATIC){
							if(fscanf(f1, "%lf", &load_CTX[5]) < 1)
								err++;
						}
						else{
							if(fscanf(f1, "%lf", &load_CTX[4]) < 1)
								err++;
							else if(fscanf(f1, "%lf", &load_CTX[5]) < 1)
								err++;
						}
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