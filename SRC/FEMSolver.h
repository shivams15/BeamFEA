#ifndef FEMSOLVER_H
#define FEMSOLVER_H

//types of analysis
#define STEADY 0
#define TRANSIENT 1
//types of boundary conditions
#define DISPLACEMENT 0
#define TRACTION 1
//types of forces
#define POINT 0
#define DISTRIBUTED 1
#define STATIC 0
#define DYNAMIC_LINEAR 1
#define DYNAMIC_SINUSOIDAL 2
#define DYNAMIC_STEPPED 3

#include "Grid.h"
#include "petsc.h"
#include <fstream>

struct Solution{
	Vec w;	//displacements and rotations
	Vec a;	//accelerations	
	Vec v;	//velocities
	Vec M;	//bending moments
	Vec V;	//shear forces
};

class FEMSolver{
private:
	int solverType;	//type of analysis: steady or transient
	Solution *prevSol;	//solution at previous time step
	Solution *nextSol;	//solution for the next time step
	int bcType_w[2] = {-1, -1};	//boundary condition types for transverse displacement w
	int bcType_phi[2] = {-1, -1};	//boundary condition types for rotation phi
	double bc_CTX[2][2];	//data required for implementation of boundary conditions
	double load_CTX[6];	//information on applied loads
	Mat M, K, K0;	//mass and stiffness matrices
	KSP wSolver;	//solver for displacements
	KSP fSolver;
	Grid* grid;
	double E;	//Young's modulus
	double I;	//moment of inertia
	double rho;	//density of material
	double A;	//cross sectional area
	double mu;	//linear mass density
	double dt;	//time step
	double J;	//Jacobian
	int saveIter = 50;	//solution is written after every saveIter timesteps
	double finalTime;	//end time for the simulation
	double beta, gamma;	//parameters for Newmark's method of time integration
	void AssembleMatrices();
	void ComputeForceVector(Vec *f0, double t);
	void CalculateShearForce_BendingMoment();
	void EvaluateForceAtGaussPoints(double* y, double t, int j);
	void EvaluateAtGaussPoints(double** x, double* (FEMSolver::*f)(double));
	double IntegralAB(double** A, double** B, int j, int k);
	double IntegralNq(double** A, double* q, int j);
	void ConfigureKSPSolver(KSP* solver, Mat* A);
	void CreateSolution(Solution** sol);
	void CreateMatrix(Mat* A, int n);
	void SolverSetup();
	double* N(double xi);
	double* N_xx(double xi);
	double* N_xxx(double xi);
	double* N_x(double xi);
	double F(double x, double t);
	void ExportData(int iter);
	int ParseDataFile(FILE* f1);
	int SolverInitialize();
public:
	void Solve();
	bool setup = false;
	FEMSolver(Grid* grid, char* fname);
};

#endif