/**************************************************************************
The main class of the walk. Some importances are listed:

constructor Walk --- initialize everything, paremeters are (total length, initial file name, pivot each step, final file name)
where "0" for initial file name means a 45 degree straight line

GoOneStep(n) --- go 1 step in MC chain with n trial, return -1 if all failed, return the number of total trial in doing pivoting if successful.
Note that all pivot will be carried out at the end of this.

run(n) --- go n steps in MC chain, which has 10 trial by default. Will call Record() once at the end.

Record() --- Write down the SAW in file "FinalWalk" by default. Or other data in future.

class proposal ---a proposal contains a pivot location, an operation, the inverse of that opeation.
Calling by using the constructor "proposal(this)" in the class Walk will initialize everything.
**************************************************************************/

#pragma once
#include "global.h"
#include "point.h"

class Walk {
	// stores a single walk
private:
	Sphere* steps;
	const int nsteps; // number of steps in the walk  
	int generation;
	int n_inner;
	int max_npivot;
	int SAWaccept;
	int no_saw;
	double old_energy;
	int npivot; // number of accepted pivots not yet applied to walk
	int* ptime; // array of pivot times 
	OPERATION_NAME* igroup; // array of group elements  
	GPoint<double>* shift; // array of shifts

	const char *init_walk_fname, *final_walk_fname, *data_fname;
	int nsimplify;

/**********************************************************

***********************************************************/
	class Proposal
	{
	public:
		int pivot_loc;
		OPERATION_NAME op;
		OPERATION_NAME invop;

		Walk* pw;
		Proposal(Walk* pw) :pivot_loc(int(pw->GetLength()*RNG_NAME())), op("rand"), invop(op.inv()){}
	};
public:
	// Function to go n step 
	void run(int outer_steps, int discard);

	int GetAutocorrelation(int n, unsigned long int intersteps = 10);

	// Central function of initialization
	Walk(int tnsteps = 1000, const char* tinit_walk_fname = "0", int tn_inner = 1000, const char* tfinal_walk_fname = "FinalWalk", int tnsimplify = 0, int tno_saw = 0, int tmax_npivot = 0);
	
	// Function to go one step in markov chain
	void GoOneStep(int stepnum, bool isrecord = true);

	// Add the pivot to the array
	void add_pivot(int pivot_loc, OPERATION_NAME* isym, GPoint<double> trans);

	// Carry out the pivots, so npivot -> 0
	void simplify(); 

	// The core, where Kennedy check is used.
	int pivot_strictly_saw(Proposal* prop);

	// A parameter to judge thermalization. Here is the fraction of 'number of steps that the angle larger than 120 degree' over 'the length
	double turn_frac(); // fraction of sites with a turn in the walk

	// Find the index of the segment of i (itime).
	inline int find_segment(int itime, int npivot, int* ptime);

	// delete every pointer
	void deallocate();

	// Initializes walk to be straight line
	void line_initialize(int direction);

	// Clean the pivot cache
	void clean_pivot();

	// A example code to calculate the a energy, lower is better. now is a empty function. 
	double GetEnergyForAll();

	double GetEnergyForOne(int ipivot, OPERATION_NAME* op);

	// Function to judge whether accept the pivots, now is empty
	bool AcceptOrNot(double newE, double oldE);

	// Called once at the end of run().
	void Record();
	
	void Writedown();

	bool Sokal(int pivot_loc);

	bool CheckAll();

	// Get the current i step
	Sphere GetStepi(int i);
	
	double GetRg2();

	int GetLength(){ return(nsteps); }
	int GetPivotNum(){ return(npivot); }
	void print(FILE *fptr);
	void scan(FILE *fptr);
// 	void plot(FILE *fptr);
// 	void plot_projected(FILE *fptr, float rho);
/*	Walk& operator=(Walk& w);*/
};