#pragma once
#include "global.h"
#include "point.h"

class Walk {
	// stores a single walk
private:
	MODEL_NAME* steps;
	const int nsteps; // number of steps in the walk  

	// number of attempted pivots the walk has been through since it was a line 
	// in multiples of INNER_LOOP
	// So if INNER_LOOP=1,000,000, this is number in millions
	int generation;
	int n_inner;
	int max_npivot;
	int no_saw;
	double old_energy;
	int npivot; // number of accepted pivots not yet applied to walk
	int* ptime; // array of pivot times 
	OPERATION_NAME* igroup; // array of group elements  
	MODEL_NAME* shift; // array of shifts

	char *init_walk_fname, *final_walk_fname, *data_fname;
	int nsimplify;

	class Proposal
	{
	public:
		int pivot_loc;
		OPERATION_NAME op;
		OPERATION_NAME invop;

		Walk* pw;
		Proposal(Walk* pw) :pivot_loc(int(floor(pw->GetLength()*RNG_NAME()))), op(int(floor(NUM_SYM * RNG_NAME() + 1))), invop(op.inv()){}
		void Rand(){ pivot_loc = int(floor(pw->GetLength()*RNG_NAME())); op = OPERATION_NAME(int(floor(NUM_SYM * RNG_NAME() + 1))); invop = op.inv(); }
	};
public:

	// Central function of initialization
	Walk(int tnsteps = 1000, char* tinit_walk_fname = "0", int tnsimplify = 3, char* tfinal_walk_fname = "FinalWalk", int tn_inner = 0, int tno_saw = 0, int tmax_npivot = 0);
	
	// Function to go one step in markov chain
	void run();

	void add_pivot(int pivot_loc, OPERATION_NAME* isym, MODEL_NAME trans);

	// Carry out the pivots, so npivot -> 0
	void simplify(); 
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
	double GetEnergy();

	// Function to judge whether accept the pivots, now is empty
	bool AcceptOrNot(double newE, double oldE);

	// A step means at the end this function will be called once
	void Record();

	// Get the current i step
	MODEL_NAME GetStepi(int i);

	int GetLength(){ return(nsteps); }
	int GetPivotNum(){ return(npivot); }
	void print(FILE *fptr);
	void scan(FILE *fptr);
// 	void plot(FILE *fptr);
// 	void plot_projected(FILE *fptr, float rho);
/*	Walk& operator=(Walk& w);*/
};