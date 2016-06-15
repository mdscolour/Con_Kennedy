#pragma once
#include "global.h"
#include "point.h"

class Walk {
	// stores a single walk
private:
	MODELNAME* steps;
	const int nsteps; // number of steps in the walk  

	// number of attempted pivots the walk has been through since it was a line 
	// in multiples of INNER_LOOP
	// So if INNER_LOOP=1,000,000, this is number in millions
	int generation=0;
	int n_inner;
	int max_npivot;
	int no_saw;
	double old_energy;
	int npivot; // number of accepted pivots not yet applied to walk
	int* ptime; // array of pivot times 
	OpMatrix* igroup; // array of group elements  
	MODELNAME* shift; // array of shifts

	char *init_walk_fname, *final_walk_fname, *data_fname;
	int nsimplify;

	class Proposal
	{
	public:
		int pivot_loc;
		OpMatrix op;
		OpMatrix invop;

		Walk* pw;
		//Prososal(){ pw = (Walk*)((char*)this - offsetof(Walk, (*this))); pivot_loc = int(floor(pw->get_nsteps()*RNG_NAME())); op = OpMatrix(floor(NUM_SYM * RNG_NAME() + 1)); invop = op.inv(); }
		Proposal(Walk* pw) :pivot_loc(int(floor(pw->get_nsteps()*RNG_NAME()))), op(int(floor(NUM_SYM * RNG_NAME() + 1))), invop(op.inv()){}
		void Rand(){ pivot_loc = int(floor(pw->get_nsteps()*RNG_NAME())); op = OpMatrix(int(floor(NUM_SYM * RNG_NAME() + 1))); invop = op.inv(); }
	};
public:
	Walk(int tnsteps = 1000, char* tinit_walk_fname = "0", int tnsimplify = 3, int tn_inner = 1000, int tno_saw = 0, int tmax_npivot = 0, int tprint_freq = 0);
	void run();
	
	void deallocate();
	void line_initialize(int direction);
	void initialize();

	double GetEnergy();
	void Record();
	int get_nsteps(){ return(nsteps); }
	int get_npivot(){ return(npivot); }

	MODELNAME step_rval(int i);
	void print(FILE *fptr);
	void scan(FILE *fptr);
// 	void plot(FILE *fptr);
// 	void plot_projected(FILE *fptr, float rho);
	
/*	Walk& operator=(Walk& w);*/

	void add_pivot(int pivot_loc, OpMatrix* isym, MODELNAME trans);
	void simplify(); //carry out the pivots, so npivot -> 0
	int pivot_strictly_saw(Proposal* prop);

	double turn_frac(); // fraction of sites with a turn in the walk
	inline int find_segment(int itime, int npivot, int* ptime);

};