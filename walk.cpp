#include "walk.h"

inline int Walk::find_segment(int itime, int npivot, int* ptime)
// Finds the segment number for a pivot time. 
// segment number iseg corresponds to times [ptime[iseg],ptime[iseg+1])
// Note that the left endpt, ptime[iseg], is included in segment iseg.
// Presently, routine is a slow search. This could be done better, but
// I am not sure it is worth it. 
{
	int iseg, isegl, isegu;

	if (itime >= ptime[npivot]) return(npivot);

	isegl = 0; isegu = npivot; iseg = 0;

	while (isegu > isegl + 1)
	{
		iseg = (isegl + isegu) / 2;
		if (itime < ptime[iseg]) isegu = iseg;
		else isegl = iseg;
	}
	return(isegl);
} // find_segment()

// void Walk::setup_symmetry()
// // this routine computes group_product[][] and group_inverse[] entries
// // by trial and error 
// {
// 	int i, j, k;
// 	Sphere origin, p1, p2, q1, q2, temp;
// 	Sphere p3, q3;
// 
// 	origin.zero();
// 
// 	for (i = 0; i <= NUM_SYM; i++)
// 	{
// 		cout << "product tensor doing " << i << endl;
// 		for (j = 0; j <= NUM_SYM; j++)  ///////////////////
// 		{
// 			
// 			group_product[i][j] = -1;
// 			/*p1.assign(0, 0, 1);
// 			temp = p1; p1.euclidean_op(&temp, &origin, j);
// 			temp = p1; p1.euclidean_op(&temp, &origin, i);
// 			p2.assign(0, 1, 0);
// 			temp = p2; p2.euclidean_op(&temp, &origin, j);
// 			temp = p2; p2.euclidean_op(&temp, &origin, i);
// 			p3.assign(1, 0, 0);
// 			temp = p3; p3.euclidean_op(&temp, &origin, j);
// 			temp = p3; p3.euclidean_op(&temp, &origin, i);
// 			for (k = 0; k <= NUM_SYM; k++)
// 			{
// 			q1.assign(0, 0, 1);
// 			temp = q1; q1.euclidean_op(&temp, &origin, k);
// 			q2.assign(0, 1, 0);
// 			temp = q2; q2.euclidean_op(&temp, &origin, k);
// 			q3.assign(1, 0, 0);
// 			temp = q3; q3.euclidean_op(&temp, &origin, k);
// 			if (p1 == q1 && p2 == q2 && p3 == q3) group_product[i][j] = k;
// 			}*/
// 			group_product[i][j] = RMatrix(RMatrix(i).dot(RMatrix(j))).GetIndex();
// 			if (group_product[i][j] == -1)
// 			{
// 				printf("ERROR in setup_symmetry() %ld %ld\n", i, j);
// 				system("pause");
// 				exit(1);
// 			}
// 		}
// 	}
// 	cout << "product tensor done" << endl;
// 	for (i = 0; i <= NUM_SYM; i++)
// 	{
// 		group_inverse[i] = -1;
// 		/*for (k = 0; k <= NUM_SYM; k++)
// 		{
// 			p1.assign(0, 0, 1);
// 			temp = p1; p1.euclidean_op(&temp, &origin, i);
// 			temp = p1; p1.euclidean_op(&temp, &origin, k);
// 			q1.assign(0, 0, 1);
// 			p2.assign(0, 1, 0);
// 			temp = p2; p2.euclidean_op(&temp, &origin, i);
// 			temp = p2; p2.euclidean_op(&temp, &origin, k);
// 			q2.assign(0, 1, 0);
// 			p3.assign(1, 0, 0);
// 			temp = p3; p3.euclidean_op(&temp, &origin, i);
// 			temp = p3; p3.euclidean_op(&temp, &origin, k);
// 			q3.assign(1, 0, 0);
// 			if (p1 == q1 && p2 == q2 && p3 == q3) group_inverse[i] = k;
// 		}*/
// 		group_inverse[i] = RMatrix(i).GetInverseIndex();
// 		if (group_inverse[i] == -1)
// 		{
// 			printf("ERROR in setup_symmetry()_inverse\n");
// 			printf("%ld\n", RMatrix(i).GetInverseIndex());
// 			system("pause");
// 			exit(1);
// 		}
// 	}
// 	cout << "inverse tensor done" << endl;
// 	return;
// 
// } // end setup_symmetry()

MODELNAME Walk::step_rval(int i)
// This routine computes the ith point on the walk from the complicated 
// data structure. Often we duplicate this code for speed rather than call this
{
	MODELNAME temp;
	int iseg;
	iseg = find_segment(i, npivot, ptime);
	temp.euclidean_op(steps + i, shift + iseg, &igroup[iseg]);
	return(temp);
} // end step_rval


void Walk::initialize()
// NB this does not create a walk. It only initializes things used by the 
// implicit pivot data structure
{
	npivot = 0;
	igroup[0] = 0;
	shift[0].zero();
	// segment number iseg corresponds to times [ptime[iseg],ptime[iseg+1])
	// So when npivot=0 we should have : 
	ptime[0] = 0;
	ptime[1] = nsteps + 1;
}

void Walk::line_initialize(int direction)
{ // initializes walk to be straight line
	// 1 is horizontal 
	// 2 is 45 degs
	// 3 is vertical
	int i;
	double i1, i2, i3;
	initialize();
	i1 = 0; i2 = 0; i3 = 0;
	for (i = 0; i <= nsteps; i++)
	{
		steps[i].assign(i1, i2, i3);
		switch (direction) {
		case 1: i1++; break;
		case 2: if (i % 2 == 0) i2++; else i1++; break;
		case 3: i2++; break;
		case 4: {int resd = i % 3;
			if (resd == 0) i1 += 1.47;
			else if (resd == 1) i1 += 1.44;
			else if (resd == 2) i1 += 1.37;
			else printf("line init, wrong");
		}break;
		default: printf("bad case in line_initialize() d\n"); exit(0); break;
		}
	}
	old_energy = GetEnergy();
	//printf("old energy: %lf\n", old_energy);
} // end line_initialize()

void Walk::deallocate()
{
	delete[] steps;
	delete[] ptime;
	delete[] igroup;
	delete[] shift;
};

void Walk::print(FILE *fptr)
// this prints the coordinates of the walk wrt the two lattice vectors
// not the actual x,y coords
// If npivot is > 0, it always prints the arrays ptime, igroup, shift.
// This is for debugging. It is incompatible with walk::scanf()
{
	int i;
	fprintf(fptr, "%d %d %lf \n", nsteps, generation, old_energy);
	for (i = 0; i <= nsteps; i++) steps[i].print(fptr);
	if (npivot>0)
	{
		fprintf(fptr, "npivot=%ld\n", npivot);
// 		fprintf(fptr, "iseg     ptime   igroup    shift \n");
// 		for (i = 0; i <= npivot; i++)
// 		{
// 			fprintf(fptr, "%4ld  %8ld   %2ld    ", i, ptime[i], igroup[i]);
// 			shift[i].print(fptr);
// 		} // end loop on i 
	} // end if (npivot>0)
} // end walk::print()

void Walk::scan(FILE *fptr)
{
	int i;
	fscanf(fptr, "%d %d %lf", &nsteps, &generation, &old_energy);
	for (i = 0; i <= nsteps; i++) steps[i].scan(fptr);
	initialize();
} // end walk::scan()

// void Walk::plot(FILE *fptr)
// // this prints out the actual x,y,z,... coordinates of the walk so it 
// // can be plotted
// {
// 	int i;
// 	if (npivot != 0)
// 	{
// 		printf("ERROR: walk::print() called with npivot!=0 \n");
// 		exit(1L);
// 	}
// 	for (i = 0; i <= nsteps; i++)
// 	fprintf(fptr, PRINTSTYLE, steps[i].coord_x(),steps[i].coord_y(), steps[i].coord_z());
// } // end walk::plot()

// Walk& Walk::operator=(Walk& w)
// {
// 	int i;
// 	if (this->get_npivot() != 0)
// 	{
// 		printf("ERROR in walk::=, npivot!=0 \n");
// 		exit(1);
// 	}
// 	(*this).initialize();
// 	nsteps = this->get_nsteps();
// 	niter = this->get_niter();
// 	for (i = 0; i <= nsteps; i++) steps[i] = this->steps[i];
// 	return(*this);
// }
// 
// float Walk::turn_frac()
// {
// 	// Finds fraction of times at which walk turns.
// 	// For hexagonal lattice we look at next nearest neighbors
// 	int i, count;
// 	count = 0;
// 
// 	for (i = 1; i<nsteps; i++)
// 	if (fabs((steps[i - 1] - steps[i + 1]).distance() - 2.) > 1.e-5) count++;
// 
// 	return(float(count) / float(nsteps - 1));
// } // end turn_frac

int Walk::pivot_strictly_saw(Proposal* prop)
// This is the version that changes i and j "simultaneously"
// If the pivot is accepted, this routine carries it out. Otherwise it is
// leaves the walk unchanged. 
// Routine returns count if the new walk  is self-avoiding, 
// returns -count if the pivot produces a self-intersection or 
// intersects the excluded region, 
// where count is the number of distance computations done.
// This return value is only used to study how long the routine takes. 
{
	OpMatrix *igroup_jseg, *igroup_iseg;
	int i, j, ip, jp, iseg, jseg, imin, imax, jmin, jmax;
	int separation, min_separation, sep_mod;
	MODELNAME origin, transi, transj, pp, stepsp, stepsi, stepsj, shift_jseg, shift_iseg;
	int count, changei_flag;
	int pivot_loc = prop->pivot_loc;
	OpMatrix* isym = &(prop->op);
	OpMatrix* invisym = &(prop->invop);
	// pivot is given by w[t] -> g (w[t]-w[pivot_loc])+w[pivot_loc]
	// this is equivalent to w[t] -> g w[t] + trans 
	// where trans= w[pivot_loc] - g w[pivot_loc]
	// stepsp=w[pivot_loc]

	origin.zero();
	iseg = find_segment(pivot_loc, npivot, ptime);
	stepsp.euclidean_op(steps + pivot_loc, shift + iseg, &igroup[iseg]);
	transi.euclidean_op(&stepsp, &origin, isym);
	transi = stepsp - transi;
	transj.euclidean_op(&stepsp, &origin, invisym);
	transj = stepsp - transj;
	count = 0;

	/////////////////////////////////////////////////////////////////////////////
	// 
	// "Simultaneous" changing of i and j. 
	// Given j<pivot_loc<i, we assume omega(jp)!=omega(ip) for j<jp<pivot_loc<ip<i
	// We then either decrease j or increase i.
	// 
	/////////////////////////////////////////////////////////////////////////////

	jmin = 0; jmax = pivot_loc - 1;
	imin = pivot_loc + 1; imax = nsteps;

	j = pivot_loc - 2;
	i = pivot_loc + 1;
	// if pivot_loc==0, j=-2 which would allow jp=-1
	if (j<jmin - 1) j = jmin - 1;

	if (pivot_loc>nsteps) i = imax + 1;

	if (!no_saw) while (i <= imax || j >= jmin)
	{
		// changei_flag=1 means we will increase i, =0 means we will increase j
		// We change the index that is closer to pivot_loc (if allowed).
		if (i - pivot_loc>pivot_loc - j) changei_flag = 0;
		else changei_flag = 1;
		if (i>imax) changei_flag = 0;
		if (j<jmin) changei_flag = 1;

		if (changei_flag)
		{// increase i. Need lower bound on distance from pivoted omega[i] to 
			// {omega[jp]: j<jp<pivot_loc}
			// This lower bound will be min_separation
			// stepsi is w[i]   
			iseg = find_segment(i, npivot, ptime);
			stepsi.euclidean_op(steps + i, shift + iseg, &igroup[iseg]);
			// pp is w[i] after pivot
			pp.euclidean_op(&stepsi, &transi, isym);
			min_separation = nsteps;
			jseg = npivot;   // can change to jseg=iseg; ?
			shift_jseg = shift[jseg] - pp;
			igroup_jseg = &igroup[jseg];
			for (jp = jmax; jp>j;) // note that j is decreased
			{
				if (ptime[jseg]>jp)
				{
#ifdef USE_FIND_SEGMENT
					jseg = find_segment(jp, npivot, ptime);
#else 
					while (ptime[jseg]>jp) jseg--;
#endif
					shift_jseg = shift[jseg] - pp;
					igroup_jseg = &igroup[jseg];
				}
				// stepsj is w[jp]   
				stepsj.euclidean_op(steps + jp, &shift_jseg, igroup_jseg);
				separation = int(stepsj.WellSeparate());
				count++;
				if (separation == 0) return(-count);
				if (separation >= min_separation)
				{
					jp -= separation - min_separation + 1;
				}
				else
				{
#ifdef USE_THIRD
					sep_mod = separation % 3;
					min_separation = (2 * separation + sep_mod) / 3;
					jp -= 1 + (separation - sep_mod) / 3;
#else 
					sep_mod = separation % 2;
					min_separation = (separation + sep_mod) / 2;
					jp -= 1 + (separation - sep_mod) / 2;
#endif
				}
			} // end loop on jp
			i += min_separation;
			if (i>imax + 1) i = imax + 1;
		} // end increase i 

		else
		{// decrease j. Need lower bound on distance from omega[j] to 
			// pivoted {omega[ip]: pivot_loc<ip<i}
			// Equivalently we can use a lower bound on distance from inverse
			// pivoted omega[j] and {omega[ip]: pivot_loc<ip<i}
			// This lower bound will be min_separation
			// stepsj is w[j]   
			jseg = find_segment(j, npivot, ptime);
			stepsj.euclidean_op(steps + j, shift + jseg, &igroup[jseg]);
			pp.euclidean_op(&stepsj, &transj, invisym);
			min_separation = nsteps;
			iseg = 0;
			shift_iseg = shift[iseg] - pp;
			igroup_iseg = &igroup[iseg];
			for (ip = imin; ip<i;) // note that i is increased
			{
				if (ptime[iseg + 1] <= ip) // check this
				{
#ifdef USE_FIND_SEGMENT
					iseg = find_segment(ip, npivot, ptime);
#else 
					while (ptime[iseg + 1] <= ip) iseg++;
#endif
					shift_iseg = shift[iseg] - pp;
					igroup_iseg = &igroup[iseg];
				}
				stepsi.euclidean_op(steps + ip, &shift_iseg, igroup_iseg);
				separation = int(stepsi.WellSeparate());
				count++;
				if (separation == 0) return(-count);
				if (separation >= min_separation)
				{
					ip += separation - min_separation + 1;
				}
				else
				{
#ifdef USE_THIRD
					sep_mod = separation % 3;
					min_separation = (2 * separation + sep_mod) / 3;
					ip += 1 + (separation - sep_mod) / 3;
#else 
					sep_mod = separation % 2;
					min_separation = (separation + sep_mod) / 2;
					ip += 1 + (separation - sep_mod) / 2;
#endif
				}
			} // end loop on ip
			j -= min_separation;
			if (j<jmin - 1) j = jmin - 1;
		} // end decrease j 
	} // end while 

	// If we reach this point the walk is self-avoiding. 
	// The pivot operations were not done to the walk. So we must do them.
	add_pivot(pivot_loc, isym, transi);
	return(count);

} // end pivot_strictly_saw()

void Walk::add_pivot(int pivot_loc, OpMatrix* isym, MODELNAME trans)
{
	int iseg, ipivot;
	MODELNAME pp;

	if (npivot>max_npivot - 1)
	{
		printf("number of implicit pivots in walk exceeds MAX_NPIVOT \n");
		exit(1);
	}

	// first, we add the pivot time to the list of pivot times. 
	iseg = find_segment(pivot_loc, npivot, ptime);
	if (pivot_loc != ptime[iseg])
	{
		ptime[npivot + 2] = ptime[npivot + 1];
		for (ipivot = npivot; ipivot >= iseg; ipivot--)
		{
			ptime[ipivot + 1] = ptime[ipivot];
			igroup[ipivot + 1] = igroup[ipivot];
			shift[ipivot + 1] = shift[ipivot];
		} // end loop on ipivot
		ptime[iseg + 1] = pivot_loc;
		npivot++;
		iseg++; // pivot will be applied to segments iseg to npivot
	}

	// second, we update igroup and shift
	for (ipivot = iseg; ipivot <= npivot; ipivot++)
	{
		pp = shift[ipivot];
		shift[ipivot].euclidean_op(&pp, &trans, isym);
		igroup[ipivot] = isym->dot(igroup[ipivot]);
	} // end loop on ipivot

} // add_pivot()

void Walk::simplify()
// carry out the pivot operations implicit in the walk, so npivot -> 0 
{
	OpMatrix*  igroup_ipivot;
	int ipivot, itime;
	MODELNAME shift_ipivot, pp;

	// even on the 0th segment there may be something to do 
	for (ipivot = 0; ipivot<npivot; ipivot++)
	{
		shift_ipivot = shift[ipivot];
		igroup_ipivot = &igroup[ipivot];
		for (itime = ptime[ipivot]; itime<ptime[ipivot + 1]; itime++)
		{
			pp = steps[itime];
			steps[itime].euclidean_op(&pp, &shift_ipivot, igroup_ipivot);
		}
	} // end loop on ipivot

	// npivot segment is different
	for (itime = ptime[npivot]; itime <= nsteps; itime++)
	{
		pp = steps[itime];
		steps[itime].euclidean_op(&pp, shift + npivot, &igroup[npivot]);
	}

	initialize();
} // end walk::simplify()

void Walk::run()
{
	generation++;
	//old_energy = GetEnergy();
	//printf("old energy: %lf\n", old_energy);
	int inner;
	// sign of accept flag indicates walk was accepted (>0) or rejected (<0) 
	// |accept_flag| is the number of tests done in the call to pivot routine
	int accept_flag = 0;
	double accept_ratio = 0.;

	////////////////////////////////////////////////////////////
	//      initialize running totals to zero                 //
	////////////////////////////////////////////////////////////
	int naccept = 0;
	
	int MCtrial = 0;
	while (MCtrial < 10)
	{
		MCtrial++;
		for (inner = 1; inner <= n_inner; inner++)
		{
			Proposal prop(this);
			accept_flag = pivot_strictly_saw(&prop);

			if (this->get_npivot() == nsimplify)
			{
				break;
			}
		} // end loop on inner

		double et = GetEnergy();
		if (et<old_energy)
		//if (RNG_NAME() <= (exp(et / old_energy) - 1) / (exp(1) - 1))
		{
			printf("The %d step accepted. new energy: %lf, old energy:%lf  \n", generation,et,old_energy);
			old_energy = et;
			simplify();
			naccept++;
			break;
		}
		else
		{
			//printf("The %d step denied. new energy: %lf, old energy:%lf  \n", generation, et, old_energy);
			initialize();
			naccept = 0;
		}		
	}
	// print the accept ratio and turn fraction.
	if (false) 
	{
		accept_ratio = 100 * double(naccept) / double(MCtrial);
		printf("%ld iterations, turn=%lf, MC accept ratio=%lf\n",
			inner, this->turn_frac(), accept_ratio);
	}
	// record the walk
	Record(); 
}

Walk::Walk(int tnsteps, char* tinit_walk_fname, int tnsimplify, int tn_inner, int tno_saw, int tmax_npivot, int tprint_freq) :
n_inner(tn_inner),
max_npivot(tmax_npivot),
no_saw(tno_saw),
nsteps(tnsteps), 
init_walk_fname(tinit_walk_fname), 
nsimplify(tnsimplify),
final_walk_fname("FinalWalk"),
data_fname("data"),
old_energy(-1)
{
	if(nsimplify == 0) nsimplify = int(sqrt(double(nsteps / 40)));
	if(max_npivot == 0) max_npivot = nsimplify+100;
#ifdef linux
	struct timeval tpstart;
	gettimeofday(&tpstart, NULL);
	srand48(tpstart.tv_usec);
	//srand48(floor(mytime()));
#endif
#ifdef _WIN32
	srand(unsigned int(time(NULL)));
#endif

	steps = new MODELNAME[nsteps + 1];

	ptime = new int[max_npivot + 2];// this may be one larger than needed 
	igroup = new OpMatrix[max_npivot + 1];
	shift = new MODELNAME[max_npivot + 1];

	////////////////////////////////////////////////////////////
	//  initialize walk to a line or  read walk from file     //
	////////////////////////////////////////////////////////////

	// if filename for initial walk starts with a 0 we generate a line
	// for the initial walk. direction is its direction:
	// 1=horizontal, 2=45 degs, 3=vertical, 4=60 degs 
	int direction = 4;
	if (init_walk_fname[0] == '0') this->line_initialize(direction);
	else
	{
		FILE *fptr = fopen(init_walk_fname, "r");
		if (fptr == NULL)
		{
			printf("FAILURE to open walk file %s, exiting \n", init_walk_fname);
			exit(0);
		}
		this->scan(fptr);
		fclose(fptr);
	}

	// checks that initial walk is ok
 	printf("Initial walk has turn fraction %lf \n", this->turn_frac());
}

double Walk::turn_frac()
{
	// Finds fraction of times at which walk turns.
	// For hexagonal lattice we look at next nearest neighbors
	long i, count;
	count = 0;
	for (i = 1; i<nsteps; i++)
	if (fabs((steps[i - 1] - steps[i + 1]).distance() - 2.) > 0.267) count++;
	return(double(count) / double(nsteps - 1));
} // end turn_frac

double Walk::GetEnergy()
{
	double r = 99;
	double etemp = 0;
	for (int i = 0; i <= nsteps; i++)
	{	
		if (i % 3 == 2) i++;
		for (int j = 2; j <= nsteps;j=j+3)
		{
			r = (step_rval(i) - step_rval(j)).distance();
			if (1<r && r<10)
			{
				etemp += -1.0 / r;
			}
			if (r < 1) etemp += -1;
		}
	}
	return etemp;
}

void Walk::Record()
{
	FILE *fptr;
	// 			// record the "data"
	// 			fptr = fopen(data_fname, "a");
	// 			fprintf(fptr, "%14.10f\n", 0.0);
	// 			fclose(fptr);

	// record the walk itself
	fptr = fopen(final_walk_fname, "w");
	this->print(fptr);
	fclose(fptr);
}
