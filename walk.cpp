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

Sphere Walk::GetStepi(int i) //the old step_rval
// This routine computes the ith point on the walk from the complicated 
// data structure. Often we duplicate this code for speed rather than call this
{
	Sphere temp;
	int iseg;
	iseg = find_segment(i, npivot, ptime);
	temp.euclidean_op(steps + i, shift + iseg, &igroup[iseg]);
	return(temp);
} //

void Walk::clean_pivot()
// NB this does not create a walk. It only initializes things used by the 
// implicit pivot data structure
{
	npivot = 0;
	igroup[0].identity();
	shift[0].zero();
	// segment number iseg corresponds to times [ptime[iseg],ptime[iseg+1])
	// So when npivot=0 we should have : 
	ptime[0] = 0;
	ptime[1] = nsteps + 1;
}

void Walk::line_initialize(int direction)
{ 
	int i;
	double i1, i2, i3;
	clean_pivot();
	i1 = 0; i2 = 0; i3 = 0;
	for (i = 0; i <= nsteps; i++)
	{
		steps[i].assign(i1, i2, i3, RADIUS, RIGID);
		//if(i == 5) steps[i].assign(i1, i2, i3, 0.2, 20);
		switch (direction) {
		case 1: i1++; break;
		case 2: if (i % 2 == 0) i2++; else i1++; break;
		case 3: i2++; break;
		default: printf("bad case in line_initialize() d\n"); exit(0); break;
		}
	}
	old_energy = GetEnergy();
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
		printf("Error, npivot=%ld while printing\n", npivot);
	}
} // end walk::print()

void Walk::scan(FILE *fptr)
{
	int i;
	fscanf(fptr, "%d %d %lf", &nsteps, &generation, &old_energy);
	for (i = 0; i <= nsteps; i++) steps[i].scan(fptr);
	clean_pivot();
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
	OPERATION_NAME *pgroup_jseg, *pgroup_iseg;
	int i, j, ip, jp, iseg, jseg, imin, imax, jmin, jmax;
	int separation, min_separation, sep_mod;
	Sphere  pp, stepsp, stepsi, stepsj;
	GPoint<double> origin, transi, transj, shift_jseg, shift_iseg;
	int count, changei_flag;
	int pivot_loc = prop->pivot_loc;
	OPERATION_NAME* poper = &(prop->op);
	OPERATION_NAME* pinvoper = &(prop->invop);
	// pivot is given by w[t] -> g (w[t]-w[pivot_loc])+w[pivot_loc]
	// this is equivalent to w[t] -> g w[t] + trans 
	// where trans= w[pivot_loc] - g w[pivot_loc]
	// stepsp=w[pivot_loc]

	origin.zero();
	iseg = find_segment(pivot_loc, npivot, ptime);
	stepsp.euclidean_op(steps + pivot_loc, shift + iseg, &igroup[iseg]);
	transi= GPoint<double>(stepsp.x,stepsp.y,stepsp.z).rotation(poper);
	transi.x = stepsp.x - transi.x;
	transi.y = stepsp.y - transi.y;
	transi.z = stepsp.z - transi.z;
	transj= GPoint<double>(stepsp.x,stepsp.y,stepsp.z).rotation(pinvoper);
	transj.x = stepsp.x - transj.x;
	transj.y = stepsp.y - transj.y;
	transj.z = stepsp.z - transj.z;
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
			pp.euclidean_op(&stepsi, &transi, poper);
			min_separation = nsteps;
			jseg = npivot;   // can change to jseg=iseg; ?
			shift_jseg = shift[jseg];
			pgroup_jseg = &igroup[jseg];
			for (jp = jmax; jp>j;) // note that j is decreased
			{
				if (ptime[jseg]>jp)
				{
					//jseg = find_segment(jp, npivot, ptime);
					while (ptime[jseg]>jp) jseg--;
					shift_jseg = shift[jseg];
					pgroup_jseg = &igroup[jseg];
				}
				// stepsj is w[jp] - w[i]   
				stepsj.euclidean_op(steps + jp, &shift_jseg, pgroup_jseg);
				separation = stepsj.WellSeparate(pp);
				count++;
				if (separation == 0) return(-count);
				if (separation >= min_separation)
				{
					jp -= separation - min_separation + 1;
				}
				else
				{
					//sep_mod = separation % 3;
					//min_separation = (2 * separation + sep_mod) / 3;
					//jp -= 1 + (separation - sep_mod) / 3; 
					sep_mod = separation % 2;
					min_separation = (separation + sep_mod) / 2;
					jp -= 1 + (separation - sep_mod) / 2;
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
			// stepsj is w[j] before this pivot
			jseg = find_segment(j, npivot, ptime);
			stepsj.euclidean_op(steps + j, shift + jseg, &igroup[jseg]);
			pp.euclidean_op(&stepsj, &transj, pinvoper);
			min_separation = nsteps;
			iseg = 0;
			shift_iseg = shift[iseg];
			pgroup_iseg = &igroup[iseg];
			for (ip = imin; ip<i;) // note that i is increased
			{
				if (ptime[iseg + 1] <= ip) // check this
				{
					//iseg = find_segment(ip, npivot, ptime);
					while (ptime[iseg + 1] <= ip) iseg++;
					shift_iseg = shift[iseg];
					pgroup_iseg = &igroup[iseg];
				}
				stepsi.euclidean_op(steps + ip, &shift_iseg, pgroup_iseg);
				//separation = int(stepsi.WellSeparate());
				separation = stepsi.WellSeparate(pp);
				count++;
				if (separation == 0) return(-count);
				if (separation >= min_separation)
				{
					ip += separation - min_separation + 1;
				}
				else
				{
					//sep_mod = separation % 3;
					//min_separation = (2 * separation + sep_mod) / 3;
					//ip += 1 + (separation - sep_mod) / 3;
					sep_mod = separation % 2;
					min_separation = (separation + sep_mod) / 2;
					ip += 1 + (separation - sep_mod) / 2;
				}
			} // end loop on ip
			j -= min_separation;
			if (j<jmin - 1) j = jmin - 1;
		} // end decrease j 
	} // end while 

	// If we reach this point the walk is self-avoiding. 
	double new_energy = GetEnergy();
	if (AcceptOrNot(new_energy,old_energy))
	{
		//printf("The %d step accepted. new energy: %lf, old energy:%lf  \n", generation,et,old_energy);
		old_energy = new_energy;
		add_pivot(pivot_loc, poper, transi);
		return(count);
	}
	else
	{
		//printf("The %d step denied. new energy: %lf, old energy:%lf  \n", generation, et, old_energy);
		//if (MCtrial == MaxTrial) printf("In generation %d, totally %d trial done and all denied.", generation, MaxTrial);
		return(-count);
	}	
} // end pivot_strictly_saw()

void Walk::add_pivot(int pivot_loc, OPERATION_NAME* poper, GPoint<double> trans)
{
	int iseg, ipivot;
	GPoint<double> pp;

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
		shift[ipivot] = pp.rotation(poper)+trans;
		igroup[ipivot] = poper->dot(igroup[ipivot]);
	} // end loop on ipivot

} // add_pivot()

void Walk::simplify()
// carry out the pivot operations implicit in the walk, so npivot -> 0 
{	
	OPERATION_NAME*  pOper_ipivot;
	int ipivot, itime;
	Sphere pp;
	GPoint<double> shift_ipivot;

	// even on the 0th segment there may be something to do 
	for (ipivot = 0; ipivot<npivot; ipivot++)
	{
		shift_ipivot = shift[ipivot];
		pOper_ipivot = &igroup[ipivot];
		for (itime = ptime[ipivot]; itime<ptime[ipivot + 1]; itime++)
		{
			pp = steps[itime];
			steps[itime].euclidean_op(&pp, &shift_ipivot, pOper_ipivot);
		}
	} // end loop on ipivot

	// npivot segment is different
	for (itime = ptime[npivot]; itime <= nsteps; itime++)
	{
		pp = steps[itime];
		steps[itime].euclidean_op(&pp, shift + npivot, &igroup[npivot]);
	}

	clean_pivot();
} // end walk::simplify()

void Walk::GoOneStep(int stepnum, bool isrecord)
{
	generation++;
	int inner = 0, accept_flag=0;
	for (inner = 1; inner <= stepnum; inner++)
	{
		Proposal prop(this);
		//printf("%d \n",prop.pivot_loc);
		accept_flag = pivot_strictly_saw(&prop);
	   if (accept_flag >= 0) SAWaccept++;
		if (npivot >= nsimplify)
		{
			simplify();
		}
	} // end loop on inner then #nsimplify successful pivot is prepared	
	simplify();
	if (npivot != 0)printf("Error. Npivot is not zero at the end.\n");

	// print the accept ratio and turn fraction
	if (false) 
	{
		//printf("%d generation, turn=%lf, MC accept ratio=%lf\n",generation, turn_frac(), 100.0 / double(MCtrial));
	}
	// record the walk 
	if(isrecord) Record(); 
}

Walk::Walk(int tnsteps, const char* tinit_walk_fname, int tn_inner, const char* tfinal_walk_fname, int tnsimplify, int tno_saw, int tmax_npivot):
n_inner(tn_inner),
max_npivot(tmax_npivot),
no_saw(tno_saw),
nsteps(tnsteps), 
init_walk_fname(tinit_walk_fname), 
nsimplify(tnsimplify),
final_walk_fname(tfinal_walk_fname),
data_fname("data"),
old_energy(-1),
generation(0),
SAWaccept(0)
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

	steps = new Sphere[nsteps + 1];

	ptime = new int[max_npivot + 2];// this may be one larger than needed 
	igroup = new OPERATION_NAME[max_npivot + 1];
	shift = new GPoint<double>[max_npivot + 1];

	////////////////////////////////////////////////////////////
	//  initialize walk to a line or  read walk from file     //
	////////////////////////////////////////////////////////////

	// if filename for initial walk starts with a 0 we generate a line
	// for the initial walk. direction is its direction:
	// 1=horizontal, 2=45 degs, 3=vertical
	int direction = 2;
	if (init_walk_fname[0] == '0') line_initialize(direction);
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
 	//printf("Initial walk has turn fraction %lf \n", this->turn_frac());
}

double Walk::turn_frac()
{
	// Finds fraction of times at which walk turns.
	// For hexagonal lattice we look at next nearest neighbors
	long i, count;
	count = 0;
	for (i = 1; i<nsteps; i++)
	if (fabs(steps[i - 1].distance(steps[i + 1]) - 2.) > 0.267) count++;
	return(double(count) / double(nsteps - 1));
} // end turn_frac

void Walk::Record()
{
	FILE *fptr;
	char buffer[50]; // <- danger, only storage for 256 characters.
	sprintf(buffer, "%s_%d", data_fname,nsteps);
	// record the "data"
	double endnorm = GetStepi(nsteps).topoint().norm();
	fptr = fopen(buffer, "a");
    fprintf(fptr,"%14.10f    %14.10f \n",GetRg2(),endnorm*endnorm);
    fclose(fptr);

//	// record the walk itself
//	fptr = fopen(final_walk_fname, "w");
//	this->print(fptr);
//	fclose(fptr);
}

void Walk::run(int outer_steps, int discard)
{
	GoOneStep(discard,false);
	for (int i = 0; i < outer_steps; i++)
	{
		GoOneStep(n_inner,true);
	}
	printf("%d + %d * %d MCSs, turn=%lf, accept ratio=%lf\n", discard, outer_steps, n_inner, turn_frac(), SAWaccept*100.0/double(outer_steps*n_inner+discard));
	//printf("%d\n",SAWaccept);
	if(npivot!=0) printf("error, npivot at the end of outer run.\n");
}

double Walk::GetRg2()
{
  GPoint<double> rc = GetStepi(0).topoint();
  double dis = 0;
  double rg2 = 0;
  
  for (int i=1;i<=nsteps;i++)
  {
  	rc = rc + GetStepi(i).topoint();
  } 
  rc /= (double)(nsteps+1);
  
  for (int i=0;i<=nsteps;i++)
  {
  	dis = (GetStepi(i).topoint()-rc).norm();
  	rg2 += dis*dis;
  } 
  rg2 /= (nsteps+1);
  return(rg2);
}

int Walk::GetAutocorrelation(int n, unsigned long int intersteps)
{
	std::string datname;
	std::stringstream n_temp;
	n_temp<<nsteps;
	
	datname = "savedata/auto/ac_fkt_"+n_temp.str();
	std::ofstream fout(datname.c_str(), std::ios::trunc);
	
	datname = "savedata/auto/Parameters_"+n_temp.str();
	std::ofstream fout1(datname.c_str(),std::ios::trunc);
	
	//unsigned long int intersteps = 1;
	double tmpRg2;
	std::vector<double> Rg2;
	double Rg2m = 0.;
	double C0 = 0.;
	double C = 0.;
	double tao_int = 0.5;
	unsigned M = 0;
	
	std::string line;
	
	// Read in the values for Rg2
	for(int k=0;k<n;k++)
	{
		GoOneStep(intersteps,false);
		tmpRg2 = GetRg2();		
		Rg2.push_back(tmpRg2);
		Rg2m += tmpRg2;
	}
	
	if (n != Rg2.size()) {std::cerr << "Alert! " << n << std::endl;}
	
	//std::cout << "Rg2[0] = " << Rg2[0] << std::endl;
	fout1 <<"Rg2[0]      " << Rg2[0] <<"\n";
	
	unsigned nu = n - n/4;
	
	Rg2m /= n;
	
	// Calculate C(0)
	for (int i=0; i<n; ++i)
		C0 += (Rg2[i]-Rg2m)*(Rg2[i]-Rg2m);
	
	C0 /= n;
	
	//std::cout << "Number of data: " << n << std::endl;
	fout1 << "NumberOfData      " << n <<"\n";
	//std::cout << "Mean radius of gyration: " << Rg2m << std::endl;
	fout1 << "MeanRadiusOfGyration      " << Rg2m<<"\n";
	fout << "# Autocorrelation function C(t); C(0)="<< C0 <<"\n";
	fout << "# t\tC(t)\ttao_int\n";
	fout << "0    1    0.5\n";
	//fout << "# t\tC0*C(t)\tC(t)\ttao_int\t10*tao_int\n";
	//fout << "0 " << C0 << " 1 0.5 5.0\n";
	
	
	bool tmp = false;
	
	// Calculate C(t) for t>0
	for (unsigned t=1; t<nu; ++t) {
		
		C = 0;
		
		// Go through all the data
		for (unsigned i=0; i<(n-t); ++i)
			C += (Rg2[i]-Rg2m)*(Rg2[i+t]-Rg2m);
			
		C /= (n-t)*C0;
		
		// Integrate w/ very simple rule
		tao_int += intersteps*C;
		
		//fout << t*intersteps << "\t" << C*C0 << "\t" << C << "\t" << tao_int << "\t" << 10*tao_int << std::endl;
		fout << t*intersteps << "    " << C << "    " << tao_int << std::endl;
		if ( (tmp==false) && ((t*intersteps) >= 10*tao_int) ) {
			M = (unsigned int)tao_int;
			//std::cout << "the number of step truncated " <<t*intersteps<<endl;
			//std::cout << "tao_int = " << M << std::endl;
			fout1<<  "tao_int      " << M << "\n";
			tmp = true;
// 			break;
		}
	}
	
	fout.close();
	fout1.close();
	
	return M;
}

//////////////////////////////////////////////////////////////////////////
// This is the part for MC
//////////////////////////////////////////////////////////////////////////

double Walk::GetEnergy()
{
/*	double r = 99;
	double etemp = 0;
	for (int i = 0; i <= nsteps; i++)
	{
		if (i % 3 == 2) i++;
		for (int j = 2; j <= nsteps; j = j + 3)
		{
			r = (GetStepi(i) - GetStepi(j)).distance();
			if (1 < r && r < 10)
			{
				etemp += -1.0 / r;
			}
			if (r < 1) etemp += -1;
		}
	}
	return etemp;*/

	return -1; // now is a empty function
// 	double etemp = 0;
// 	for (int i = 1; i < nsteps; i++)
// 	{
//         etemp += steps[i].k*((GetStepi(i-1)-GetStepi(i)).dot(GetStepi(i+1)-GetStepi(i)));
// 		//printf("%lf\n",((GetStepi(i-1)-GetStepi(i)).dot(GetStepi(i+1)-GetStepi(i))));
// 	}
//     return etemp;
}

bool Walk::AcceptOrNot(double newE, double oldE)
{
	return true;
	//return(newE < oldE);
	//return(RNG_NAME() <= (exp(newE / oldE) - 1) / (exp(1) - 1));
    //return(RNG_NAME() <= exp(oldE-newE));
}
