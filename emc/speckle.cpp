#include "pol_montecarlo.h"
#include "netcdfcpp.h"
#include <math.h>
#include <iostream.h>
#include <unistd.h>

// Remember all the length and time units are in mfp=lsa.
const double n_water=1.33;
const double n_sct=1.59;
const double lambda=0.515;	// micro
const double dia=0.46;		// micro
const double xx=pi*dia/lambda*n_water; // size parameter of the Mie sphere
const double mr=n_sct/n_water, mi=0.0;	// m=mr-i*mi is the refractive index
const double L=40;		// slab thickness in lsa
const double HW=10;		// the half width of the window in lsa
const int NN=100;		// NN grids on the backscattering window
const double KK=45;		// the wavenumber in lsa
const int MAXSCT=1000;
double S[NN][NN][4];
double EE[NN][NN][4];
static long prephotons=0;	// number of photons used in previous simulations
static long photons=10000;      // number of photons to use in the current simulation
				// 45mins in the mode 2

// convert parallel E1 and perpendicular E2 to Stokes vector
static inline void EtoS(dcmplx E1, dcmplx E2, double S[4]) {
    // we comment this out if we only collect Ex and Ey in the detector
    // assert( fabs(norm(E1) + norm(E2) - 1) < 1e-8 );
    S[0] = norm(E1) + norm(E2);
    S[1] = norm(E1) - norm(E2);
    S[2] = 2*(conj(E1)*E2).real();
    S[3] = 2*(conj(E1)*E2).imag();
}

static int write_results(const char* fname)
{
    NcFile nc (fname, NcFile::Replace); // Create and leave in define mode

    // Check if the file was opened
    if (! nc.is_valid()) {
	cerr << "can't create netCDF file " << fname << "\n";
	return 0;
    }

    // Create dimensions
    NcDim* nn = nc.add_dim("NN", NN);
    NcDim* four = nc.add_dim("FOUR", 4);

    // Create variables and their attributes
    NcVar* c = nc.add_var("S", ncDouble, nn, nn, four);
    c->add_att("long_name", "Stokes vector");

    NcVar* c1 = nc.add_var("EE", ncDouble, nn, nn, four);
    c1->add_att("long_name", "ExRe ExIm EyRe EyIm");

    // Global attributes
    nc.add_att("NN", NN);
    nc.add_att("HW", HW);
    nc.add_att("xx", xx);
    nc.add_att("mr", mr);
    nc.add_att("mi", mi);
    nc.add_att("L", L);
    nc.add_att("photons", photons);

    // Start writing dcata, implictly leaves define mode
    c->put(&S[0][0][0], c->edges());
    c1->put(&EE[0][0][0], c1->edges());

    // close of nc takes place in destructor
    return 1;
}

static int load_results(const char* fname)
{
    NcFile nc (fname, NcFile::ReadOnly); // Create and leave in define mode

    // Check if the file was opened
    if (! nc.is_valid()) {
	cerr << "can't open netCDF file " << fname << "\n";
	return 0;
    }

    // Get attribute
    NcAtt* a=nc.get_att("photons");
    prephotons = a->values()->as_long(0);

    // Get dimensions
    NcDim* nn = nc.get_dim("NN");
    long n = nn->size();
    long four = 4;

    // Get variables
    NcVar* c = nc.get_var("S");
    NcVar* c1 = nc.get_var("EE");

    // Push values
    c->get(&EE[0][0][0], n, n, four);
    c1->get(&S[0][0][0], n, n, four);

    // display the current state
    printf("Previous simulations used %ld photons.\n", prephotons);

    nc.close();
    return 1;
}

int main (int argc, char** argv)
{
    dcmplx Ei1, Ei2;		// incident E field
    double Si[4];
    dcmplx Ex, Ey;		// outgoing E field
    double So[4];
    dcmplx Ed1, Ed2, phase;     // used by pointEstimator
    double td, deposit, Q[3][3];
    int ix, iy;
    double T, R, Lost;		// total transmission, reflection and lost photons
    char *seedfname=new char[50];
    char *ncfname=new char[50];

    dcmplx m(mr, mi);
    scatterer sct(xx, m);
    photonPacket *ph;

    if (argc == 1)
        sprintf(seedfname, "speckle.seed");
    else
      	sprintf(seedfname, "speckle%s.seed", argv[1]);

    if ( access(seedfname, R_OK) == 0 )
	ph = new photonPacket(&sct, seedfname);
    else
	ph = new photonPacket(&sct);

    for (int cno = 0; cno < 4; cno++) {
	T = R = Lost = 0;	// reset counters
	for (ix=0; ix<NN; ix++)
	    for (iy=0; iy<NN; iy++)
		for (int i=0; i<4; i++)
		    S[ix][iy][i] = 0;

	if (cno == 0) {
	    Ei1 = 1; Ei2 = 0;	// Horizontal polarized
	    if (argc == 1)
                sprintf(ncfname, "speckleH.nc");
            else
            	sprintf(ncfname, "speckleH%s.nc", argv[1]);
	} else if (cno == 1) {
	    Ei1 = 0; Ei2 = 1;	// Vertical polarized
	    if (argc == 1)
                sprintf(ncfname, "speckleV.nc");
            else
            	sprintf(ncfname, "speckleV%s.nc", argv[1]);
	} else if (cno == 2) {
	    Ei1 = 1/sqrt(2); Ei2 = 1/sqrt(2); // 45 degree polarized
	    if (argc == 1)
                sprintf(ncfname, "speckleP.nc");
            else
            	sprintf(ncfname, "speckleP%s.nc", argv[1]);
	} else {
	    Ei1 = 1/sqrt(2); Ei2 = dcmplx(0, 1/sqrt(2)); // right circular
	    if (argc == 1)
                sprintf(ncfname, "speckleR.nc");
            else
            	sprintf(ncfname, "speckleR%s.nc", argv[1]);
	}
	if ( access(ncfname, R_OK) == 0 ) {
	    printf("Loaded %s and continue...\n", ncfname);
	    load_results(ncfname);
	} else
	    printf("Start over...\n");

	EtoS(Ei1, Ei2, Si);
	printf("Incindent Stokes vector: [%f %f %f %f]\n", Si[0], Si[1], Si[2], Si[3]);

	for (int i = 1; i <= photons; i++) {
	    ph->launch (Ei1, Ei2);

	    if (i % 10000 == 0) fprintf(stderr, "%d\n", i);

	    while (1) {
		ph->move ();

		if (ph->nsct > MAXSCT) { Lost += ph->weight; break; }
		if (ph->z > L) { T += ph->weight; break; }
		if (ph->z < 0) { R += ph->weight; break; }

		// find when and where light crosses the boundary
		ph->pointestimator(&td, &deposit, Q, &Ed1, &Ed2, 0, 0, 0, -1);
		ix = (int)floor((ph->x+HW)/HW*NN/2);
		iy = (int)floor((ph->y+HW)/HW*NN/2);
		if (ix < 0) ix = 0;
		if (ix >= NN) ix = NN - 1;
		if (iy < 0) iy = 0;
		if (iy >= NN) iy = NN - 1;

		// compute Ex and Ey from Ed
		phase = dcmplx(cos(KK*td), sin(KK*td));
		Ex = (Ed1*Q[0][0] + Ed2*Q[1][0])*phase;
		Ey = (Ed1*Q[0][1] + Ed2*Q[1][1])*phase;

		// accumulate Ex and Ey
		EE[iy][ix][0] += Ex.real()*sqrt(deposit);
		EE[iy][ix][1] += Ex.imag()*sqrt(deposit);
		EE[iy][ix][2] += Ey.real()*sqrt(deposit);
		EE[iy][ix][3] += Ey.imag()*sqrt(deposit);

		//
		// IMPORTANT: the Stokes vector of the outgoing light refers to a local coordinate system (-x, y, -z) in reflection
		//            The field (Ex, Ey, Ez) in the global system must be transformed to (-Ex, Ey, -Ez)
		EtoS(-Ex, Ey, So);

		// We use [iy][ix] agreeing with the convention in Jessica's iquv
		S[iy][ix][0] += So[0]*deposit;
		S[iy][ix][1] += So[1]*deposit;
		S[iy][ix][2] += So[2]*deposit;
		S[iy][ix][3] += So[3]*deposit;

		ph->scatter();
	    }
	}

	printf("Case#, T, R, Lost=%d %f %f %f\n", cno, T, R, Lost);
	photons += prephotons;
	write_results(ncfname);

	photons -= prephotons;	// revert photons to the inital value for the next case
    }

    ph->saveRand(seedfname);
    delete ph;
    delete seedfname;
    delete ncfname;
}
