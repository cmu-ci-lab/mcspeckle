#include "pol_montecarlo.h"
#include "netcdfcpp.h"
#include <math.h>
#include <iostream.h>
#define ranf mtrand->drand

// Remember all the length and time units are in mfp=lsa.

const double lambda=0.6328;	// micro
const double dia=2.0;		// micro
const double n_water=1.33;
const double n_sct=1.59;
const double xx=pi*dia/lambda*n_water; // size parameter of the Mie sphere
const double mr=n_sct/n_water, mi=0.0; // m=mr-i*mi is the refractive index
const double L=4;		// slab thickness in lsa
const double HW=10;		// the half width of the window in lsa
const int NN=100;		// NN grids on the backscattering window
const int MAXSCT=500;		// kill the photon after this number of scattering
double M[NN][4][4];
const int photons=10000000*30;	// 10 mins for 10000000.


// convert parallel E1 and perpendicular E2 to Stokes vector
static inline void EtoS(dcmplx E1, dcmplx E2, double S[4]) {
    dcmplx E12;
    S[0] = 1;
    S[1] = norm(E1) - norm(E2);
    E12 = conj(E1)*E2;
    S[2] = 2*E12.real();
    S[3] = 2*E12.imag();
}


// M(rho, phi=0) = So*Si*D
// with D=[[1,0,0,0]; [0,3-cos4phi,sin4phi,0]; [0,sin4phi,3+cos4phi,0,]; [0,0,0,4]]
static inline void mullerDeposit(double So[4], double Si[4], double phi, double weight, double M[4][4]) {
    double sin4phi=sin(4*phi), cos4phi=cos(4*phi);
    double D00=1, D11=3-cos4phi, D12=sin4phi, D21=sin4phi, D22=3+cos4phi, D33=4;
    for (int i=0; i<4; i++) {
	M[i][0] += So[i]*Si[0]*D00*weight;
	M[i][1] += So[i]*(Si[1]*D11 + Si[2]*D21)*weight;
	M[i][2] += So[i]*(Si[1]*D12 + Si[2]*D22)*weight;
	M[i][3] += So[i]*Si[3]*D33*weight;
    }
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
    NcVar* c = nc.add_var("M", ncDouble, nn, four, four);
    c->add_att("long_name", "Muller matrix at phi=0");

    // Global attributes
    nc.add_att("NN", NN);
    nc.add_att("HW", HW);
    nc.add_att("xx", xx);
    nc.add_att("mr", mr);
    nc.add_att("mi", mi);
    nc.add_att("L", L);
    nc.add_att("photons", photons);

    // Start writing dcata, implictly leaves define mode
    c->put(&M[0][0][0], c->edges());

    // close of nc takes place in destructor
    return 1;
}


int main ()
{
    dcmplx Eix, Eiy;		// incident E field
    double alpha, beta, Si[4];	// the incident Stokes vector
    dcmplx Ed1, Ed2, Eox, Eoy;	// outgoing E field
    double td, deposit, Q[3][3], So[4];
    dcmplx Ei1, Ei2, Eo1, Eo2;	// incident and outgoing E with respect to the phi=0 plane
    double dt, x, y, rho, phi;
    int irho, nsct;
    double T=0, R=0, Lost=0;	// total transmission and reflection
    mt_prng *mtrand=new mt_prng();

    dcmplx m(mr, mi);
    scatterer sct(xx, m);

    photonPacket *ph = new photonPacket(&sct);

    for (int i = 1; i <= photons; i++) {
	alpha = ranf()*pi/2;
	beta = ranf()*pi;
	// Eix, Eiy are the parallel (x) and perpendicular (y) electric fields
	Eix = dcmplx(cos(alpha)*cos(beta), -cos(alpha)*sin(beta));
	Eiy = dcmplx(sin(alpha)*cos(beta), sin(alpha)*sin(beta));
	ph->launch (Eix, Eiy);

	if (i % 10000 == 0) fprintf(stderr, "%d\n", i);

	nsct = 0;
	while (1) {
	    ph->move ();

	    // use pointEstimator
	    if (nsct > MAXSCT) break;
	    if (ph->z > L) { T += ph->weight; break; }
	    if (ph->z < 0) { R += ph->weight; break; }

	    dt = (ph->t-ph->told)/(ph->z-ph->zold)*(0.-ph->zold);
	    x = ph->xold + ph->P[2][0]*dt;
	    y = ph->yold + ph->P[2][1]*dt;
	    rho = sqrt(x*x + y*y);

	    irho = (int)floor(rho/HW*NN);
	    if (irho >= NN) irho = NN - 1;

	    // crossing the z=0 boundary
	    ph->pointestimator(&td, &deposit, Q, &Ed1, &Ed2, 0, 0, 0, -1);

	    // compute Eox and Eoy from Ed
	    Eox = Ed1*Q[0][0] + Ed2*Q[1][0];
	    Eoy = Ed1*Q[0][1] + Ed2*Q[1][1];

	    // We need to rotate both (Eix, Eiy) and (Eox, Eoy) such that the electric field is with respect to the phi=0 plane in the gloabl coordinate system
	    phi = atan2(y, x);
	    Ei1 = Eix*cos(phi) + Eiy*sin(phi);
	    Ei2 = -Eix*sin(phi) + Eiy*cos(phi);
	    Eo1 = Eox*cos(phi) + Eoy*sin(phi);
	    Eo2 = -Eox*sin(phi) + Eoy*cos(phi);

	    //
	    // IMPORTANT: outgoing stokes vector is defined in the local coordinate system
	    //
	    EtoS(Ei1, Ei2, Si);
	    EtoS(-Eo1, Eo2, So);

	    // accumulate Muller matrix.
	    mullerDeposit(So, Si, phi, deposit, M[irho]);

	    ph->scatter ();
	    nsct ++;
	    }
    }

    printf("T, R, Lost=%f %f %f\n", T, R, Lost);
    write_results("rmuller.nc");
    delete ph;
}
