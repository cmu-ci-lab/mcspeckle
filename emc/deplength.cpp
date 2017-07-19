#include "pol_montecarlo.h"
#include <math.h>

// Remember all the length and time units are in mfp=lsa.
//
// A test case for slabs to obtain the depolarization length
// incident CW plane wave
//
// Polystyrene sphere diamter 0.5micro
// refractive index: 1.59 with background refractive index 1.33
// incident wavelength 632.8nm
//                          X
//                          |              |
//			    |              |
//			    |              |
//			    |              |
//		   Analyzer |              |
//			+   |              |
//			+   |O             |
//	 incident-------+---|--------------+------ Z
//			+   |              |
//			+   |       L      |
//
//
//-

const double lambda=0.633;	// micro
const double n_water=1.33;
const double dia=0.94;		// diameter in micro
const double x=4; //pi*dia/lambda*n_water;
const double mr=1.2, mi=0.0;	// m=mr-i*mi is the refractive index

int main ()
{
    dcmplx Ex, Ey, Ez, Ed1, Ed2;
    double td, deposit, Q[3][3];
    dcmplx m(mr, mi);
    scatterer sct(x, m);

    double L, Ipar, Iper, IxyR, IxyI;
    int photons = 100000*10;

    photonPacket *ph = new photonPacket(&sct);

    printf("# L, photons, Ipar, Iper, IxyR, IxyI, DOP, DOCP\n");

    for (L=1; L<10; L+=2) {
	Ipar = Iper = IxyR = IxyI = 0;

	for (int i = 1; i <= photons; i++) {
	    // launch a linear polarized light
	    ph->launch ();
	    // launch a circular polarized light
	    //ph->launch (dcmplx(1/sqrt(2), 0), dcmplx(0, 1/sqrt(2)));

	    //if (i % 1000 == 0) fprintf(stderr, "%d\n", i);

	    while (1) {
	    // either choice is OK.

#if 0 	

		ph->move ();
		if (ph->z < 0) break;
		if (ph->z > L) {
		    // NOTE: 0.99 corresponds to a detecting half angle of 8 degree.
		    if (ph->P[2][2] > 0.995 && ph->nsct > 0) {
			Ex = ph->E1*ph->P[0][0] + ph->E2*ph->P[1][0];
			Ey = ph->E1*ph->P[0][1] + ph->E2*ph->P[1][1];
			Ez = ph->E1*ph->P[0][2] + ph->E2*ph->P[1][2];

			Ipar += norm(Ex)*ph->weight;
			Iper += norm(Ey)*ph->weight;
			IxyR += (Ex*conj(Ey)).real()*ph->weight;
			IxyI += (Ex*conj(Ey)).imag()*ph->weight;
		    }
		    break;
		}

#else

		if (ph->z < 0 || ph->z > L) break;
		if (ph->nsct >= 0) {
		    ph->pointestimator(&td, &deposit, Q, &Ed1, &Ed2, L, 0, 0, 1);
		    Ex = Ed1*Q[0][0] + Ed2*Q[1][0];
		    Ey = Ed1*Q[0][1] + Ed2*Q[1][1];
		    Ez = Ed1*Q[0][2] + Ed2*Q[1][2];

		    Ipar += norm(Ex)*deposit;
		    Iper += norm(Ey)*deposit;
		    IxyR += (Ex*conj(Ey)).real()*deposit;
		    IxyI += (Ex*conj(Ey)).imag()*deposit;
		}
		ph->move ();

#endif

		ph->scatter ();
	    }
	}
	printf("%f %d %f %f %f %f %f %f\n", L, photons, Ipar, Iper, IxyR, IxyI, (Ipar-Iper)/(Ipar+Iper), -2*IxyI/(Ipar+Iper));
    }
    delete ph;
}
