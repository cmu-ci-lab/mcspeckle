#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include "dmiev.h"

typedef std::complex<double> dcmplx;
#define pi M_PI

class scatterer {
public:
    scatterer(double x, dcmplx m, int N=10001, int nslot=10000) : x(x), m(m), N(N), nslot(nslot), invtab(new double[nslot+1]) {
	double *mulist = new double[N];
	s1tab = new dcmplx[N];
	s2tab = new dcmplx[N];
	for (int i=0; i<N; i++)
	    mulist[i] = cos(pi*i/(N-1.));
	miev(x, m, N, mulist, s1tab, s2tab);
	delete[] mulist;

	psca0 = (norm(s1tab[0]) + norm(s2tab[0]))/2;
	mievinfo(x, m, &qext, &qsca, &g);
	geninvtab();

	fprintf(stderr, "xx, m=%f %f %f\n", x, m.real(), m.imag());
	fprintf(stderr, "qext, qsca, g=%f %f %f\n", qext, qsca, g);
	fprintf(stderr, "psca0, csca=%g %g\n", psca0, csca());
    }

    ~scatterer() {
	delete[] s1tab;
	delete[] s2tab;
	delete[] invtab;
    }

    double albedo() const { return  qsca/qext; }

    /* =int ( |s2|^2*cos(phi)^2 + |s1|^2*sin(phi)^2 ) sin(theta) dtheta dphi */
    double csca() const { return pi*qsca*x*x; }

    void phasef(int n, double* mu, dcmplx* s1, dcmplx* s2) const {
	miev(x, m, n, mu, s1, s2);
    }

    void phasef(double mu, dcmplx* s1, dcmplx* s2) const {
	amiev(x, m, mu, s1, s2);
    }

    /* the version that uses the stored look up table  */
    void phasef_lu(double mu, dcmplx* s1, dcmplx* s2) const {
	double theta = acos(mu);
	int p = (int)floor(theta/pi*(N-1));
	*s1 = s1tab[p];
	*s2 = s2tab[p];
    }

    /* yield the scattering angle given the probability within (0,1) */
    double draw_mu(double p) const {
	int P = (int)floor(p*nslot);
	if ( p >= 1 ) return invtab[nslot];
	else {
	    assert (P >= 0 && P < nslot);
	    return invtab[P] + nslot*(p-P/(0.+nslot))*(invtab[P+1] - invtab[P]);
	}
    }

    double x, qext, qsca, g, psca0; /* psca0=(|S2|^2+|S1|^2)/2 with theta=0 */
    dcmplx m;

private:
    int N;			/* the number of angles in the forward lookup table */
    int nslot;			/* the precision of mu is 2./nslot */
    dcmplx *s1tab, *s2tab;      /* the forward lookup table storing s1 and s2 */
    double *invtab;		/* the inverse lookup table */

    void geninvtab() {
	int N = 360001;		/* use a fixed number of angles */
	double *tab = new double[N];
	double *mulist = new double[N];
	dcmplx *s1list = new dcmplx[N];
	dcmplx *s2list = new dcmplx[N];
	double f = pi/(N-1.)/(qsca*x*x);
	double sq1, sq2, p;
	int j;

	for (int i=0; i<N; i++)
	    /* the equal angle slicing gives much better accuracy as it samples more densely in the forward and backward directions */
	    mulist[i] = cos(pi*(1 - i/(N-1.)));

	miev(x, m, N, mulist, s1list, s2list);

	/* generate the forward table */
	tab[0] = 0;
	sq1 = norm(s2list[0]) + norm(s1list[0]);
	for (int i=1; i<N; i++) {
	    sq2 = norm(s2list[i]) + norm(s1list[i]);
	    tab[i] = tab[i-1] + f*(sq1 + sq2)/2.*sin(pi*(1 - (i-0.5)/(N-1)));
	    sq1 = sq2;
	}

	fprintf(stderr, "The last entry of the forward table (should be 1): %12.10g \n", tab[N-1]);

	/* generate the inverse table */
	invtab[0] = -1;
	invtab[nslot] = 1;
	j = 0;
	for (int i=1; i<nslot; i++) {
	    p = i/(0. + nslot);
	    if (tab[j] > p) {
		fprintf(stderr, "Warning! Please increase the number N of sampling angles or decrease nslot.\nj=%d p=%f tab[j-1:j+1]=%f %f %f\nDo not meet the required precision: %f\n", j, p, tab[j-1], tab[j], tab[j+1], 1./nslot);
		exit(1);
	    }
	    else
		while (1) {
		    j ++;
		    if (tab[j] > p) break;
		}
	    /* find p to be within tab[j-1] and tab[j] */
	    invtab[i] = mulist[j-1] + (mulist[j] - mulist[j-1])
				     *(p-tab[j-1])/(tab[j] - tab[j-1]);
	}

	/*
	for (int i=0; i<nslot+1; i++) 
	    fprintf(stderr, "i, invtab[i]: %d %12.10g\n", i, invtab[i]); 
	*/

	delete[] tab;
	delete[] mulist;
	delete[] s1list;
	delete[] s2list;
    }
};


