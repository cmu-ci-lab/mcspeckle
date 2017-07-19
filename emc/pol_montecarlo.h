#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <fstream.h>
#include "mtwist.h"
#include "scatterer.h"

typedef std::complex<double> dcmplx;

#define pi M_PI
#define ranf mtrand->drand

/*
  In the Monte Carlo code, the length unit is lsa = (mu_a+mu_s)^(-1)
  the time is also measured in lsa such that the real value of the
  speed of light is irrevelant.  To convert the time in lsa to the
  time in ps, use the following formula:
    [time in ps] = [time in lsa] * (lsa*n/c)
  where n is the refractive index, c=0.2998mm/ps (in vacuum)

  Remember time and length _is_ measured in lsa.
*/
class photonPacket
{
public:
    photonPacket(const scatterer* sct, char* fname=NULL, unsigned long seed=0) : sct(sct), albedo(sct->albedo()), csca(sct->csca()), psca0(sct->psca0) {
	mtrand = new mt_prng();
	if (fname != NULL) {
	    ifstream stateIn(fname);
	    if( stateIn ) {
		stateIn >> *mtrand;
		stateIn.close();
		fprintf(stderr, "Loaded the seed file %s.\n", fname);
	    }
	    else {
		fprintf(stderr, "The seed file %s does not exist.\n", fname);
		if (seed != 0) {
		    fprintf(stderr, "Use the supplied seed %ld.\n", seed);
		    mtrand->seed32new(seed);
		} else {
		    fprintf(stderr, "Use the goodseed.\n");
		    mtrand->goodseed();
		}		    
	    }
	}
	else {
	    if (seed != 0) {
	    	fprintf(stderr, "Use the supplied seed %ld.\n", seed);
	     	mtrand->seed32new(seed);
	    } else {
	    	fprintf(stderr, "Use the goodseed.\n");
	    	mtrand->goodseed();
	    }		    
	}
    }

    void saveRand(char* fname) const {
	// A stream is convenient for saving to a file.
	ofstream stateOut(fname);
	if( stateOut ) {
	    stateOut << *mtrand;
	    stateOut.close();
	    fprintf(stderr, "Save the seed file %s.\n", fname);
	}
	else
	    fprintf(stderr, "Can not open %s for writing.\n", fname);
    }

    void launch(dcmplx E1=1, dcmplx E2=0, double l=1, double m=0, double n=0, double u=0, double v=0, double w=1, double x=0, double y=0, double z=0, double t=0) {
	double p, q, r;
	xold = x; yold = y; zold = z; told = t;
	this->x = x; this->y = y; this->z = z; this->t = t;
	this->E1 = E1; this->E2 = E2;
	this->E1old = E1; this->E2old = E2;
	assert( fabs(norm(E1) + norm(E2) - 1) < 1e-12 );
	p = n*v - m*w; q = l*w - n*u; r = m*u - l*v;
	P[0][0] = l; P[0][1] = m; P[0][2] = n;
	P[1][0] = p; P[1][1] = q; P[1][2] = r;
	P[2][0] = u; P[2][1] = v; P[2][2] = w;
	Pold[0][0] = l; Pold[0][1] = m; Pold[0][2] = n;
	Pold[1][0] = p; Pold[1][1] = q; Pold[1][2] = r;
	Pold[2][0] = u; Pold[2][1] = v; Pold[2][2] = w;
	weight = 1.0;
	nsct = 0;
#if 0
	printf("Enter:\n");
	for (int i=0; i<3; i++)
	    printf("    P:  %d [%f %f %f]\n", i, P[i][0], P[i][1], P[i][2]);
	printf("    E:  %f+(%fi)  %f+(%fi)\n", E1.real(), E1.imag(), E2.real(),E2.imag());
#endif
    }


    void move() {  /* move to next scattering or absorption event */
	double d=0;
	while (d == 0) 
	    d=-log(mtrand->drand());
	xold = x;
	yold = y;
	zold = z;
	told = t;
	x += d * P[2][0];
	y += d * P[2][1];
	z += d * P[2][2];
	t += d;			/* total path length in mfp */
    }


    void absorb() {		/* Absorb light in the medium */
	weight *= albedo;
	if (weight < 0.001){	/* Roulette */
	    bit -= weight;
	    if (ranf() > 0.1) weight = 0; else weight /= 0.1;
	    bit += weight;
	}
    }


    void scatter() {  /* Scatter photon and establish new direction */
	/* get the new scattering mu=cos(theta) and phi */
	dcmplx s1, s2, e;
	double A[3][3];
	double mu, nu, phi, cosphi, sinphi, s1sq, s2sq, e1sq, e2sq, e12, F, a, b, c, Fmax;
	ncall ++;

	/* Record the values of E1 & E2 before scattering */
	E1old = E1;
	E2old = E2;

#if 0	/* More accurate than the second choice. Will slow down the simulation for larger particles tremendously. Required for the situations where the angular distribution is important or critical, for example, in rmuller.cpp. */ 
	/* For x=1, m=1.12, photons=100000, it takes 6.5min */
	/* rejection method for both mu1 and phi */
	e1sq = norm(E1);
	e2sq = norm(E2);
	e12 = (E1*conj(E2)).real();

	/* rejection for mu. Estimated rejection per success is 4*pi*psca0/csca */
	while (1) {
	    mu = ranf()*2 - 1.;
	    /* sct->phasef(mu, &s1, &s2); */
	    /* use a shortcut */
	    sct->phasef_lu(mu, &s1, &s2);
	    s2sq =  norm(s2);	/* norm yields |s2|^2 */
	    s1sq =  norm(s1);
	    rejmu ++;
	    assert( (s1sq + s2sq)/2./psca0 < 1 + 1e-8 );
	    if ( ranf()*psca0 < (s1sq + s2sq)/2 ) break;
	}
	// printf("Rejection mu: %f\n", rejmu/ncall);

	a = s2sq*e1sq + s1sq*e2sq;
	b = s1sq*e1sq + s2sq*e2sq;
	c = 2*(s2sq-s1sq)*e12;
	Fmax = (a+b)/2 + sqrt( (a-b)*(a-b) + c*c )/2;

	/* rejection for phi */
	while (1) {
	    phi = ranf()*2*pi;
	    cosphi = cos(phi);
	    sinphi = sin(phi);
	    F = a*cosphi*cosphi + b*sinphi*sinphi + c*cosphi*sinphi;
	    rejphi ++;
	    assert( F/Fmax < 1 + 1e-8 );
	    if ( ranf()*Fmax < F ) break;
	}
	// printf("          phi: %f\n", rejphi/ncall);

#else   /* Less accurate in angular sampling than the first method. Much faster than the first choice for larger particles. */
	/* For x=1, m=1.12, photons=100000, it takes 5.2min */
	/* draw mu */
	mu = sct->draw_mu(ranf());
	/* sct->phasef(mu, &s1, &s2); */
	sct->phasef_lu(mu, &s1, &s2);
	s2sq =  norm(s2);	/* norm yields |s2|^2 */
	s1sq =  norm(s1);

	e1sq = norm(E1);
	e2sq = norm(E2);
	e12 = (E1*conj(E2)).real();

	a = s2sq*e1sq + s1sq*e2sq;
	b = s1sq*e1sq + s2sq*e2sq;
	c = 2*(s2sq-s1sq)*e12;
	Fmax = (a+b)/2 + sqrt( (a-b)*(a-b) + c*c )/2;

	/* rejection method for phi */
	while (1) {
	    phi = ranf()*2*pi;
	    cosphi = cos(phi);
	    sinphi = sin(phi);
	    F = a*cosphi*cosphi + b*sinphi*sinphi + c*cosphi*sinphi;
	    assert( F/Fmax < 1 + 1e-8 );
	    if ( ranf()*Fmax < F ) break;
	}
#endif

	nu = sqrt(1 - mu*mu);

	/* rotation for direction */
	A[0][0] = mu*cosphi;
	A[0][1] = mu*sinphi;
	A[0][2] = -nu;
	A[1][0] = -sinphi;
	A[1][1] = cosphi;
	A[1][2] = 0;
	A[2][0] = nu*cosphi;
	A[2][1] = nu*sinphi;
	A[2][2] = mu;

	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++)
		Pold[i][j] = P[i][j];

	/* P = A.Pold */
	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++) {
		P[i][j] = 0;
		for (int k=0; k<3; k++)
		    P[i][j] += A[i][k]*Pold[k][j];
	    }

	/* update E1 and E2 */
	F = sqrt(F);
	e = (cosphi*E1 + sinphi*E2)*s2/F;
	E2 = (-sinphi*E1 + cosphi*E2)*s1/F;
	E1 = e;
	assert( fabs(norm(E1) + norm(E2) - 1) < 1e-10);
	nsct ++;

	/* scattering amplitude matrix elements */
	S1 = s1/F;
	S2 = s2/F;

#if 0
	/* check */
	for (int i=0; i<3; i++)
	    for (int j=i; j<3; j++) {
		double P2 = (P[i][0]*P[j][0] + P[i][1]*P[j][1] + P[i][2]*P[j][2]);
		assert( fabs(P2 - (i==j)) < 1e-8 );
	    }
	printf("    Scattering with (%f, %f):\n", mu, phi);
	printf("         E field:     [%g+(%gi) %g+(%gi)]\n", E1.real(), E1.imag(), E2.real(), E2.imag());
#endif
    }


    // next event estimator: (Ed1, Ed2) with the direction vector Q at time td crossing the boundary z=zd in the direction (ud, vd, wd)
    void pointestimator(double *td, double *deposit, double Q[3][3], dcmplx *Ed1, dcmplx *Ed2, double zd, double ud=0, double vd=0, double wd=1) const {
	double mu, nu, sinphi, cosphi, p[3], A[3][3], s1sq, s2sq, e1sq, e2sq, e12, F;
	dcmplx s1, s2;

	assert (fabs(ud*ud + vd*vd + wd*wd - 1) < 1e-12);
	/* the incident propagating direction is Pold[2] and detected in (ud, vd, wd) */
	mu = Pold[2][0]*ud + Pold[2][1]*vd + Pold[2][2]*wd;
	nu = sqrt(1 - mu*mu);

	if (fabs(1 - mu) < 1e-11) { /* photon is already in the right direction */
	    for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
		    Q[i][j] = Pold[i][j];

	    sct->phasef_lu(1., &s1, &s2);
	    s2sq = norm(s2);
	    s1sq = norm(s1);
	    F = (s1sq + s2sq)/2;

	    *Ed1 = E1old*s2/sqrt(F);
	    *Ed2 = E2old*s1/sqrt(F);
	}
	else if (fabs(1 + mu) < 1e-11) { /* photon is in the opposite direction */
	    // In this case, consider a large sphere of m>1, then 
	    // S matrix = diag(S2, S1) = diag(-1,1)
	    // E1'=-E1, E2'=E2
	    // The requirement E1'*m' + E2'*n' + E1*m + E2*n = 0 should be met
	    // The answer is m'=m and n'=-n
	    for (int j=0; j<3; j++) {
		Q[0][j] = Pold[0][j];
		Q[1][j] = -Pold[1][j];
		Q[2][j] = -Pold[2][j];
	    }

	    sct->phasef_lu(-1., &s1, &s2);
	    s2sq = norm(s2);
	    s1sq = norm(s1);
	    F = (s1sq + s2sq)/2;

	    *Ed1 = E1old*s2/sqrt(F);
	    *Ed2 = E2old*s1/sqrt(F);
	}
	else {
	    /* p is the unit vector (u,v,w) crossproduct (ud, vd, wd) */
	    p[0] = (Pold[2][1]*wd - Pold[2][2]*vd)/nu;
	    p[1] = (Pold[2][2]*ud - Pold[2][0]*wd)/nu;
	    p[2] = (Pold[2][0]*vd - Pold[2][1]*ud)/nu;

	    cosphi = Pold[1][0]*p[0] + Pold[1][1]*p[1] + Pold[1][2]*p[2];
	    sinphi = -(Pold[0][0]*p[0] + Pold[0][1]*p[1] + Pold[0][2]*p[2]);

	    /* rotation for direction */
	    A[0][0] = mu*cosphi;
	    A[0][1] = mu*sinphi;
	    A[0][2] = -nu;
	    A[1][0] = -sinphi;
	    A[1][1] = cosphi;
	    A[1][2] = 0;
	    A[2][0] = nu*cosphi;
	    A[2][1] = nu*sinphi;
	    A[2][2] = mu;

	    /* Q = A.Pold */
	    for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
		    Q[i][j] = 0;
		    for (int k=0; k<3; k++)
			Q[i][j] += A[i][k]*Pold[k][j];
	    }

	    /* compute Ed */
	    sct->phasef_lu(mu, &s1, &s2);
	    s2sq = norm(s2);	/* norm yields |s2|^2 */
	    s1sq = norm(s1);

	    e1sq = norm(E1old);
	    e2sq = norm(E2old);
	    e12 = (E1old*conj(E2old)).real();

	    F = (s2sq*e1sq + s1sq*e2sq)*cosphi*cosphi + (s1sq*e1sq + s2sq*e2sq)*sinphi*sinphi + 2*(s2sq-s1sq)*e12*cosphi*sinphi;
	    *Ed1 = (cosphi*E1old + sinphi*E2old)*s2/sqrt(F);
	    *Ed2 = (-sinphi*E1old + cosphi*E2old)*s1/sqrt(F);

	    assert( fabs(norm(*Ed1) + norm(*Ed2) - 1) < 1e-10);
	}

	assert ( fabs(Q[2][0] - ud) < 1e-3 );
	assert ( fabs(Q[2][1] - vd) < 1e-3 );
	assert ( fabs(Q[2][2] - wd) < 1e-3 );

	/* compute deposit and detecting time. detector is located very far from the cell */
	if ( nsct == 0 )
	    if ( fabs(1 - mu) < 1e-11 )           /* ballistic light */
		*deposit = weight*exp(-fabs((z-zd)/wd));
	    else
		*deposit = 0;
	else			/* the collecting solid angle is 1. The introdution of this collection angle DOES ONLY make a difference in the exact forward scattering direction. */
	    *deposit = weight*F/csca*exp(-fabs((z-zd)/wd));  /* F/csca is the probability density to scatter into the outgoing direction */
	*td = t + fabs((z-zd)/wd);

#if 0
	printf("    Pold: [[%f %f %f]; [%f %f %f]; [%f %f %f]]\n", Pold[0][0], Pold[0][1], Pold[0][2], Pold[1][0], Pold[1][1], Pold[1][2], Pold[2][0], Pold[2][1], Pold[2][2]);
	printf("    Eold: [%f+(%fi) %f+(%fi)]\n", E1old.real(), E1old.imag(), E2old.real(), E2old.imag());
	printf("    Scattering angles (mu, cosphi, sinphi): %f %f %f\n", mu, cosphi, sinphi);
	printf("    Rotation matrix: [[%f %f %f]; [%f %f %f]; [%f %f %f]]\n", A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]);
	printf("    Q: [[%f %f %f]; [%f %f %f]; [%f %f %f]]\n", Q[0][0], Q[0][1], Q[0][2], Q[1][0], Q[1][1], Q[1][2], Q[2][0], Q[2][1], Q[2][2]);
	printf("   Ed: [%f+(%fi) %f+(%fi)]\n", Ed1->real(), Ed1->imag(), Ed2->real(), Ed2->imag());
	printf("    F: %f\n\n", F/csca);
#endif
    }


    int terminated() const {
	return weight == 0;
    }


public:
    /* P[3][3]:
       (l, m, n) the unit vector m
       (p, q, r) the cross product of m and n
       (u, v, w) the propagation direction n
    */
    /* Pold, E1old and E2old are the values of P, E1, and E2 before the scattering */
    double xold, yold, zold, told, x, y, z, t, weight, P[3][3], Pold[3][3];
    dcmplx E1, E2, E1old, E2old; 
    dcmplx S1, S2;	        /* S1, S2 the Hulst scattering amplitude element used in scattering */
    int nsct;			/* number of scattering events */

private:
    const scatterer* sct;
    double albedo, csca, psca0;	/* psca0 = (|S2|^2 + |S1|^2)/2 with theta=0 */
    double bit;			/* store the extra bits from roulette */
    double rejmu, rejphi;
    int ncall, N;
    mt_prng *mtrand;
};

#undef ranf
