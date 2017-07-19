====================================================
Release 0.5 of Electric Field Monte Carlo code (EMC)

Dr. Min Xu <minxu@sci.ccny.cuny.edu>
see COPYRIGHT
====================================================

1. Pre-requirements
EMC depends on the installation of certain softwares.

(a) a C++ compiler
   GNU C++ should be fine. The compiler should support templates.
 
(b) MIEV package by Dr. Warren J. Wiscombe
(wiscombe@climate.gsfc.nasa.gov).
   This is an excellent code for Mie calculations.

(c) mtwist-0.5 package by Prof. Geoff Kuenning at Harvey Mudd College
(geoff@cs.hmc.edu). 
   This is an Mersenne Twist Pseudorandom Number Generator
   Package. The Mersenne Twist method for generating pseudorandom 
   numbers is an extremely fast, flexible, and desirable approach to
   random-number generation. It has superb statistical properties and
   a ridiculously long period (2^19937-1).

(d) Netcdf package (optional)
   Netcdf API is used in some examples coming with EMC to store
   simulation outputs. The benefit of using netcdf files is the
   platform-independence of these files.

I have included parts of (a-c) in the distribution only for the purpose 
your convenience. Full credit should be given to respective authors.

2. EMC
EMC package is implemented in C++ and contained in three headers
files. We use x as the size parameter of the particle, m the complex
relative refractive index, S1 and S2 are the diagonal elements of the
amplitude scattering matrix (whose off-diagonal elements are zero) in 
convention of Hulst.

dmiev.h:	a C interface for MIEV package, used by scatterer.h

scatterer.h:	a C++ header implementing a class "scatterer". The key
		methods include: 
		 
		 // N and nslot can be increased for a higher accuracy
		 // of the forward and inverse table.
		 scatterer(double x, dcmplx m, int N=10001, int nslot=10000);
		 
		 // obtain S1 and S2 at a list of cosines of angles of
		 // length n
		 phasef(int n, double* mu, dcmplx* s1, dcmplx* s2);
		 
		 // obtain S1 and S2 at the cosine of one angle mu
		 phasef(double mu, dcmplx* s1, dcmplx* s2);

		 // a quicker version using the lookup table pre-computed
		 phasef_lu(double mu, dcmplx* s1, dcmplx* s2);

		 // yield the scattering angle given the probability
		 // within (0,1) using the inverse table pre-computed
		 draw_mu(double p);


pol_montecarlo.h:	the core header file of EMC implementing a
		class "photonPacket". The key methods include:

		 // initialization
		 photonPacket(const scatterer* sct, char* fname=NULL,
		 unsigned long seed=0); 
		 
		 // launch the photons. The incident electric field is
		 // given by (E1, E2) where light is propagating in the
		 // direction (u, v, w), E1 is in the direction (l, m,
		 // n), and E2 is in the direction specified by the cross 
		 // product of the above two directions. Light is
		 // incident upon the position (x,y,z) at time t. The
		 // incident light intensity must be unity (|E1|^2 +
		 // |E2|^2=1) and the directions must be unit vectors.
		 void launch(dcmplx E1=1, dcmplx E2=0, double l=1,
		 double m=0, double n=0, double u=0, double v=0,
		 double w=1, double x=0, double y=0, double z=0,
		 double t=0); 

		 // move to next scattering or absorption event
		 void move();

		 // absorption event
		 void absorb();

		 // scattering event
		 void scatter();

		 // next event estimator of (Ed1, Ed2) with the
                 // direction vector Q at time td crossing the boundary
                 // z=zd in the direction (ud, vd, wd) 
		 void pointestimator(double *td, double *deposit,
                 double Q[3][3], dcmplx *Ed1, dcmplx *Ed2, double zd,
                 double ud=0, double vd=0, double wd=1); 

                
		The photon status is stored in the public data members
                of the class, they are:

		 xold, yold, zold, told:
		 the position and time of the previous scattering
		 event

		 x, y, z, t, weight:
		 the current position and time of photon and its weight

		 E1, E2:
		 electric field components

	         P[3][3]:
		 (l, m, n) the unit vector for E1 direction
		 (p, q, r) the unit vector for E2 direction
		 (u, v, w) the propagation direction

		 nsct:                   
		 number of scattering events
