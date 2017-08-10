# use Intel Fortran
#FC = ifort -r8 -g
#FLIB = -L/usr/local/intel_parallel_studio_xe_2017_update2/compilers_and_libraries_2017.2.174/linux/compiler/lib/intel64_lin/ 
# use g77 and fort77 
#FC = fort77 -r8 -g
#FLIB = -lg2c
FC = gfortran -fdefault-real-8 -fdefault-double-8 -g
FLIB = -lgfortran
CC = g++ -O3 -g -Wall -Wunused-variable -I.

mtwist.o: mtwist.c
	${CC} -c mtwist.c

deplength: deplength.cpp MIEV0noP.o ErrPack.o pol_montecarlo.h mtwist.o
	${CC} -o deplength deplength.cpp mtwist.o MIEV0noP.o ErrPack.o -L/usr/lib ${FLIB}

speckle: mtwist.o MIEV0noP.o ErrPack.o speckle.cpp pol_montecarlo.h
	${CC} -o speckle speckle.cpp mtwist.o MIEV0noP.o ErrPack.o -L/usr/lib -lnetcdf_c++ -lnetcdf ${FLIB}

rmuller: mtwist.o MIEV0noP.o ErrPack.o rmuller.cpp pol_montecarlo.h
	${CC} -o rmuller rmuller.cpp mtwist.o MIEV0noP.o ErrPack.o -L/usr/lib -lnetcdf_c++ -lnetcdf ${FLIB}

##################
## dependencies ##
##################
MIEV0noP.o: MIEV0noP.f
	${FC} -c MIEV0noP.f

ErrPack.o: ErrPack.f
	${FC} -c ErrPack.f

pol_montecarlo.h: scatterer.h dmiev.h
