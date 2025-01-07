# Global Variables controlling the compilation 
# h5fc is the fortran compiler with implied hdf5 linking
FC = h5fc

# With OpenMP
FFLAGS = -march=native -Wall -Ofast -fPIC  -fopenmp #-fdefault-real-16

# Without OpenMP
#FFLAGS = -march=native -Wall -Ofast -fPIC -g

# Tests the vanilla DG solver on a simple two-dimensional Riemann problem with circular or translation invariance
dsdgtest: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 dsdgtest.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 dsdgtest.f90 -o dsdgtest

# An accuracy test. Prints the discrepance from the reference solution.
dsdgacctest: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 acctest.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 acctest.f90 -o acctest

# an unstructured proof of concept solver
udsdg: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 udsdg.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 udsdg.f90 -o udsdg

# Tests the forward facing step problem
ffs: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 ffs.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 ffs.f90 -o ffs

# writes out the collocation points for plotting reasons
triplts: qpla.f90 quadshp.f90 opshp.f90 trianghp.f90 param.f90 csv.f90 bf.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 trianghp.f90 opshp.f90 csv.f90 triplts.f90 -o triplts


# tests the Kelvin-Helmholtz instability
dgkhi: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 dgkhi.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 dgkhi.f90 -o dgkhi

# tests a two-dimensional Riemann problem
twdrp: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 twdrp.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 twdrp.f90 -o twdrp

# tests the double mach reflection problem
dmr: qpla.f90 bf.f90 quadshp.f90 opshp.f90 testcases.f90 trgrid.f90 trianghp.f90 euler.f90 sdgsol.f90 dmr.f90 param.f90 Makefile
	$(FC) $(FFLAGS) param.f90 qpla.f90 bf.f90 quadshp.f90 opshp.f90 trianghp.f90 trgrid.f90 euler.f90 testcases.f90 sdgsol.f90 dmr.f90 -o dmr


