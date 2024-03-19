# Global Variables controlling the compilation 
# h5fc is the fortran compiler with implied hdf5 linking
FC = h5fc

# With OpenMP
FFLAGS = -march=native -Wall -Ofast -fPIC -g -fopenmp

# Without OpenMP
#FFLAGS = -march=native -Wall -Ofast -fPIC -g

F2PYFL = --f90flags=-fopenmp -lgomp

# This target tests the double precision linear algebra module
latest: dpla.f90 latest.f90 Makefile
	$(FC) $(FFLAGS) dpla.f90 latest.f90 -o latest

# Tests the operators constructed in ops.f90
opstest: ops.f90 triang.f90 matrices.f90 opstest.f90 quad.f90 param.f90 dpla.f90 Makefile 
	$(FC) $(FFLAGS) dpla.f90 matrices.f90 param.f90 triang.f90  ops.f90 opstest.f90 quad.f90 -o opstest

# Tests the construction of the G-Matrix
gtest:  dpla.f90 quad.f90 matrices.f90 ops.f90 Gtest.f90 triang.f90 Makefile
	$(FC) $(FFLAGS) dpla.f90 quad.f90 triang.f90 matrices.f90 ops.f90 Gtest.f90 -o gtest

# A accuracy/convergence test. Prints the discrepance from the reference solution.
dsdgacctest: dpla.f90 testcases.f90 triang.f90 euler.f90 sdgsol.f90 matrices.f90 dsdgtest.f90 param.f90 Makefile
	$(FC) $(FFLAGS) dpla.f90 matrices.f90 param.f90 triang.f90 euler.f90 quad.f90 ops.f90 testcases.f90 sdgsol.f90 acctest.f90 -o acctest

# Tests the vanilla DG solver on a simple two-dimensional Riemann problem with circular or translation invariance
dsdgtest: testcases.f90 triang.f90 euler.f90 sdgsol.f90 matrices.f90 dsdgtest.f90 param.f90 dpla.f90 Makefile
	$(FC) $(FFLAGS) dpla.f90 matrices.f90 param.f90 triang.f90 euler.f90 quad.f90 ops.f90 testcases.f90 sdgsol.f90 dsdgtest.f90 -o dsdgtest

# A unstructured stabilizes solver on triangles - used for the forward facing step and airfoil tests
udsdg: testcases.f90 triang.f90 euler.f90 sdgsol.f90 ops.f90 matrices.f90 udsdg.f90 param.f90 dpla.f90 Makefile
	$(FC) $(FFLAGS) dpla.f90 matrices.f90 param.f90 triang.f90 euler.f90 quad.f90 ops.f90 testcases.f90 sdgsol.f90 udsdg.f90 -o udsdg


# A dummy target to carry out all tests for a given order.
alltests: udsdg dsdgtest dsdgacctest gtest opstest latest
	./latest
	./opstest
	./gtest
	#ulimit -s unlimited		# We need more stack space
	./acctest grids/acc1.1.node grids/acc1.1.ele acc1 1.0
	./acctest grids/acc2.1.node grids/acc2.1.ele acc2 1.0
	ulimit -s unlimited & ./acctest grids/acc3.1.node grids/acc3.1.ele acc3 1.0
	./udsdg grids/Naca0012.1.node grids/Naca0012.1.ele Naca0012 3.0
	./udsdg grids/NacaLR.1.node grids/NacaLR.1.ele Naca0012LR 3.0
	./dsdgtest grids/empty.1.node grids/empty.1.ele sedov 0.2

allfigures:
	python3 naca0012.py
	python3 euler.py 
	python3 ffs.py 
	python3 sedov.py

clean:
	rm *.o *.mod *.out vgcore.* acctest dgtest gtest dsdgtest cdtest dualgfv dualgtest latest llftest ddgtest specdiag sdgtest opstest udddg uddg udg udsdg ullf usok usw
