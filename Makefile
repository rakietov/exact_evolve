bh_exact: bh_exact.cpp Makefile
	 g++ *.cpp -o bh_exact -O3 -std=c++11 -march=native -I ~/armadillo-6.400.3/include -DARMA_DONT_USE_WRAPPER -DARMA_USE_ARPACK -fopenmp -L/opt/intel/Compiler/13.1/mkl/lib/intel64 -lmkl_gnu_thread -lmkl_core -lmkl_intel_lp64 -lmkl_blacs_openmpi_lp64 
 
   
fdiag_dense: fdiag_dense.f90 Makefile
	 gfortran fdiag_dense.f90 -o fdiag_dense -L/opt/intel/Compiler/13.1/mkl/lib/intel64 -lmkl_gnu_thread -lmkl_core -lmkl_intel_ilp64 -lmkl_blacs_openmpi_ilp64 -lpthread -lgomp -larpack -m64

