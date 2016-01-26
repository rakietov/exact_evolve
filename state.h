#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "basis.h"

#ifndef STATE_H_
#define STATE_H_

class Tstate{
public:

    Tstate ( Tbasis bas, Thamiltonian ham, std::vector< std::complex<double> > coeff );

	void evolve( double Dt );   
 
	double calc_norm();
	double calc_av_ni( int i );	
	double calc_bosonic_even_odd_imbalance();

	void print();



private:
	Tbasis basis;	
	Thamiltonian hamiltonian;
	std::vector< std::complex<double> > coeffs;

      
};

#endif 
