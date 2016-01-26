#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "basis.h"

#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

class Thamiltonian{
public:

    Thamiltonian ( Tbasis bas, double t_hop, std::vector<double> U_int, std::vector<double> eps_i, int n_eig );
    
    double ada( int i, int j, std::vector<int> &A);
    double adadaa(int i,int j,int k, int l, std::vector<int> &A);
    
    int find_vector( std::vector<int> v );

    void is_symmetric();
    
    void diagonalize();

    void write_eigenvalues( std::string fname );
    std::vector<double> get_eigenvalues(){ return eigenvalues;};
	arma::mat get_eigvec(){ return eigvec;};
	arma::vec get_eigval(){ return eigval;};

	std::vector < std::complex<double> > evolve( std::vector< std::complex<double> > vec, double  Dt);


private:
     arma::Mat< double > hamiltonian;
     double t_hop;
     std::vector<double> U_int;
     std::vector<double> eps_i;
     
     std::vector<double> eigenvalues;
     
     Tbasis basis;
     
     arma::vec eigval;
     arma::mat eigvec;
     
     int no_of_eigenvalues;
         
};


#endif 
