#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>

#include "ind.h"

#ifndef BASIS_H_
#define BASIS_H_

class Tbasis{
public:

    Tbasis(){ n_sites = 0; n_part = 0; n_vec = 0; };
    Tbasis( int n_s, int n_p );
    ~Tbasis(){ /*std::cout << __func__ << std::endl; */};
    
    int comp_n_vec();
    
    int get_n_sites(){ return n_sites; };
    int get_n_part(){ return n_part; };
    int get_max_n_part(){ return max_n_part; };
    int get_n_vec(){ return n_vec; };
    
    std::vector<int> get_ith_basis_vector( int i ){return basis_vectors[i]; };
    std::vector< Tind > get_hash_table(){ return hash_table; };
    
    void write_vectors( std::string fname );
    void fill_hash_table();
    void write_hash_table();

private:
    int n_sites;
    int n_part;
    int max_n_part;
    int n_vec;
    
	std::vector< std::vector<int> > basis_vectors;
	std::vector< Tind > hash_table;
};


#endif 
