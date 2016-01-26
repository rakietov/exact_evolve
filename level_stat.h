#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "hamiltonian.h"

#ifndef LEVEL_STAT_H_
#define LEVEL_STAT_H_

class Tlevel_stat{
public:

    Tlevel_stat ( std::vector<double> eigv, int first, int last, double del);
    void write_eigenvalues();
    void write_unfolded_spectrum();
    void write_unfolded_level_spacings( std::string fname );
    void write_avRBSL( std::string fname );



    double sigma_0( double E );
    double sigma_d( double E, double delta );
    double mean_sq_deviation( double delta, double Emin, double Emax, double n_points);

    void unfold_spectrum();
    
    
    std::vector<double> get_unfolded_level_spacings(); 
    

private:
     std::vector<double> eigenvalues;
     int first_eig, last_eig;
     double delta;
     std::vector<double> unfolded_spectrum;

};


#endif 
