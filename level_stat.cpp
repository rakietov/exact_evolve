#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "level_stat.h"
#include "hamiltonian.h"


// ----------------------------------------------------------------------------------------
Tlevel_stat::Tlevel_stat (  std::vector<double> eigv, int first, int last, double del )
    : eigenvalues( eigv ), first_eig( first ), last_eig( last ), delta( del )
{


}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tlevel_stat::write_eigenvalues( )
{
    for( int i = 0; i < eigenvalues.size(); ++i)
    {
        std:: cout<<i<<" "<<eigenvalues[i]<<std::endl;
    }

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tlevel_stat::write_unfolded_spectrum( )
{
    for( int i = 0; i < unfolded_spectrum.size(); ++i )
    {
        std:: cout<<i<<" "<<unfolded_spectrum[i]<<std::endl;
    }

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
std::vector<double> Tlevel_stat::get_unfolded_level_spacings()
{
    std::vector<double> spacings;
    
    for( int i = 0; i < unfolded_spectrum.size()-1; ++i)
    {
        spacings.push_back( unfolded_spectrum[i+1] - unfolded_spectrum[i] );
    }

    return spacings;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tlevel_stat::write_unfolded_level_spacings( std::string fname )
{
    std::fstream fs;
    
    fname += std::string( "_spacings.txt");
    fs.open( fname.c_str(), std::fstream::app);
    
    for( int i = 0; i < unfolded_spectrum.size()-1; ++i)
    {
        fs << unfolded_spectrum[i+1] - unfolded_spectrum[i] << std::endl;
    }

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
double Tlevel_stat::sigma_0( double E)
{
    std::vector< double >::iterator p1 = std::lower_bound ( eigenvalues.begin() + first_eig, eigenvalues.begin() + last_eig - 1, E);
    
    return (p1 - eigenvalues.begin())*(1./(last_eig - first_eig)) ;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
double Tlevel_stat::sigma_d( double E, double delta )
{
    double sum = 0.;
    for( int i = first_eig; i < last_eig; ++i)
    {
        sum += 0.5 + 0.5* std::erf( (E - eigenvalues[i])/delta);
    }
    return sum*(1./(last_eig-first_eig) ) ;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
double Tlevel_stat::mean_sq_deviation( double delta, double Emin, double Emax, double n_points)
{
    double sum = 0.;
    double dE = (Emax-Emin)/n_points ;
    for( int i = 0; i < n_points; ++i)
    {
        double Eact = Emin +  dE * i;
        sum += dE * ( sigma_d( Eact, delta) - sigma_0( Eact) )*( sigma_d( Eact, delta) - sigma_0( Eact) );
    }
    return sum;
}
// ========================================================================================



// ----------------------------------------------------------------------------------------
void Tlevel_stat::unfold_spectrum()
{


    double min_MSD = 1000.;
    double res_delta;
/*
    for( int id = 1; id < 30; ++id)
    {
        double tmp_d = id*0.1;
        double tmp_MSD = mean_sq_deviation( tmp_d, eigenvalues[id] - 2., eigenvalues[id] + 2., 10000);
        std::cout << __func__ <<" "<<tmp_d << " " << tmp_MSD << std::endl;
        if(tmp_MSD < min_MSD)
        {
            min_MSD = tmp_MSD;
            res_delta = tmp_d;
        }
    }


    //std::cout << __func__ <<min_MSD << " " << res_delta << std::endl;
*/
 //   res_delta = 36.3;
    //delta = res_delta;
    
    for(int i = first_eig; i < last_eig; ++i)
    {
        unfolded_spectrum.push_back( (last_eig-first_eig)*sigma_d( eigenvalues[i], delta) );
        //unfolded_spectrum.push_back(eigenvalues[i] );
    }
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tlevel_stat::write_avRBSL( std::string fname )
{
    double avR = 0.;
    
    for( int i = first_eig; i < last_eig; ++i)
    {
        double nom, denom;
        
        nom = eigenvalues[i-1] - eigenvalues[i-2];
        denom = eigenvalues[i-1] - eigenvalues[i-2];    

        if( eigenvalues[i] - eigenvalues[i-1] < nom)
        {
            nom = eigenvalues[i] - eigenvalues[i-1];
        }
        if( eigenvalues[i] - eigenvalues[i-1] > denom)        
        {
            denom = eigenvalues[i] - eigenvalues[i-1];
        }
        
        avR +=  nom/denom;
    }    
    avR = avR /(last_eig - first_eig) ;

    
    std::fstream fs;
    
    fname += std::string( "_avRBSL.txt");
    fs.open( fname.c_str(), std::fstream::app);
    
    fs << avR << std::endl;
    
}
// ========================================================================================













