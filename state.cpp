#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "basis.h"
#include "hamiltonian.h"
#include "state.h"



// ----------------------------------------------------------------------------------------
Tstate::Tstate ( Tbasis bas,  Thamiltonian h, std::vector<std::complex<double> > coe)
	: basis( bas ), coeffs (coe), hamiltonian( h )
{
  

}
// =======================================================================================

 
// ----------------------------------------------------------------------------------------
double Tstate::calc_norm ()
{
	double norm_sq = 0.0;
	for(int i = 0; i < coeffs.size(); ++i)
	{
		norm_sq += std::pow(std::abs( coeffs[i] ),2.);
	}
	
	return std::sqrt(norm_sq);
}
// =======================================================================================

// ----------------------------------------------------------------------------------------
void Tstate::print()
{
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		auto v = basis.get_ith_basis_vector( i );		
		for(int j = 0; j < basis.get_n_sites(); ++j)
		{
			std::cout << v[j] << " ";
		}
		std::cout << " : " << coeffs[i] << std::endl;
	}

}
// =======================================================================================



// ----------------------------------------------------------------------------------------
void Tstate::evolve(double  Dt)
{
	arma::Mat< std::complex<double> > alpha, evol;
	arma::mat eigve;
	arma::vec eigva;
	alpha.resize( basis.get_n_vec(), basis.get_n_vec() );
	evol.resize( basis.get_n_vec(), basis.get_n_vec() );
	
	eigve = hamiltonian.get_eigvec();
	eigva = hamiltonian.get_eigval();
	//eigvec.print();
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		for(int j = 0; j < basis.get_n_vec(); ++j)
		{
			//std::cout << eigvec(i,j) << std::endl;
			alpha(i,j) = std::complex<double> ( eigve(i,j) , 0.);
		}
	
	}
	
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		evol(i, i) = std::exp(  std::complex<double> ( 0.,-eigva[i]*Dt ) );    
	}

	arma::Mat< std::complex<double> > res_U = alpha*evol*alpha.t();

	//auto asdf = res_U * res_U.t();
	//asdf.print();
	//res_U.print();

	std::vector < std::complex<double> > res_vec;
	
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		res_vec.push_back( std::complex<double> ( 0.,0.) );
	} 
	
	/*
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		std::cout << coeffs[i]<< std::endl;
	}
	*/

	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		for(int j = 0; j < basis.get_n_vec(); ++j)
		{
			res_vec[i] += coeffs[j]*res_U(j,i);
		}
	}     
	
	for( int i = 0; i < basis.get_n_vec(); ++i) 
	{
		coeffs[i] = res_vec[i];
	}

}
// ========================================================================================



// ----------------------------------------------------------------------------------------
double Tstate::calc_av_ni( int i )
{
	double av_ni = 0.0;
	for(int j = 0; j < coeffs.size(); ++j)
	{
		av_ni += std::pow(std::abs( coeffs[j] ),2.) * basis.get_ith_basis_vector( j )[i];
	}
	
	return av_ni;
}
// =======================================================================================




// ----------------------------------------------------------------------------------------
double Tstate::calc_bosonic_even_odd_imbalance()
{
	double nom = 0.;
	double denom = 0.;

	for(int i = 0; i < basis.get_n_sites(); i+=2)
	{
		nom += this -> calc_av_ni( i );
		denom += this -> calc_av_ni( i ); 
	}
	for(int i = 1; i < basis.get_n_sites(); i+=2)
	{
		nom -= this -> calc_av_ni( i );
		denom += this -> calc_av_ni( i );
	}

	return nom / denom;
}
// =======================================================================================



