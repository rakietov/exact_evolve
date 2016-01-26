#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "basis.h"
#include "ind.h"
#include "hamiltonian.h"


// ----------------------------------------------------------------------------------------
bool comp_ind2(Tind a, Tind b){
	return (a.V<b.V);
}
// ========================================================================================

// ----------------------------------------------------------------------------------------
Thamiltonian::Thamiltonian ( Tbasis bas, double t, std::vector<double> U, std::vector<double> eps, int n_eig  )
    : t_hop(t), U_int(U), eps_i(eps), basis(bas), no_of_eigenvalues( n_eig )
{
    hamiltonian.resize( basis.get_n_vec(), basis.get_n_vec() );

// Calculate kinetic terms
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
/*
		std::cout<<"For vector: "<<std::endl;
		for(int z = 0; z < basis.get_n_sites(); ++z)
		{
			std::cout<<basis.get_ith_basis_vector(i)[z]<<" ";
		}
			std::cout<<std::endl;
*/
        std::vector<int> v1, v2;
		for(int j = 0; j < basis.get_n_sites() - 1 ; ++j)
		{
			int u = j + 1;

// If M == 4 commented (why?)
			if(j == basis.get_n_sites() - 1)
			{
			    u = 0;
			}
//			std::cout<<"u is: "<<u<<std::endl;


			v1 = basis.get_ith_basis_vector(i);
			v2 = basis.get_ith_basis_vector(i);

       

			double tmp1 = ada(j , u , v1);
			double tmp2 = ada(u , j , v2);

/*
     		for(int z = 0; z < basis.get_n_sites(); ++z)
		    {
			std::cout<<v1[z]<<" ";
		    }   
		    std::cout<<std::endl;
		    for(int z = 0; z < basis.get_n_sites(); ++z)
		    {
			std::cout<<v2[z]<<" ";
		    }          
            std::cout << std::endl;
*/

			for(int z = 0; z < basis.get_n_sites(); ++z)
			{
				if(v1[z] > basis.get_max_n_part() )tmp1 = 0.;
				if(v2[z] > basis.get_max_n_part() )tmp2 = 0.;
			}
			
            int no1 = find_vector( v1 );
            int no2 = find_vector( v2 );

            hamiltonian( no1, i) += tmp1 * t_hop;
            hamiltonian( no2, i) += tmp2 * t_hop;
        }
    }
    
    
// Calculate interaction terms

//    std::cout << "calculate interaction terms" <<std::endl;
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
	   	std::vector<int> v1;

        for(int j = 0; j < basis.get_n_sites(); ++j)
		{
		    v1 = basis.get_ith_basis_vector( i );
/*
                    for(int z = 0; z < basis.get_n_sites(); ++z)
                    {
                        std::cout<<v1[z]<<" ";
                    } std::cout << std::endl;
*/



		    double tmp1 = adadaa(j, j,j,j , v1);
/*
	            for(int z = 0; z < basis.get_n_sites(); ++z)
		    {
			std::cout<<v1[z]<<" ";
	  	    } std::cout << std::endl;
*/
		    for(int z = 0; z < basis.get_n_sites(); ++z)
		    {
			    if(v1[z] > basis.get_max_n_part() )tmp1 = 0.;
		    }

   			int no1 = find_vector( v1 );		   			
   			hamiltonian( no1, i) += 0.5 * tmp1 * U_int[j];
			//	std::cout << no1 <<" "<<i<<" "<< 0.5 * tmp1 * U_int[j] <<std::endl;
        }
    }
    
    
// Calculate onsite energy terms
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
	   	std::vector<int> v1;
        v1 = basis.get_ith_basis_vector(i);
        for(int j = 0; j < basis.get_n_sites(); ++j)
		{
		    double tmp1 = ada(j, j, v1);


		    for(int z = 0; z < basis.get_n_sites(); ++z)
		    {
			    if(v1[z] > basis.get_max_n_part() )tmp1 = 0.;
		    }

     		int no1 = find_vector( v1 ); 			
    		hamiltonian( no1 , i) += tmp1 * eps_i[j];
        }
    }    
    std::cout << "Hamiltonian set up" << std::endl;
 //   hamiltonian.print();
/*
  hamiltonian = hamiltonian * hamiltonian;

for(int i = 0; i < basis.get_n_vec(); ++i){
	hamiltonian(i,i) -= 1000.;
}
*/
   // hamiltonian.print();    
    
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Thamiltonian::diagonalize()
{

    std::cout << "Starting diagonalization" << std::endl;    

  arma::eig_sym(eigval, eigvec, hamiltonian);
//	arma::eigs_sym(eigval, eigvec, hamiltonian, no_of_eigenvalues, "sa"); 
    for( int i = 0; i < eigval.size(); ++i)
    {
        eigenvalues.push_back(eigval[i]);
    }

    std::cout << "Hamiltonian diagonalized" << std::endl;        
    
/*    
    std::cout.precision(11);
    std::cout.setf(std::ios::fixed);

    //std::cout << eigval[3];

    eigval.raw_print( std::cout);
*/
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Thamiltonian::write_eigenvalues( std::string fname )
{
    std::fstream fs;
    
    fname += std::string( "_eigenvalues.txt");
    fs.open( fname.c_str(), std::fstream::app);

    fs.precision(15);
    for( int i = 0; i < eigenvalues.size(); ++i)
    {
        fs<<i <<" "<<eigenvalues[i]<<std::endl;
    }
    fs.close();
    
/*    
    std::cout.precision(11);
    std::cout.setf(std::ios::fixed);

    //std::cout << eigval[3];

    eigval.raw_print( std::cout);
*/
}
// ========================================================================================



// ----------------------------------------------------------------------------------------
std::vector < std::complex<double> > Thamiltonian::evolve( std::vector< std::complex<double> > vec, double  Dt)
{
	arma::Mat< std::complex<double> > alpha, evol;
	alpha.resize( basis.get_n_vec(), basis.get_n_vec() );
	evol.resize( basis.get_n_vec(), basis.get_n_vec() );
	
	//eigvec.print();
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		for(int j = 0; j < basis.get_n_vec(); ++j)
		{
			//std::cout << eigvec(i,j) << std::endl;
			alpha(i,j) = std::complex<double> ( eigvec(i,j) , 0.);
		}
	
	}
	
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		evol(i) = std::exp(  std::complex<double> ( 0.,-eigval[i]*Dt ) );    
	}

	arma::Mat< std::complex<double> > res_U = alpha*evol*alpha.t();

	//res_U.print();

	std::vector < std::complex<double> > res_vec;
	
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		res_vec.push_back( std::complex<double> ( 0.,0.) );
	} 
	
	for(int i = 0; i < basis.get_n_vec(); ++i)
	{
		for(int j = 0; j < basis.get_n_vec(); ++j)
		{
			res_vec[i] += vec[j]*res_U(j,i);
		}
	}     


	return res_vec;
}
// ========================================================================================




// ----------------------------------------------------------------------------------------
int Thamiltonian::find_vector( std::vector<int> v )
{
			Tind I1;
			I1.tag(v, 0 );
			
		    std::vector< Tind > hash_Tab = basis.get_hash_table();
  			std::vector< Tind >::iterator p1 = std::lower_bound ( hash_Tab.begin(), hash_Tab.end() - 1, I1, comp_ind2);

   			return hash_Tab[ p1 - hash_Tab.begin() ].pos;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
double Thamiltonian::ada( int i, int j, std::vector<int> &A)
{
	double tmp1 = A[i];	
	double tmp2 = A[j];


	A[i]++;

	if(A[j] > 0)A[j]--;
	else return 0.;

	if( i == j) return A[i];
	else return std::sqrt((tmp1+1.)*(tmp2));
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Thamiltonian::is_symmetric( )
{
   for(int i = 0; i < basis.get_n_vec(); ++i)
   {
       for(int j = 0; j < i; ++j)
       {
	   if(hamiltonian(i,j) != hamiltonian(j,i))
	   {
 	       std::cout<<"Hamiltonian not symmetric at "<<i<<" "<<j<<std::endl;
	       return;
	   }
       }

   } 



}
// ========================================================================================




// ----------------------------------------------------------------------------------------
double Thamiltonian::adadaa(int i,int j,int k, int l, std::vector<int> &A)
{
	double tmp1 = A[l];
	if(A[l] > 0)A[l]--;
	else return 0;
	
	double tmp2 = A[k];
	if(A[k] > 0)A[k]--;
	else return 0;

	double tmp3 = A[j];
	A[j]++;
	double tmp4 = A[i];
	A[i]++;

	return sqrt(tmp1*tmp2*(tmp3+1)*(tmp4+1));
}
// ========================================================================================












