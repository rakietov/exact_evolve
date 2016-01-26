#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include "basis.h"
#include "ind.h"


// ----------------------------------------------------------------------------------------
bool comp_ind(Tind a, Tind b){
	return (a.V<b.V);
}
// ========================================================================================

// ----------------------------------------------------------------------------------------
Tbasis::Tbasis( int n_s, int n_p  ) 
    : n_sites(n_s), n_part(n_p), max_n_part(n_part)
{
    this -> comp_n_vec();
        
    basis_vectors = std::vector<std::vector<int>>  (n_vec + 1,  std::vector<int>(n_sites + 1));

    //std::cout << __func__ << basis_vectors.size() << basis_vectors[0][0] << std::endl;
  
	basis_vectors[0][0] = n_part;

// Generate basis 
	for(int i=0; i < n_vec; ++i)
	{
		int k=0,s=0;
		
		for(int j=0; j < n_sites; ++j)
		{
			bool flag=false;
			if(basis_vectors[i][j] == 0) {
			    k=j;
			}
			for(int l=k;l<n_sites-1;l++) {
			    
			    if(basis_vectors[i][l] != 0) {
			        flag=true;
			    }
			}
			
			k--;
			
			if(!flag)
			{
			    break;
			}
			
			if(j==n_sites-1) {
			    k=j-1;
			}
			
		}

		for(int j=0; j <= k; j++){
			basis_vectors[i+1][j] = basis_vectors[i][j];
			s += basis_vectors[i][j];
		}

		s--;
//		cout<<i+1 << " " << D <<" " <<k+1<<" "<<M<<endl;
        if( k >= 0){
            //std::cout<< __func__ << " " << k << std::endl;
		    basis_vectors[i+1][k] = basis_vectors[i][k]-1;
		    basis_vectors[i+1][k+1] = n_part-s;
        }
	}

	
    int Dact = 0;
    
	for(int i=0;i< n_vec;i++){
		bool flag = true;
		for(int j=0;j< n_sites;j++)
		{
			if(basis_vectors[i][j]>max_n_part)flag = false;
        }
		if(flag){
			for(int j=0;j< n_sites;j++)
				basis_vectors[Dact][j] = basis_vectors[i][j];
			Dact ++;
		}
	}
	n_vec = Dact++;   


	this -> fill_hash_table();
    std:: sort ( hash_table.begin(), hash_table.end(), comp_ind );
	//this -> write_hash_table();

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tbasis::fill_hash_table()
{
	for(int i=0;i< n_vec; ++i )
	{
	    Tind ind1;
	    ind1.tag( basis_vectors[i], i );
		hash_table.push_back( ind1 );
	}

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tbasis::write_hash_table()
{
	for(int i=0;i<n_vec;i++)
		std::cout<<hash_table[i].pos<<" "<<hash_table[i].V<<std::endl;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Tbasis::write_vectors( std::string fname )
{
    std::fstream fs;
    
    fname += std::string( "_basis.txt");
    fs.open( fname.c_str(), std::fstream::out);
    
	for(int i = 0 ; i < n_vec; i++)
	{
		fs<<i<<":    ";
		for(int j=0; j < n_sites; j++)
		{
			fs<<basis_vectors[i][j]<<" ";
	    }
	    fs<<std::endl;
	}
	
	fs.close();
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
int Tbasis::comp_n_vec( ) 
{
    int n = n_part;
    int m = n_sites;

    long double d=1.,n1,n2;
	int l=m-1.;
	for(int i=n+1.;i<n+m;i++){
	d*=i;
	if(l)d/=i-n;
	}
	
	n_vec = d;	

}
// ========================================================================================





























