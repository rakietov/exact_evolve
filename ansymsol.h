/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ACompSol.cc
   Template functions that exemplify how to print information
   about complex standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ACOMPSOL_H
#define ACOMPSOL_H

#include <cmath>
#include "arcomp.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arlnsmat.h"
#include <algorithm>
#include <iostream>

using namespace std;

template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnz, arcomplex<ARFLOAT> A[], 
              ARINT irow[], ARINT pcol[], arcomplex<ARFLOAT> EigVal[], 
              arcomplex<ARFLOAT>* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric eigen-problems
  on standard "std::cout" stream.
*/

{

  ARINT                                         i;
  arcomplex<ARFLOAT>*                           Ax;
  ARFLOAT*                                      ResNorm;
  ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  std::cout<< std::endl << std::endl << "Testing ARPACK++ function AREig" << std::endl;
  std::cout<< "complex standard eigenvalue problem: A*x - lambda*x \n \n";

  std::cout<< "Dimension of the system            : " << n     << std::endl;
  std::cout<< "Number of 'converged' eigenvalues  : " << nconv << std::endl << std::endl;

  // Printing eigenvalues.

  std::cout<< " with smallest amplitudes:" << std::endl;

  for (i=0; i<nconv; i++) {
    std::cout<< "  lambda[" << (i+1) << "]: " << EigVal[i] << std::endl;
  }
  std::cout<< std::endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/lapy2(real(EigVal[i]),imag(EigVal[i]));
    }

    for (i=0; i<nconv; i++) {
      std::cout<< "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout<< ")*x(" << (i+1) << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout<< std::endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution.


template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnzA, arcomplex<ARFLOAT> A[], 
              ARINT irowA[], ARINT pcolA[], ARINT nnzB, 
              arcomplex<ARFLOAT> B[], ARINT irowB[], ARINT pcolB[],
              arcomplex<ARFLOAT> EigVal[], arcomplex<ARFLOAT>* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problem on standard "std::cout" stream.
*/

{

  ARINT                                        i;
  arcomplex<ARFLOAT>                           *Ax;
  arcomplex<ARFLOAT>                           *Bx;
  ARFLOAT                                      *ResNorm;
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  std::cout<< std::endl << std::endl << "Testing ARPACK++ function AREig" << std::endl;
  std::cout<< "Complex generalized eigenvalue problem: A*x - lambda*B*x \n \n";

  std::cout<< "Dimension of the system            : " << n     << std::endl;
  std::cout<< "Number of 'converged' eigenvalues  : " << nconv << std::endl << std::endl;

  // Printing eigenvalues.

  std::cout<< "Eigenvalues with smallest magnitudes:" << std::endl;

  for (i=0; i<nconv; i++) {
    std::cout<< "  lambda[" << (i+1) << "]: " << EigVal[i] << std::endl;
  }
  std::cout<< std::endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new arcomplex<ARFLOAT>[n];
    Bx      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/lapy2(real(EigVal[i]),imag(EigVal[i]));
    }

    for (i=0; i<nconv; i++) {
      std::cout<< "||A*x(" << i << ") - lambda(" << i;
      std::cout<< ")*B*x(" << i << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout<< std::endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution.


#endif // ACOMPSOL_H

void sum(double nzval1[], double nzval2[], int irow1[], int irow2[],
		int pcol1[],int pcol2[], int nnz1, int nnz2, int d){

	for(int i = d;i>0;i--){
		pcol1[i]-=pcol1[i-1];
		pcol2[i]-=pcol2[i-1];
	}

	int l=0, nnzt=0;

	int *pcolt = new int [d+1];
	double *nzvalt = new double [nnz1+nnz2+2];
	int *irowt = new int [nnz1+nnz2+2];


	for(int i=0;i<d+1;i++)pcolt[i] = 0;
	for(int i=0;i<nnz1+nnz2+2;i++){
		nzvalt[i] = 0;
		irowt[i] = 0;
	}

	int r1 = 0, r2 = 0;
while(l<=d  && r1<nnz1 && r2 <nnz2){
/*
		cout<<"l: "<<l<<endl;
	cout<<endl<<"pcol1"<<endl;
	for(int i = 0;i<=d;i++){
		cout<<pcol1[i]<<" ";

	}
	cout<<endl<<"pcol2"<<endl;
	for(int i = 0;i<=d;i++){
		cout<<pcol2[i]<<" ";
	}
	cout<<endl<<"r1 "<<r1<<endl<<"r2 "<<r2<<endl;

	cout<<endl<<"irow1"<<endl;
	for(int i = 0;i<nnz1;i++){
		cout<<irow1[i]<<" ";

	}
	cout<<endl<<"irow2"<<endl;
	for(int i = 0;i<nnz2;i++){
		cout<<irow2[i]<<" ";
	}
cout<<"tu1"<<endl;
*/
		while((pcol1[l] != 0 || pcol2[l] != 0) && r1<nnz1 && r2 <nnz2){
//		cout<<"wchodze"<<endl;
			bool flag = false;
			if(irow1[r1]<irow2[r2] && pcol1[l]!=0 && !flag){
	//			cout<<1<<endl;
				nzvalt[nnzt] = nzval1[r1];
				irowt[nnzt] = irow1[r1];
				pcolt[l]++;
				nnzt++;
				r1++;
				pcol1[l]--;
				flag = true;
			}
			if(irow1[r1]>irow2[r2]&&pcol2[l]!=0 && !flag){
//				cout<<2<<endl;
				nzvalt[nnzt] = nzval2[r2];
				irowt[nnzt] = irow2[r2];
				pcolt[l]++;
				nnzt++;
				r2++;
				pcol2[l]--;
				flag = true;
			}
			if(irow1[r1]<irow2[r2]&&pcol1[l]==0 && !flag){
//				cout<<3<<endl;
				nzvalt[nnzt] = nzval2[r2];
				irowt[nnzt] = irow2[r2];
				pcolt[l]++;
				nnzt++;
				r2++;
				pcol2[l]--;
				flag = true;
			}
			if(irow1[r1]>irow2[r2] && pcol2[l]==0 && !flag){
//				cout<<4<<endl;
				nzvalt[nnzt] = nzval1[r1];
				irowt[nnzt] = irow1[r1];
				pcolt[l]++;
				nnzt++;
				r1++;
				pcol1[l]--;
				flag = true;
			}
			if(irow1[r1] == irow2[r2] && pcol1[l]!=0&&pcol2[l]!=0 && !flag){
//				cout<<l<<" "<<r1<<" "<<r2<<" "<<nnzt<<" "<<pcol1[l]<<" "<<pcol2[l]<<endl;
//				cout<<5<<endl;cout<<l<<endl;
				nzvalt[nnzt] = nzval1[r1]+nzval2[r2];
				irowt[nnzt] = irow2[r2];
				pcolt[l]++;
				nnzt++;
				r1++;
				r2++;
				pcol2[l]--;
				pcol1[l]--;
				flag = true;
			}
			if(irow1[r1] == irow2[r2] && pcol1[l]==0&&pcol2[l]!=0 && !flag){
//					cout<<6<<endl;
				nzvalt[nnzt] = nzval2[r2];
				irowt[nnzt] = irow2[r2];
				pcolt[l]++;
				nnzt++;
				r2++;
				pcol2[l]--;
				flag = true;
			}
			if(irow1[r1] == irow2[r2] && pcol1[l]!=0&&pcol2[l]==0 && !flag){
//				cout<<7<<endl;
				nzvalt[nnzt] = nzval1[r1];
				irowt[nnzt] = irow1[r1];
				pcolt[l]++;
				nnzt++;
				r1++;
				pcol1[l]--;
				flag = true;
			}
		}
	l++;
	}

l=0;
while(r1<nnz1 && r2 >=nnz2){
	if(l>d)l--;
/*
		cout<<"l1: "<<l<<endl;
	cout<<endl<<"pcol1"<<endl;
	for(int i = 0;i<=d;i++){
		cout<<pcol1[i]<<" ";

	}
	cout<<endl<<"pcol2"<<endl;
	for(int i = 0;i<=d;i++){
		cout<<pcol2[i]<<" ";
	}
	cout<<endl<<"r1 "<<r1<<endl<<"r2 "<<r2<<endl;

	cout<<endl<<"irow1"<<endl;
	for(int i = 0;i<nnz1;i++){
		cout<<irow1[i]<<" ";

	}
	cout<<endl<<"irow2"<<endl;
	for(int i = 0;i<nnz2;i++){
		cout<<irow2[i]<<" ";
	}
cout<<endl;
*/

		while((pcol1[l] != 0 || pcol2[l] != 0) && r1<nnz1){
//				cout<<1<<endl;
				nzvalt[nnzt] = nzval1[r1];
				irowt[nnzt] = irow1[r1];
				pcolt[l]++;
				nnzt++;
				r1++;
				pcol1[l]--;
		}
	l++;
	}

l=0;
while(r1 >= nnz1 && r2 <nnz2){
	if(l>d)l--;
/*
	cout<<endl<<"l2: "<<l;
	cout<<endl<<"pcol1"<<endl;
	for(int i = 0;i<=d;i++){
		cout<<pcol1[i]<<" ";

	}
	cout<<endl<<"pcol2"<<endl;
	for(int i = 0;i<=d;i++){
		cout<<pcol2[i]<<" ";
	}
	cout<<endl<<"r1 "<<r1<<endl<<"r2 "<<r2<<endl;

	cout<<endl<<"irow1"<<endl;
	for(int i = 0;i<nnz1;i++){
		cout<<irow1[i]<<" ";

	}
	cout<<endl<<"irow2"<<endl;
	for(int i = 0;i<nnz2;i++){
		cout<<irow2[i]<<" ";
	}
cout<<endl;
*/

		while((pcol1[l] != 0 || pcol2[l] != 0) && r2<nnz2){
//				cout<<2<<endl;
				nzvalt[nnzt] = nzval2[r2];
				irowt[nnzt] = irow2[r2];
				pcolt[l]++;
				nnzt++;
				r2++;
				pcol2[l]--;
		}
	l++;
	}


//cout<<endl;
	for(int i=0;i<nnzt;i++){
		nzval1[i] = nzvalt[i];
		irow1[i] = irowt[i];
//		cout<<"nzval1: "<<nzval1[i]<<" irow1: "<<irow1[i]<<endl;
	}
	for(int i=0;i<=d;i++){
		pcol1[i] = pcolt[i];

	}
	for(int i=0;i<=d;i++){
		pcol1[i+1] += pcol1[i];
	}
	for(int i=0;i<=d;i++){
//		cout<<"pcol1: "<<pcol1[i]<<endl;
	}



	delete[] pcolt;
	delete[] nzvalt;
	delete[] irowt;

}
//Testy na dodawanie macierzy:
/*

	int irow1[] = {0, 3, 1, 0, 2, 0, 3};
	int pcol1[] = {0, 2, 3, 5, 7};
	double nzval1[] = {-1,1,3,2,-4,5,2};

	int irow2[] = {1, 3, 1, 0, 2, 0, 3};
	int pcol2[] = {0, 2, 3, 5, 7};
	double nzval2[] = {-2,9,3,2,-4,5,2};

*/

/*
	int irow1[] = {3,3};
	int pcol1[] = {0,1,1,1,2};
	double nzval1[] = {1,4};

	int irow2[] = {0,0};
	int pcol2[] = {0, 1,1,1,2};
	double nzval2[] = {1,2};
	sum(nzval1, nzval2,irow1,  irow2,
		 pcol1, pcol2, 2, 2, 4);
*/




struct ind {
	double V;
	int pos;
};

bool comp(ind a, ind b){
	return (a.V<b.V);
}

double tag (int* A,int m ){
	double tmp = 0;
	for(int j=0;j<m;j++)tmp+=sqrt(100*j+3)*A[j];
	return tmp;
}

double anih( int i, int *A){
	double tmp1 = A[i];	

	if(A[i] > 0)A[i]--;
	else return 0;

	return sqrt(tmp1);
}





double ada( int i, int j, int *A){
	double tmp1 = A[i];	
	double tmp2 = A[j];


	A[i]++;

	if(A[j] > 0)A[j]--;
	else return 0;

	if( i == j) return A[i];
	else return sqrt((tmp1+1)*(tmp2));
}

double adadaa(int i,int j,int k, int l, int *A){

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


int comp_dim(int n,int m){
	double d=1,n1,n2;
	int l=m-1;
	for(int i=n+1;i<n+m;i++){
	d*=i;
	if(l)d/=i-n;
	}
	return d;	
}


struct vir {
	double V;
	int row;
};

bool comp_vir (vir a,vir b){
	return (a.row<b.row);
}



void przesortuj(double nzval[], int irow[],int pcol[],int d, int nnz){

	vir * A = new vir[10*nnz];
	
	for(int i=0;i<nnz;i++){
		A[i].V = nzval[i];
		A[i].row = irow[i];//cout<<irow[i]<<endl;
	}
	
	int l = 0;
	while(l<d){
//		cout<<"Sortuje od "<<pcol[l]<<" do "<<pcol[l+1]<<endl;
//		cout<<"czyli liczby o rowach: ";
//		for(int i=pcol[l];i<pcol[l+1];i++)cout<<irow[i]<<" ";cout<<endl;
		sort(A+pcol[l],A+pcol[l+1],comp_vir);
//		cout<<"otrzymując: ";
//		for(int i=pcol[l];i<pcol[l+1];i++)cout<<A[i].row<<" ";cout<<endl;
		l++;
	}

	for(int i=0;i<nnz;i++){
		nzval[i] =A[i].V ;
		irow[i]=A[i].row  ;
	}

	delete[] A;

}


int reduce(double nzval[], int irow[],int pcol[],int d, int nnz){

	int * irowO = new int[2*nnz];
	int * pcolO = new int[2*d+2];
	double * nzvalO = new double[2*nnz];
	int nnzO = 0;

	for(int i=0;i<nnz+2;i++){
		irowO[i] = 0;
		nzvalO[i] = 0.;
	}
	for(int i=0;i<d+2;i++)pcolO[i] = 0;


	for(int i = d;i>0;i--){
		pcol[i]-=pcol[i-1];
	}
	

	
	int l = 0,r=0,rO=0;
	while(l<=d && r<nnz){
		while(pcol[l+1]!=0){
			int n=0;
			while(pcol[l+1]!=0 && irow[r] == irow[r+n+1]){
				nzval[r] +=nzval[r+n+1];
				n++;
				pcol[l+1]--;					
				}
			nzvalO[rO] = nzval[r]; 
			pcol[l+1]--;
			irowO[rO] = irow[r];
//			cout<<"dodaje na rO: "<<rO<<" nzvalO: "<<nzvalO[rO]<<" irowO: "<<irowO[rO]<<endl;
			pcolO[l+1]++;
	//		cout<<"pcolO: "<<pcolO[l+1]<<"dla l+1 = "<<l+1<<endl;
			rO++;
			r+= 1+n;
			nnzO++;
		}
		l++;
	}
	
	for(int i=0;i<nnzO;i++){
		nzval[i] = nzvalO[i];
		irow[i] = irowO[i];
	}
	for(int i=0;i<d+2;i++){
		pcolO[i+1]+=pcolO[i];

	}

	for(int i=0;i<d+2;i++){
		pcol[i]+=pcolO[i];
	}

	delete[] pcolO;
	delete[] nzvalO;
	delete[] irowO;


	return nnzO;
	

}

//void wypisz(ofstream &out, double nzval[], int irow[], int pcolO[], int nnz, int d){
//  out << 5.6;
  
//}

/*

wypisz(cout,....)

ofstream file;
file.open("matrix.dat)"
wypisz(file,...)
file.close();


*/

void wypisz(double nzval[], int irow[], int pcolO[], int nnz, int d){

int l=0,kol=0;
float zero = 0.;
	
	int * pcol = new int[d+5];
	
	for(int i=0;i<d+3;i++)pcol[i] = pcolO[i];

	for(int i = d;i>0;i--){
		pcol[i]-=pcol[i-1];
	}

for(int j=0;j<d;j++){
	bool flag = true;
	for (int i=0;i<d;i++){
		if(irow[l] == i && flag){
			printf("  %.10f   ",nzval[l]);
			l++;
			pcol[j+1]--;
		}
		else printf("   %.10f   ",zero);
		if(pcol[j+1] == 0)flag = false;
		}
	cout<<endl;
	}


	delete[] pcol;

}


void wypisz2(double nzval[], int irow[], int pcolO[], int nnz, int d){

int l=0,r=0;
float zero = 0.;
	
	int * pcol = new int[d+5];
	
	for(int i=0;i<d+3;i++)pcol[i] = pcolO[i];

	for(int i = d;i>0;i--){
		pcol[i]-=pcol[i-1];
	}

while( r<=d ){
	while(pcol[r] != 0){
		printf("%d %d %.6f\n",irow[l]+1,r,nzval[l]);
		pcol[r]--;
		l++;
	}
	r++;
}
	delete[] pcol;
}









