#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>



#ifndef IND_H_
#define IND_H_

class Tind 
{
public:
    double V;
    int pos;
    void tag ( std::vector< int > vec, int i );
};
/*
// ----------------------------------------------------------------------------------------
bool comp_ind(Tind a, Tind b){
	return (a.V<b.V);
}
// ========================================================================================
*/

#endif 
