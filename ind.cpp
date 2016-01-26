#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>

#include "ind.h"

// ----------------------------------------------------------------------------------------
void Tind::tag ( std::vector< int > vec, int i )
{
	double tmp = 0.;
	for(int j = 0;j < vec.size(); ++j)tmp+=sqrt(100.*j+3.)*vec[j];
    Tind ind1;
    ind1.pos = i;
    ind1.V = tmp;
    
	*this = ind1;
}
// ========================================================================================












