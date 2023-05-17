#include <iostream>
#include <valarray>

#include "phys_par.h"	
#include "grid_par.h"	

using std::cout;
using std::cin;
using std::valarray;

//extern int K;
//extern const double h_r = (r0 - rb0)/K;
//extern valarray<double> Pro_r_int;			//keep in case Pro_r_int and 
//extern valarray<double> Pro_r_semiint;		//Pro_r_semiint are read from a file
//extern void CreateWeightMatrix(char*);
//extern double weight_matrix[4][3];

//extern void CreateWeightMatrix(char* type_of_scheme);
// /*extern*/ struct scheme_type_error;

/********************************************************************************/
/*				Generate uniform grids on the interval [rb0, r0]				*/

struct GridsUni
{
	std::valarray<double> Pro_r_int;
	std::valarray<double> Pro_r_semiint;
	int P;
	GridsUni();
	~GridsUni();
};


GridsUni::GridsUni():Pro_r_int(K + 2),Pro_r_semiint(K + 2)
{
	cout << "constructor for GridsUni"<<"\n";

/************************************************************************************/

	cout << "K = "<<K <<"\n";
	double h_r = (r0 - rb0)/K;
	cout << "h_r = "<< h_r <<"\n\n";
	Pro_r_int = 0;
	Pro_r_semiint = 0;

//	cout << "Pro_r_int[1] = " << Pro_r_int[1] <<"\n\n";
	for(int k = 0; k < K + 2; k++)
	{
		Pro_r_int[k] = rb0 - (h_r)/2 + (h_r)*k;
	};

	for(int k = 0; k < K + 2; k++)
	{
		Pro_r_semiint[k] = rb0 + (h_r)*k;
	};
//	cout << "Pro_r_int[1] = " << Pro_r_int[1] <<"\n\n";
//	cout << "grid_rho(rb0) = " << grid_rho(rb0) <<'\n';

/************************************************************************************/
};
GridsUni::~GridsUni(){cout<<"destructor for GridsUni"<<"\n\n";}


GridsUni Grids;

extern const std::valarray<double> Tab_r_int = Grids.Pro_r_int;
extern const valarray<double>Tab_r_semiint = Grids.Pro_r_semiint;

/*																				*/
/********************************************************************************/