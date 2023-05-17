#include <iostream>
#include <valarray>
#include <cmath>

#include "grid_gen.h"
#include "phys_par.h"
#include "grid_par.h"

using std::cout;
using std::cin;
using std::valarray;


extern void print_valarray_length(const valarray<double>&, int l =10);
extern void CreateWeightMatrix(const char*);

extern double weight_matrix[4][3];
/*extern*/ struct scheme_type_error;

double f_alpha(double r)
{
	return grid_kappa(r)*pow(r,NN - 1)*p1(r);
};

double f_beta(double r)
{
	double h_x = use_h_x();			//without this, we'd get h_x = 0 
									//instead of the correct value of h_x 
//	cout << "\n\tcontrol: in f_beta h_x = " << h_x <<'\n'; 
	return p1(r)/(pow(r, NN -1))/h_x/h_x;
};

double f_pow(double r){return pow(r,NN - 1);};

valarray<double>LambdaX(const valarray<double> U,
	const valarray<double> L0,
	const valarray<double> Lpls,
	const valarray<double> Lmns
)
{
	return L0 * (U)+
		Lpls * (U.shift(1)) +
		Lmns * (U.shift(-1));
};


/********************************************************************************/
/*				Create grid operators for the interval [rb0, r0]				*/

class GridOperators
{
public:
	valarray<double>Epsilon;

	valarray<double>LF00;
	valarray<double>LF10;
	valarray<double>LF20;

	valarray<double>LF0pls;
	valarray<double>LF1pls;
	valarray<double>LF2pls;

	valarray<double>LF0mns;
	valarray<double>LF1mns;
	valarray<double>LF2mns;
/****************************/
/*		For tests only		*/
	valarray<double>TestL0;
	valarray<double>TestLpls;
	valarray<double>TestLmns;
/*							*/
/****************************/
	GridOperators();
	~GridOperators();
	
};

GridOperators::GridOperators():
	Epsilon(K+2),
	LF00(K+2),   LF10(K+2),   LF20(K+2),
	LF0pls(K+2), LF1pls(K+2), LF2pls(K+2),
	LF0mns(K+2), LF1mns(K+2), LF2mns(K+2)
	,TestLpls(K+2),TestLmns(K+2),TestL0(K+2)
	{
	cout << "constructor for GridOperators" <<'\n';
	double h_x = use_h_x();
	double h_t = use_h_t();
	cout << "\tcontrol: in GridOperators h_x = " << h_x << "\n"; 
	CreateWeightMatrix(type_of_scheme);
	const valarray<double> Alpha = Tab_r_semiint.apply(f_alpha);
//	cout << "Alpha: "; print_valarray_length(Alpha);
	const valarray<double> Beta = Tab_r_int.apply(f_beta);
	const valarray<double> Delta = 
		Tab_r_int.apply(grid_sigma)/
		h_t/
		(weight_matrix[1][0] - weight_matrix[1][2]);
	const valarray<double> Epsilon = 
		h_t*h_t/Tab_r_int.apply(grid_rho);

	const valarray<double> Lpls = Beta*Alpha;
	const valarray<double> Lmns = Beta*(Alpha.shift(-1));
	const valarray<double> L0 = - Lpls - Lmns;
/****************************/
/*		For tests only		*/
	TestL0 = L0;
	TestLpls = Lpls;
	TestLmns = Lmns;

	for (int ii = 0; ii<=3; ii++)
	{
		for (int jj = 0; jj<=2; jj++)
		{
			cout << weight_matrix[ii][jj] <<'\t';
		};
		cout<< '\n';
	};
/*							*/
/****************************/


//	print_valarray_length(Lpls);
//	print_valarray_length(Lmns);
//	print_valarray_length(L0);

	LF00 =  Epsilon*(
		weight_matrix[0][0] * L0 -
		weight_matrix[1][0] * Delta -
		weight_matrix[2][0] * gamma/(Tab_r_int.apply(f_pow)) -			// !!!! See the original
		weight_matrix[3][0] * (Tab_r_int).apply(grid_q)
		);
	LF10 = (  valarray<double>(2, K+2) +
			Epsilon*(
		weight_matrix[0][1] * L0 -
		weight_matrix[1][1] * Delta -
		weight_matrix[2][1] * gamma/(Tab_r_int.apply(f_pow)) -			// !!!! See the original
		weight_matrix[3][1] * (Tab_r_int).apply(grid_q)
		)
		);
	LF20 =(  valarray<double>(-1, K+2) +
			Epsilon*(
		weight_matrix[0][2] * L0 -
		weight_matrix[1][2] * Delta -
		weight_matrix[2][2] * gamma/(Tab_r_int.apply(f_pow)) -			// !!!! See the original
		weight_matrix[3][2] * (Tab_r_int).apply(grid_q)
		)
		);

	LF0pls = Epsilon*
		weight_matrix[0][0] * Lpls;
	LF1pls = Epsilon*
		weight_matrix[0][1] * Lpls;
	LF2pls = Epsilon*
		weight_matrix[0][2] * Lpls;

	LF0mns = Epsilon*
		weight_matrix[0][0] * Lmns;
	LF1mns = Epsilon*
		weight_matrix[0][1] * Lmns;
	LF2mns = Epsilon*
		weight_matrix[0][2] * Lmns;
	

//	print_valarray_length(Lmns);
//	print_valarray_length(L0);

	

};
GridOperators::~GridOperators(){
	cout << "destructor for GridOperators" <<'\n';
};


