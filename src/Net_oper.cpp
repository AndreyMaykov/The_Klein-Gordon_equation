#include <iostream.h>
#include <valarray>
#include "Net_Gen.h"
#include "phys_par.h"
#include "Net_par.h"
#include <cmath>


using std::valarray;


extern void print_valarray_length(const valarray<double>&, int l =10);
extern void CreateWeightMatrix(char*);

extern double weight_matrix[4][3];
/*extern*/ struct scheme_type_error;

double f_alpha(double r)
{
	return net_kappa(r)*pow(r,NN - 1)*p1(r);
};

double f_beta(double r)
{
	double h_x = use_h_x();			//без этого приема получим h_x = 0 
									//вместо "правильного" значения шага 
//	cout << "\n\tcontrol: in f_beta h_x = " << h_x <<'\n'; 
	return p1(r)/(pow(r, NN -1))/h_x/h_x;
};

double f_pow(double r){return pow(r,NN - 1);};


/********************************************************************************/
/*				Генерация сеточных операторов на отрезке [rb0, r0]				*/

class NetOperators
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
/*	Только для теста		*/
	valarray<double>TestL0;
	valarray<double>TestLpls;
	valarray<double>TestLmns;
/*							*/
/****************************/
	NetOperators();
	~NetOperators();
	
};

NetOperators::NetOperators():
	Epsilon(K+2),
	LF00(K+2),   LF10(K+2),   LF20(K+2),
	LF0pls(K+2), LF1pls(K+2), LF2pls(K+2),
	LF0mns(K+2), LF1mns(K+2), LF2mns(K+2)
	,TestLpls(K+2),TestLmns(K+2),TestL0(K+2)
	{
	cout << "constructor for NetOperators" <<'\n';
	double h_x = use_h_x();
	double h_t = use_h_t();
	cout << "\tcontrol: in NetOperators h_x = " << h_x << "\n"; 
	CreateWeightMatrix(type_of_scheme);
	const valarray<double> Alpha = Tab_r_semiint.apply(f_alpha);
//	cout << "Alpha: "; print_valarray_length(Alpha);
	const valarray<double> Beta = Tab_r_int.apply(f_beta);
	const valarray<double> Delta = 
		Tab_r_int.apply(net_sigma)/
		h_t/
		(weight_matrix[1][0] - weight_matrix[1][2]);
	const valarray<double> Epsilon = 
		h_t*h_t/Tab_r_int.apply(net_rho);

	const valarray<double> Lpls = Beta*Alpha;
	const valarray<double> Lmns = Beta*(Alpha.shift(-1));
	const valarray<double> L0 = - Lpls - Lmns;
/****************************/
/*	Только для теста		*/
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
		weight_matrix[2][0] * gamma/(Tab_r_int.apply(f_pow)) -			///!!!! ошибка !!!! потеряна степень r
		weight_matrix[3][0] * (Tab_r_int).apply(net_q)
		);
	LF10 = (  valarray<double>(2, K+2) +
			Epsilon*(
		weight_matrix[0][1] * L0 -
		weight_matrix[1][1] * Delta -
		weight_matrix[2][1] * gamma/(Tab_r_int.apply(f_pow)) -			///!!!! ошибка !!!! потеряна степень r
		weight_matrix[3][1] * (Tab_r_int).apply(net_q)
		)
		);
	LF20 =(  valarray<double>(-1, K+2) +
			Epsilon*(
		weight_matrix[0][2] * L0 -
		weight_matrix[1][2] * Delta -
		weight_matrix[2][2] * gamma/(Tab_r_int.apply(f_pow)) -			///!!!! ошибка !!!! потеряна степень r
		weight_matrix[3][2] * (Tab_r_int).apply(net_q)
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
NetOperators::~NetOperators(){
	cout << "destructor for NetOperators" <<'\n';
};


