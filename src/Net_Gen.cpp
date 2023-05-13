#include <iostream.h>
#include <valarray>
#include "phys_par.h"	
#include "Net_par.h"	

using std::valarray;

//extern int K;
//extern const double h_r = (r0 - rb0)/K;
//extern valarray<double> Pro_r_int;			//оставить на тот случай, если Pro_r_int и 
//extern valarray<double> Pro_r_semiint;		//Pro_r_semiint будут считыватьс€ из файла
//extern void CreateWeightMatrix(char*);
//extern double weight_matrix[4][3];

//extern void CreateWeightMatrix(char* type_of_scheme);
// /*extern*/ struct scheme_type_error;

/********************************************************************************/
/*				√енераци€ равномерных сеток на отрезке [rb0, r0]				*/

struct NetsUni
{
	std::valarray<double> Pro_r_int;
	std::valarray<double> Pro_r_semiint;
	int P;
	NetsUni();
	~NetsUni();
};


NetsUni::NetsUni():Pro_r_int(K + 2),Pro_r_semiint(K + 2)
{
	cout << "constructor for NetsUni"<<"\n";

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

	for(k = 0; k < K + 2; k++)
	{
		Pro_r_semiint[k] = rb0 + (h_r)*k;
	};
//	cout << "Pro_r_int[1] = " << Pro_r_int[1] <<"\n\n";
//	cout << "net_rho(rb0) = " << net_rho(rb0) <<'\n';

/************************************************************************************/
};
NetsUni::~NetsUni(){cout<<"destructor for NetsUni"<<"\n\n";}


NetsUni Nets;

extern const std::valarray<double> Tab_r_int = Nets.Pro_r_int;
extern const valarray<double>Tab_r_semiint = Nets.Pro_r_semiint;

/*																				*/
/********************************************************************************/