#include <iostream>
#include <valarray>
#include <math.h>


using std::cout;
using std::cin;
using std::valarray;


extern int K;
//extern double CRight1;
//extern double CLeft1;

void F_IterStep_Pointer(const valarray<double>* P_U0Prev,
							  valarray<double>* P_U0New,
							   const valarray<double>& U12F,
							   const double& UIFRight,
							   const double& UIFLeft,
							   const valarray<double>& LF00,
							   const valarray<double>& LF0pls,
							   const valarray<double>& LF0mns,
							   const double& CRight1,
							   const double& CLeft1
							  )
{
	(*P_U0New) = 
		LF00 * (*P_U0Prev) +
		LF0pls * ( P_U0Prev ->shift( 1) ) +
		LF0mns * ( P_U0Prev ->shift(-1) ) +
		U12F;

	(*P_U0New)[K+1] =
		CRight1 * ( (*P_U0New)[K] ) +
		UIFRight;
	(*P_U0New)[0] =
		CLeft1 * ( (*P_U0New)[1] ) +
		UIFLeft;
};
/****************************************************************************/
/*	The iteration convergence criterion 
	is based on the norm of the grid space C.				 				*/

extern bool IterCrit(const valarray<double>& V1,
					 const valarray<double>& V2, 
					 const double iteration_precision
					 )
{

	return ( ( (V1 - V2).apply(fabs) ).max() < iteration_precision);
};
