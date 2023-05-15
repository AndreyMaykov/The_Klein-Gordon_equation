#include <iostream>
#include <valarray>
#include <cmath>

#include "GridGen.h"
#include "FDiffOper.h"
#include "IterSt.h"
#include "GridPar.h"
#include "PhysPar.h"
#include "PrintVal.h"

using std::valarray;
using std::cout;
using std::cin;

extern void print_valarray_length(const valarray<double>&,int = 10);
extern double weight_matrix[4][3];
char character_for_pause;


void main()
{
	cout<<"\t*** start main ***"<<'\n';

	GridOperators* P_GridOp = new GridOperators;

	const valarray<double> Epsilon = (P_GridOp->Epsilon);

	const valarray<double> LF00 = (P_GridOp->LF00);
	const valarray<double> LF10 = (P_GridOp->LF10);
	const valarray<double> LF20 = (P_GridOp->LF20);

	const valarray<double> LF0pls = (P_GridOp->LF0pls);
	const valarray<double> LF1pls = (P_GridOp->LF1pls);
	const valarray<double> LF2pls = (P_GridOp->LF2pls);

	const valarray<double> LF0mns = (P_GridOp->LF0mns);
	const valarray<double> LF1mns = (P_GridOp->LF1mns);
	const valarray<double> LF2mns = (P_GridOp->LF2mns);


	delete P_GridOp;

/********************************************************************************/
/*					Set initial conditions U2init и U1init					*/   

	valarray<double>U2init(0.0, K + 2);
	valarray<double>U1init(0.0, K + 2);
	GridData* P_GridDat = new GridData;
	U1init = P_GridDat->U1init;	//see original
	U2init = P_GridDat->U2init;	//see original
	delete P_GridDat;
/********************************************************************************/
/*							Only for testing with NN != 1						*/

//						For NN = 2;	gammma =0; K = 5; fc = 2;
//	{0.2367407537764706,-0.2367407537764706,-0.5496392412230395,
//	-0.6257641411610598,-0.4715947967519451,-0.1669374353539749,0.1669374353539749};
//					      omega_tild_Vec =  0.3043594240208501
	const double Vec[] = {0.6464850627430188,0.6464850627430188,0.5370050592066881,0.3519290866767011,
  0.1342298627898095,-0.07144432577829789,-0.2271204346311821,-0.\
3084200194700978,-0.3084200194700978};
	double omega_tild_Vec =  0.3200243368887432;
	std::valarray<double>ValVec(Vec, K+2);
//						Analytical solution: sin(omega_tild*t)*ValVec				
	
/*																				*/
/********************************************************************************/

	double SIN = cos(-omega_tild_Vec*h_t/2);
	U2init =   SIN*ValVec;
	U1init =   SIN*(ValVec);
	double SIN_nn;
/*																				*/
/********************************************************************************/
	cout << "U2init: ";
	cout << "omega_tild_Vec = "<< omega_tild_Vec <<'\n';
	cin >>character_for_pause;
	print_valarray_length(U2init,10);
	cout << "U1init: ";
	print_valarray_length(U1init,10);
	cin >> character_for_pause;

	valarray<double>U0Prev(0.0, K + 2);
	valarray<double>U0New(0.0, K + 2);
	valarray<double>VF_internal(0.0, K+2);

	valarray<double>U12F(0.0, K + 2);
	double UIFRight;
	double UIFLeft;

	valarray<double>* P_U2 = &U2init;
	valarray<double>* P_U1 = &U1init;
	valarray<double>* P_U0Prev = &U0Prev;
	valarray<double>* P_U0New = &U0New;
	valarray<double>* temp_adress = P_U0Prev;
	
/********************************************************************************************/
/*							Loop variable n: beginning of the loop							*/

	for (int nn = nn_0; nn <= nn_max ; nn++)
	{
		cout<< "nn = " << nn <<'\n';
		
		for (int k = 1; k<=K; k++)
		{
			VF_internal[k] = f_internal(Tab_r_int[k], (nn-1.5)*h_t);
		};
		U12F = LF10 * (*P_U1) + LF1pls * (P_U1->shift(1)) + LF1mns * (P_U1->shift(-1))+
			   LF20 * (*P_U2) + LF2pls * (P_U2->shift(1)) + LF2mns * (P_U2->shift(-1))+	
			   Epsilon * VF_internal;
		UIFRight = CRight2 * (*P_U1)[K+1] + CRight3 *(*P_U1)[K] + 
			CRight0 * (
			/* 		terms from the open boundary conditions (см. с. -07_6- ) 			*/ 
			f_r0( (nn - 1)*h_t )
			);
		UIFLeft = CLeft2 * (*P_U1)[0] + CLeft3 * (*P_U1)[1] + 
			CLeft0 * (
			/* 		terms from the open boundary conditions (см. с. -07_6- ) 			*/ 
			f_rb0( (nn - 1)*h_t )
			);

		F_IterStep_Pointer(P_U1, P_U0New, 
						   U12F, UIFRight, UIFLeft,
					       LF00, LF0pls, LF0mns,
					       CRight1, CLeft1
					  );

		cout << "!IterCrit(*P_U0Prev, *P_U0New, iteration_precision) =" <<
       		   ( !IterCrit(*P_U0Prev, *P_U0New, iteration_precision) ) << '\n';
			   
/****************************************************************************************/
/*						Fixed-point iteration loop -- beginning							*/

		if ( ((weight_matrix[0][0] != 0) || (weight_matrix[0][2] != 0 ))
			 &
			 ( !IterCrit(*P_U1, *P_U0New, iteration_precision))
		   )
		{
			int jj = 0;
//			cout << "!IterCrit(*P_U0Prev, *P_U0New, iteration_precision) =" <<
//				( !IterCrit(*P_U0Prev, *P_U0New, iteration_precision) ) << '\n';
			do 
			{
				temp_adress = P_U0New;
				P_U0New = P_U0Prev;
				P_U0Prev = temp_adress;

				F_IterStep_Pointer(P_U0Prev, P_U0New, 
								   U12F, UIFRight, UIFLeft,
								   LF00, LF0pls, LF0mns,
								   CRight1, CLeft1
								   );
				jj++;
//				cout << "jj = " << jj << '\n';
/*				Control = ((*P_U0New )+ (*P_U2) - 2.0 * (*P_U1))/h_t/h_t - 
//					LambdaX(U1, L0, Lpls, Lmns);
					(
					((*P_U0Prev).shift(-1) + (*P_U0Prev).shift(1) - 2.0 * (*P_U0Prev))/h_x/h_x/4.0 +
					((*P_U1).shift(-1) + (*P_U1).shift(1) - 2.0 * (*P_U1))/h_x/h_x/2.0 +
					((*P_U2).shift(-1) + (*P_U2).shift(1) - 2.0 * (*P_U2))/h_x/h_x/4.0
					); 
				cout << "     Control for iterations";
				print_valarray_length(Control,10);
				cin >> character_for_pause;
*/
			}
			while( !IterCrit(*P_U0Prev, *P_U0New, iteration_precision) );
		};

/*							Fixed-point iteration loop -- end							*/
/****************************************************************************************/

// The Control valarray was removed -- see the original

		cout << "U0 :";
		print_valarray_length(*P_U0New);
		 
		valarray<double> u_prec0(0.0, K+2);
		for (int k = 0; 
					k<=K+1; 
					k++){
			u_prec0[k] = u_static3(Tab_r_int[k],(nn-1.0/2.0)*h_t);
		};
		SIN_nn = cos(omega_tild_Vec*(nn - 0.5)*h_t);
		u_prec0 = SIN_nn*ValVec;
		valarray<double> u_prec1(0.0, K+2);
		for (int k = 0; 
					k<=K+1; 
					k++){
			u_prec1[k] = u_static3(Tab_r_int[k],(nn-3.0/2.0)*h_t);
		};
		SIN_nn = cos(omega_tild_Vec*(nn - 1.5)*h_t);
		u_prec1 = SIN_nn*ValVec;
		valarray<double> u_prec2(0.0, K+2);
		for (int k = 0; 
					k<=K+1; 
					k++){
			u_prec2[k] = u_static3(Tab_r_int[k],(nn-5.0/2.0)*h_t);
		};
		SIN_nn = cos(omega_tild_Vec*(nn - 2.5)*h_t);
		u_prec2 = SIN_nn*ValVec;
		
//		double omega;
		omega = use_omega();
		cout << "    du/dt + cos(omega*h_x/2)*du/dx|_(r=rb0):"<<'\n';
		cout << (u_prec0[1]+u_prec0[0]-u_prec1[1]-u_prec1[0])/2/h_t +
			cos(omega*h_x/2)*
			(u_prec0[1]-u_prec0[0]+u_prec1[1]-u_prec1[0])/2/h_x
			<< '\n';
		cout << "    du/dt + cos(omega*h_x/2)*du/dx|_(r=r0):"<<'\n';
		cout << (u_prec0[K+1]+u_prec0[K]-u_prec1[K+1]-u_prec1[K])/2/h_t +
			cos(omega*h_x/2)*
			(u_prec0[K+1]-u_prec0[K]+u_prec1[K+1]-u_prec1[K])/2/h_x
			<< '\n';
		cout << "    d(U0New)/dt + cos(omega*h_x/2)*d(U0New)/dx|_(r=r0):"<<'\n';
		cout << ((*P_U0New)[K+1]+(*P_U0New)[K]-(*P_U1)[K+1]-(*P_U1)[K])/2/h_t +
			cos(omega*h_x/2)*
			((*P_U0New)[K+1]-(*P_U0New)[K]+(*P_U1)[K+1]-(*P_U1)[K])/2/h_x
			<< '\n';
		cin >> character_for_pause;
//		cout << "U12F :";
//		print_valarray_length(U12F);
//		print_valarray_length(u_prec0);
		cout << "coeff" ;
		print_valarray_length((*P_U0New)/u_prec0);
		print_valarray_length((*P_U1)/u_prec1);
		print_valarray_length((*P_U2)/u_prec2);
//		print_valarray_length(U0New/u_prec0);
		cout << "    U0New : ";
		print_valarray_length(U0New);
		cin >> character_for_pause;
		temp_adress = P_U2;
		P_U2 = P_U1;
		P_U1 = P_U0New;
		P_U0New = temp_adress;
//		cout << temp_address <<'\t'<<P_U0New<<'\t'<<P_U0Prev<<'\t'<<P_U1<<'\t'<<P_U2<<'\n';
	};
	
/*								Loop variable n: end of the loop							*/
/********************************************************************************************/
	
	

	cout<<"\t*** complete main ***"<<'\n';
}	