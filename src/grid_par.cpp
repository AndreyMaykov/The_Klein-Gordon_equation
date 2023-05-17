#include <iostream>
#include <valarray>

#include "grid_gen.h"
#include "phys_par.h"

using std::cout;
using std::cin;
using std::valarray;

int K = 7;
double fc = 0.95;
int nn_0 = 2;
int nn_max = 18; 
extern const double iteration_precision = 0.000000001;

extern const double rb0;
extern const double r0;
const char* type_of_scheme = "conservative";		// See the original 
	
extern double fx(double r){return r;};				
extern double p1(double r){return 1;};				
extern const double h_x = (fx(r0) - fx(rb0))/double(K);

double& use_h_x()									// See the original
{												
	static double h_x = (fx(r0) - fx(rb0))/K;	
	return h_x;									
};

extern const double h_t = fc * h_x;

double& use_h_t()								
{												
	static double h_t = fc * (fx(r0) - fx(rb0))/K;	
	return h_t;									
};



/*
phys_par ph_p = phys_par();
phys_par ph_p = use_ph_p();
const double a_rb0 = ph_p.A_rb0;
const double b_rb0 = ph_p.B_rb0;
const double c_rb0 = ph_p.C_rb0;
const double a_r0 = ph_p.A_r0;
//const double b_r0 = ph_p.B_r0;
const double c_r0 = ph_p.C_r0;
*/


const double a_rb0 = use_a_rb0();
const double b_rb0 = use_b_rb0();
const double c_rb0 = use_c_rb0();
const double a_r0 = use_a_r0();
const double b_r0 = use_b_r0();
const double c_r0 = use_c_r0();


const double CL0 =   a_rb0/2.0/h_t - p1(rb0) * b_rb0/2.0/h_x + c_rb0/4.0;
const double CL1 = - a_rb0/2.0/h_t - p1(rb0) * b_rb0/2.0/h_x - c_rb0/4.0;
const double CL2 =   a_rb0/2.0/h_t + p1(rb0) * b_rb0/2.0/h_x - c_rb0/4.0;
const double CL3 =   a_rb0/2.0/h_t - p1(rb0) * b_rb0/2.0/h_x - c_rb0/4.0;

const double CR0 =   a_r0/2.0/h_t + p1(r0) * b_r0/2.0/h_x + c_r0/4.0;
const double CR1 = - a_r0/2.0/h_t + p1(r0) * b_r0/2.0/h_x - c_r0/4.0;
const double CR2 =   a_r0/2.0/h_t - p1(r0) * b_r0/2.0/h_x - c_r0/4.0;
const double CR3 =   a_r0/2.0/h_t + p1(r0) * b_r0/2.0/h_x - c_r0/4.0;

extern const double CRight0 = 1.0/CR0;
extern const double CRight1 = CR1/CR0;
extern const double CRight2 = CR2/CR0;
extern const double CRight3 = CR3/CR0;

extern const double CLeft0 = 1.0/CL0;
extern const double CLeft1 = CL1/CL0;
extern const double CLeft2 = CL2/CL0;
extern const double CLeft3 = CL3/CL0;

struct GridData
{
	valarray<double> U2init;
	valarray<double> U1init;
	GridData();
	~GridData();
};
GridData::GridData():
U1init(0.0,K+2),
U2init(0.0,K+2)
{
	cout << "constructor for GridData" << '\n';
	for(int k = 1; k <= K; k++)
	{
		U1init[k] = u_0_init(Tab_r_int[k]) + (h_t/2)*u_1_init(Tab_r_int[k]);
		U2init[k] = u_0_init(Tab_r_int[k]) - (h_t/2)*u_1_init(Tab_r_int[k]);
	};
	U1init[0] = 2*u_0_init(Tab_r_semiint[0]) - u_0_init(Tab_r_int[1]) +
		(h_t/2)*( 2*u_1_init(Tab_r_semiint[0]) - u_1_init(Tab_r_int[1]) );
	U2init[0] = 2*u_0_init(Tab_r_semiint[0]) - u_0_init(Tab_r_int[1]) -
		(h_t/2)*( 2*u_1_init(Tab_r_semiint[0]) - u_1_init(Tab_r_int[1]) );
	U1init[K+1] = 2*u_0_init(Tab_r_semiint[K]) - u_0_init(Tab_r_int[K]) +
		(h_t/2)*( 2*u_1_init(Tab_r_semiint[K]) - u_1_init(Tab_r_int[K]) );
	U2init[K+1] = 2*u_0_init(Tab_r_semiint[K]) - u_0_init(Tab_r_int[K]) -
		(h_t/2)*( 2*u_1_init(Tab_r_semiint[K]) - u_1_init(Tab_r_int[K]) );
};
GridData::~GridData()
{
	cout << "destructor for GridData" << '\n';
};

