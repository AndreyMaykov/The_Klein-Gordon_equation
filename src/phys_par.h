#pragma once

extern const int NN;
extern const double gamma;		//  !!!!
extern const double rb0;
extern const double r0;

extern double rho(double r);	
extern double kappa(double r);	
extern double sigma(double r);	
extern double q(double r);		

extern double grid_rho(double r);	// = (r<=r0)?rho(r):1;
extern double grid_kappa(double r);	// = (r<=r0)?kappa(r):1;
extern double grid_sigma(double r);	// = (r<=r0)?sigma(r):0;
extern double grid_q(double r);		// = (r<=r0)?q(r):0;

extern double u_0_init(double r);
extern double u_1_init(double r);


//extern const double a_rb0, b_rb0, c_rb0, a_r0, b_r0, c_r0;

extern double f_internal(double r, double t);
extern double f_rb0(double r);
extern double f_r0(double r);
/*
struct phys_par
{
	double A_rb0, B_rb0;
	double C_rb0, A_r0, B_r0, C_r0 ;
	phys_par();
	~phys_par();
};
*/
struct phys_par
{
	static double A_rb0;
	static double B_rb0;
	static double C_rb0, A_r0, B_r0, C_r0 ;
	phys_par();
	~phys_par();
};
//extern const phys_par ph_p;
extern phys_par& use_ph_p();
extern double& use_a_rb0();
extern double& use_b_rb0();
extern double& use_c_rb0();
extern double& use_a_r0();
extern double& use_b_r0();
extern double& use_c_r0();
/********************************************************************************/
/*				Exact solutions for the five-poin stencil case					*/
extern double omega;
double& use_omega();					
extern double u_static1(double r, double t);
extern double u_static2(double r, double t);
extern double u_static3(double r, double t);
extern double u_travelling(double r, double t);

/*																				*/
/********************************************************************************/