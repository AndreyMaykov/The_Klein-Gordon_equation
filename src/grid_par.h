#pragma once

extern int K, nn_0, nn_max; 
extern const double iteration_precision;
extern const double rb0;
extern const double r0;
extern const double h_x;
extern const double h_t;
extern const char* type_of_scheme;		// Choose the difference scheme type 
extern double fx(double r);
extern double p1(double r);
double& use_h_x();
double& use_h_t();

extern const double CRight0, CRight1, CRight2, CRight3,
					 CLeft0,  CLeft1,  CLeft2,  CLeft3;

extern const double iterations_precision;

using std::valarray;
struct GridData
{
	valarray<double> U2init;
	valarray<double> U1init;
	GridData();
	~GridData();
};