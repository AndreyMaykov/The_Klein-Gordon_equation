#include <cmath>
#include <iostream.h>
#include <valarray>		//Нужно только для теста с NN != 1;

const double Pi = 3.141592653589793;

extern const int NN = 3;
extern const double gamma = 3.2;/////////!!!!
extern const double rb0 = 10.0;
extern const double r0 = 21.0;
double rho(double r){return 1;};
double kappa(double r){return 1;};
double sigma(double r){return 0;};
double q(double r){return 0;};
//double p1(double r){return 1;};
/************************************************************/
/*	  "Расширения" функций, описывающих среду, за r = r0	*/

double net_rho(double r){return (r<=r0)?rho(r):1;};
double net_kappa(double r){return (r<=r0)?kappa(r):1;};
double net_sigma(double r){return (r<=r0)?sigma(r):0;};
double net_q(double r){return (r<=r0)?q(r):0;};


/********************************************************************/
/*	Коэффициенты дифф. оператора левой части граничных условий:		*/
/*	(a*du/dt + b*du/dr + c*u)|_(r = rb0 || r0) = ... + f_(rb0 || r0)*/
/*	а также неоднородности в граничных условиях и в самом дифф.		*/
/*	уравнении  -- см. в конце файла									*/


														
/********************************************************************/

/********************************************************************************/
/*							Данные для тестов									*/
/*			              Тест 1: стоячие волны									*/

int m = 2;
double omega = Pi * m /(r0 - rb0);
char* ff;
double& use_omega()								
{													//Почему-то внутри use_omega получаем
//	cout <<"omega in use_omega = "<<omega<<'\n';	//omega = 0, поэтому важно использовать именно
	static double om = Pi * m /(r0 - rb0);			// <-- такое выражение. Напрашивается такое
//	cout <<"om in use_omega = "<<om<<'\n';			//объяснение: то обращение к use_omega, из 
//	cin >> ff;										//которого происходит вывод omega, имеет место	
	return om;										//в Net_par.cpp и там правильное значение может
};													//быть только у константного выражения
													//Pi * m /(r0 - rb0); обращения же к u_static
													//имеют место в main.cpp, а туда правильно пере-
													//дается и само omega.
const double oo = use_omega();
//double h_t = 1.0; double h_x = 2.0;			//Не забыть изменить значения K и/или fc в Net_par.cpp
double h_t = 0.6; double h_x = 2.0;
/*	КОНСЕРВАТИВНАЯ СХЕМА:
	omega_tild находится с помощью  MATHEMATICA  ConsTest.nb  из уравнения
 
	tg(omega_tild * h_t / 2) = (h_t/h_x)*sin(omega * h_x/2)

	при m = 1, r0 = 20, rb0 = 10, h_x = 2 (т. е. K = 5), h_t = 1 (т. е. fc = 0.5)
			omega_tild = 0.306592585897309;
	при m = 1, r0 = 20, rb0 = 10, h_x = 2 (т. е. K = 5), h_t = 0.6 (т. е. fc = 0.3)
			omega_tild = 0.3081362764762855
	при m = 2, r0 = 20, rb0 = 10, h_x = 2 (т. е. K = 5), h_t = 0.6 (т. е. fc = 0.3)
			omega_tild = 0.5818042033204008
*/
//double omega_tild = 0.306592585897309;
//double omega_tild = 0.3081362764762855;
//double omega_tild = 0.5818042033204008;

/*	СХЕМА "КРЕСТ":
	omega_tild находится с помощью  MATHEMATICA  CrosTest.nb  из уравнения
 
	sin(omega_tild * h_t / 2) = (h_t/h_x)*sin(omega * h_x/2)

	при m = 1, r0 = 20, rb0 = 10, h_x = 2 (т. е. K = 5), h_t = 1 (т. е. fc = 0.5)
			omega_tild = ;
	при m = 1, r0 = 20, rb0 = 10, h_x = 2 (т. е. K = 5), h_t = 0.6 (т. е. fc = 0.3)
			omega_tild = ;
	при m = 2, r0 = 20, rb0 = 10, h_x = 2 (т. е. K = 5), h_t = 0.6 (т. е. fc = 0.3)
			omega_tild = 0.590874802986392
*/




double omega_tild =  0.676914121276099;


//		Точные решения:															
extern double u_static1(double r, double t){return sin(omega_tild * t) * sin(omega *(r - rb0));};
extern double u_static2(double r, double t){return cos(omega_tild * t) * sin(omega *(r - rb0));};
extern double u_static3(double r, double t){return sin(omega_tild * t) * cos(omega *(r - rb0));};
extern double u_travelling(double r, double t){return
										 sin( omega_tild * t - omega * (r-rb0) );};

/********************************************************************************/
/*						Начальные условия для  u_static1						*/
//				не забыть:	a_rb0= 0; b_rb0 = 0; c_rb0 = 1;
//							a_r0 = 0; b_r0 = 0; c_r0 = 1; 
/*
extern double u_0_init(double r){return 0;};
extern double u_1_init(double r)
	{return sin(omega_tild* h_t/2)/(h_t/2) *sin(omega*(r-rb0));};
*/
/********************************************************************************/
/*						Начальные условия для  u_static2						*/
//				не забыть:	a_rb0= 0; b_rb0 = 0; c_rb0 = 1;
//							a_r0 = 0; b_r0 = 0; c_r0 = 1; 
/*
/*
extern double u_0_init(double r)
	{return cos(omega_tild* h_t/2) *sin(omega*(r-rb0));};
extern double u_1_init(double r){return 0;};
*/
/********************************************************************************/
/********************************************************************************/
/*						Начальные условия для  u_static3						*/
//				не забыть:	a_rb0= 0; b_rb0 = 1; c_rb0 = 0;
//							a_r0 = 0; b_r0 = 1; c_r0 = 0; 

extern double u_0_init(double r){return 0;};
extern double u_1_init(double r)
	{ 
	double ff;
	if((fabs(r - rb0) 
		> 0.25*h_t)&
		(fabs(r0 - r) > 
		0.25*h_t))
			ff = sin(omega_tild* h_t/2.0)/(h_t/2.0) *cos(omega*(r-rb0));
	else
			ff = sin(omega_tild* h_t/2.0)/(h_t/2.0) *cos(omega*(r-rb0))*cos(omega*h_x/2)
	;
		return ff;
	};

/********************************************************************************/
/*			              Тест 2: Бегущие волны									*/

//	Точное решение: u(r, t) = sin( omega_tild * t - omega * (r-rb0) ) 
//	Начальные условия:
/*
extern double u_0_init(double r){return (0.0 - cos(omega_tild*h_t/2)*sin(omega*(r-rb0)));};
extern double u_1_init(double r)
	{ 
	double ff;
	if((fabs(r - rb0) 
		> 0.25*h_t)&
		(fabs(r0 - r) > 
		0.25*h_t))
			ff = sin(omega_tild* h_t/2.0)/(h_t/2.0) *cos(omega*(r-rb0));
	else
			ff = sin(omega_tild* h_t/2.0)/(h_t/2.0) *cos(omega*(r-rb0))*cos(omega*h_x/2)
	;
		return ff;
	};
*/
//Граничные условия:
//Вариант 1: a_r0 = b_r0 = a_rb0 = b_rb0 = 0; c_r0 = c_rb0 = 1;
/*
extern double f_rb0(double t)
{
	return 0.25*( 
		sin( omega_tild*(t + h_t/2) - omega * ( - h_x/2) ) +
		sin( omega_tild*(t + h_t/2) - omega * (   h_x/2) ) +
		sin( omega_tild*(t - h_t/2) - omega * ( - h_x/2) ) +
		sin( omega_tild*(t - h_t/2) - omega * (   h_x/2) )
		);
};
extern double f_r0(double t)
{
	return 0.25*(
		sin( omega_tild*(t + h_t/2) - omega * ( r0 - rb0 - h_x/2) ) +
		sin( omega_tild*(t + h_t/2) - omega * ( r0 - rb0 + h_x/2) ) +
		sin( omega_tild*(t - h_t/2) - omega * ( r0 - rb0 - h_x/2) ) +
		sin( omega_tild*(t - h_t/2) - omega * ( r0 - rb0 + h_x/2) )
		);
};
*/

/*Вариант 2: 

консервативная схема:
 
a_rb0 = 1;
b_rb0 = cos(use_omega()*h_x/2);
c_rb0 = 0;

a_r0 = 1;
b_r0 = cos(use_omega()*h_x/2);
c_r0 = 0;

схема крест:

a_rb0 = 1;
b_rb0 = cos(use_omega()*h_x/2)/cos(omega_tild*h_t/2);
c_rb0 = 0;

a_r0 = 1;
b_r0 = cos(use_omega()*h_x/2)/cos(omega_tild*h_t/2);
c_r0 = 0;
*/

/*						Окончание данных для тестов								*/
/********************************************************************************/
extern double f_internal(double r, double t){return 0;};
double f_rb0(double t)
{
	return 0;
};
double f_r0(double t)
{
	return 0;
};
/********************************************************************************/
/*
struct phys_par
{
	double A_rb0;
	double B_rb0;
	double C_rb0, A_r0, B_r0, C_r0 ;
	phys_par();
	~phys_par();
};

phys_par::phys_par()
{
	cout << "      constructor for phys_par"<<'\n';
	A_rb0 = 1; 
	B_rb0 = b_rb0; 
	C_rb0 = c_rb0;
	A_r0 = a_r0; B_r0 = b_r0; C_r0 = c_r0;
	cout << A_rb0 <<'\t'<< B_rb0<<'\t'<< C_rb0<<'\t'<< A_r0 <<'\t'<< B_r0<<'\t'<< C_r0<<'\n' ;
	char* ch = "ch";
	cin>>ch;
};
phys_par::~phys_par(){cout << "destructor for phys_par"<<'\n';};
*/
/*
struct phys_par
{
	static double A_rb0;
	static double B_rb0;
	static double C_rb0, A_r0, B_r0, C_r0 ;
	phys_par();
	~phys_par();
};

phys_par::phys_par()
{
	cout << "      constructor for phys_par"<<'\n';
	cout << A_rb0 <<'\t'<< B_rb0<<'\t'<< C_rb0<<'\t'<< A_r0 <<'\t'<< B_r0<<'\t'<< C_r0<<'\n' ;
//	char* ch = "ch";
//	cin>>ch;
};

phys_par::~phys_par(){cout << "destructor for phys_par"<<'\n';};
double phys_par::A_rb0 = 1.0;
double phys_par::B_rb0 = cos(use_omega()*h_x/2);
double phys_par::C_rb0 = 0.0;
double phys_par::A_r0 = 1.0;
double phys_par::B_r0 = cos(use_omega()*h_x/2);
double phys_par::C_r0 = 0.0;

/*extern const*///phys_par ph_p0 = phys_par();
/*
extern phys_par& use_ph_p()
{	
	static phys_par jj = phys_par();
	return jj;
};
phys_par ph_p0 = use_ph_p();
*/
extern double& use_a_rb0()
{	
	static double jj = 0.0;
	return jj;
};
extern double& use_b_rb0()
{	
//	static double jj = cos(use_omega()*h_x/2);	//Для консервативной схемы
//	static double jj = cos(use_omega()*h_x/2)/cos(omega_tild*h_t/2);	//Для схемы "крест"
	static double jj = -1.0;
	return jj;
};
extern double& use_c_rb0()
{	
	static double jj = 0.0;
	return jj;
};
extern double& use_a_r0()
{	
	static double jj = 0.0;
	return jj;
};
//char* fff;
extern double& use_b_r0()
{
//	cout << ::omega <<'\n';
//	cout << b_r0<<'\n';
//	static double jj = cos(use_omega()*h_x/2);	//Для консервативной схемы
//	static double jj = cos(use_omega()*h_x/2)/cos(omega_tild*h_t/2);	//Для схемы "крест"
	static double jj = 1.0;
//	cout << "jj in use_b_r0 ="<<jj<<'\n';
//	cin>>fff;
	return jj;
};
extern double& use_c_r0()
{	
	static double jj = 0.0;
	return jj;
};

/********************************************************************************/