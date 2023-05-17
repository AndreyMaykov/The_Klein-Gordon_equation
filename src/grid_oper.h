#pragma once

extern const std::valarray<double> Alpha;
extern const std::valarray<double> Beta;
extern double f_alpha(double);
extern double f_beta(double);
//using std::valarray;

valarray<double>LambdaX(const std::valarray<double> U,
	const std::valarray<double> L0,
	const std::valarray<double> Lpls,
	const std::valarray<double> Lmns
);

class GridOperators
{
public:
	std::valarray<double>Epsilon;

	std::valarray<double>LF00;
	std::valarray<double>LF10;
	valarray<double>LF20;

	std::valarray<double>LF0pls;
	std::valarray<double>LF1pls;
	std::valarray<double>LF2pls;

	std::valarray<double>LF0mns;
	std::valarray<double>LF1mns;
	std::valarray<double>LF2mns;

	/****************************/
	/*		For tests only		*/
	std::valarray<double>TestL0;
	std::valarray<double>TestLpls;
	std::valarray<double>TestLmns;
	/*							*/
	/****************************/

	GridOperators();
	~GridOperators();

};