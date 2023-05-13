using std::valarray;

extern const std::valarray<double> Alpha;
extern const std::valarray<double> Beta;
extern double f_alpha(double);
extern double f_beta(double);

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