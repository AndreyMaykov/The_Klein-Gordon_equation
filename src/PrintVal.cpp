/************************************************************************************/
// See the original

#include <iostream>
#include <valarray>

#include "PrintVal.h";

using std::cout;


void print_valarray_length(const std::valarray<double>& v, int length =10)
{
	int size_v =  v.size();
	int nn = 0; int nn10 = length; int nn_line = 1;
//	cout <<"\nlength = "<< length<<'\n';
	cout <<"\nline "<< nn_line<<": ";
	while((nn < size_v))
	{
		if (/*(cout<<"  nn = "<<nn<<"; "<<"  nn10 = "<<nn10<<"; ")&&*/(nn == nn10))
		{
			cout<< "\n\n";
			nn_line++;
			nn10 = nn10+length;
			cout <<"line "<< nn_line<<": ";
		};
		cout <<"  "<< v[nn];
		
		nn++;
	};
	cout << '\n';
}	
