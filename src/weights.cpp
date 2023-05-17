#include <iostream>
#include <valarray>

using std::cout;
using std::cin;
using std::valarray;

struct scheme_type_error
{
	int ii;
	scheme_type_error(int iii){ii = iii;}
};

double weight_matrix[4][3];
void CreateWeightMatrix(const char* type_of_scheme)
{
	cout <<"\ntype of scheme is \""<<type_of_scheme<<"\"\n\n" ;
	int tos = 2;
	if (type_of_scheme == "cross")
		tos = 1;
	else
		if(type_of_scheme == "conservative")
			tos = 2;
		else
			if(type_of_scheme == "unknown")
				tos = 0;
			else
			{
//				cout << "\n\n   error in definition: type_of_scheme = " <<type_of_scheme <<"\n\n";
				throw scheme_type_error(tos);
			};

	switch (tos)
	{
	case 1: 
		{
			weight_matrix[0][0]=0;			weight_matrix[0][2]=0;
			weight_matrix [1][0]=1.0;		weight_matrix [1][2]= - 1.0;
			weight_matrix [2][0]=0;			weight_matrix [2][2]=0;
			weight_matrix [3][0]=0;			weight_matrix [3][2]=0;
		};
		break;
	case 2:
		{
			weight_matrix[0][0]=1.0/4;		weight_matrix[0][2]=1.0/4;
			weight_matrix [1][0]=1.0;		weight_matrix [1][2]= - 1.0;
			weight_matrix [2][0]=1.0/4;		weight_matrix [2][2]=1.0/4;
			weight_matrix [3][0]=1.0/4;		weight_matrix [3][2]=1.0/4;
		};
		break;
	case 0:
		cout<<"\n\n   unknown\n\n";
		{
			cout << "weight_matrix[0][0] = ";	cin >> weight_matrix[0][0]; cout << '\t';
			cout << "weight_matrix[0][2] = ";	cin >> weight_matrix[0][2]; cout << '\n';
			cout << "weight_matrix[1][0] = ";	cin >> weight_matrix[1][0]; cout << '\t';
			cout << "weight_matrix[1][2] = ";	cin >> weight_matrix[1][2]; cout << '\n';
			cout << "weight_matrix[2][0] = ";	cin >> weight_matrix[2][0]; cout << '\t';
			cout << "weight_matrix[2][2] = ";	cin >> weight_matrix[2][2]; cout << '\n';
			cout << "weight_matrix[3][0] = ";	cin >> weight_matrix[3][0]; cout << '\t';
			cout << "weight_matrix[3][2] = ";	cin >> weight_matrix[3][2]; cout << '\n';
		};
		break;
	};
for (int ll = 0; ll < 4; ll++)
{
	double dd;
	switch (ll)
	{
	case 1: dd = 0; break;
	default: dd = 1;
	};
	weight_matrix[ll][1] = dd - weight_matrix [ll][0] - weight_matrix[ll][2];
};

/************************************************************************************/
/*								See the original									*/
/*
	for(ll =0; ll<4;ll++)
	{ 
		for(int mm = 0; mm<3; mm++)
		{
			cout << weight_matrix[ll][mm] << "\t" ;
		};
		cout << "\n\n";
	};
*/
/*																					*/
/************************************************************************************/
}

