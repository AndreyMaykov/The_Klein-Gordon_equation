#pragma once

struct scheme_type_error
{
	int ii;
	scheme_type_error(int iii){ii = iii;}
};

double weight_matrix[4][3];
void CreateWeightMatrix(char* type_of_scheme);
