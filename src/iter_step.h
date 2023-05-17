#pragma once

void F_IterStep_Pointer(const valarray<double>* P_U0Prev,
									 valarray<double>* P_U0New,
							   const valarray<double>& U12F,
							   const double& UIFRight,
							   const double& UIFLeft,
							   const valarray<double>& LF00,
							   const valarray<double>& LF0pls,
							   const valarray<double>& LF0mns,
							   const double& CRight1,
							   const double& CLeft1
							  );

bool IterCrit(const valarray<double>& V1,const valarray<double>& V2, const double iteration_precision);