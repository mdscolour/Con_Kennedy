#include "point.h"

Matrix Matrix::inv()
{
	double determinant = row1.x*(row2.y*row3.z - row3.y*row2.z)
		- row1.y*(row2.x*row3.z - row2.z*row3.x)
		+ row1.z*(row2.x*row3.y - row2.y*row3.x);
	double invdet = 1 / determinant;
	return Matrix(
		(row2.y*row3.z - row3.y*row2.z)*invdet,
		-(row1.y*row3.z - row1.z*row3.y)*invdet,
		(row1.y*row2.z - row1.z*row2.y)*invdet,
		-(row2.x*row3.z - row2.z*row3.x)*invdet,
		(row1.x*row3.z - row1.z*row3.x)*invdet,
		-(row1.x*row2.z - row2.x*row1.z)*invdet,
		(row2.x*row3.y - row3.x*row2.y)*invdet,
		-(row1.x*row3.y - row3.x*row1.y)*invdet,
		(row1.x*row2.y - row2.x*row1.y)*invdet);
}