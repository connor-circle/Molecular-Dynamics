#include"function.h"
void VWrapALL(float** r)
{
	int ip;

	for (ip = 1; ip <= N_atom; ip++)
	{
		r[ip][0] = r[ip][0] - floor(r[ip][0] / xboxL) * xboxL;

		r[ip][1] = r[ip][1] - floor(r[ip][1] / yboxL) * yboxL;

		r[ip][2] = r[ip][2] - floor(r[ip][2] / zboxL) * zboxL;

	}

}