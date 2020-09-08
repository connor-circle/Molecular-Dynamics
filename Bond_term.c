#include "all_def.h"
#include "function.h"
void bond_stretching_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, int** Bl, float** Bc, float** r)
{
	int ip, jp;
	float rbond, kbond;
	float dx, dy, dz, dij2, dij;
	float dr, ff, F_bond;
	float potential_Bond = 0;
	int ib;

#pragma omp parallel private(ib,ip,jp,rbond,kbond,dx, dy, dz, dij2, dij,dr, ff, F_bond) shared(r,Bl,Bc) \
    reduction(+:potential_Bond)
	{
#pragma omp for schedule(static)
		for (ib = 0; ib < N_Bond; ib++)
		{
			ip = Bl[ib][0];//one bond has two atoms
			jp = Bl[ib][1];

			rbond = Bc[ib][0];
			kbond = Bc[ib][1];

			dx = r[ip][0] - r[jp][0];
			dy = r[ip][1] - r[jp][1];
			dz = r[ip][2] - r[jp][2];

			if (dx > 0.5 * xboxL)
			{
				dx = dx - xboxL;
			}
			else if (dx < -0.5 * xboxL)
			{
				dx = dx + xboxL;
			}
			if (dy > 0.5 * yboxL)
			{
				dy = dy - yboxL;
			}
			else if (dy < -0.5 * yboxL)
			{
				dy = dy + yboxL;
			}
			if (dz > 0.5 * zboxL)
			{
				dz = dz - zboxL;
			}
			else if (dz < -0.5 * zboxL)
			{
				dz = dz + zboxL;
			}

			//there is no jugde for boundary
			dij2 = dx * dx + dy * dy + dz * dz;
			dij = sqrtf(dij2);

			dr = dij - rbond;//the difference between real bond and const bond
			ff = kbond * dr;
			potential_Bond += ff * dr;//U=Kb*(r-r0)^2
			F_bond = ff * (-2.0 / dij);//force 

#pragma omp critical
			{
				BF[ip].F_bond[0] += F_bond * dx;
				BF[ip].F_bond[1] += F_bond * dy;
				BF[ip].F_bond[2] += F_bond * dz;

				BF[jp].F_bond[0] -= F_bond * dx;
				BF[jp].F_bond[1] -= F_bond * dy;
				BF[jp].F_bond[2] -= F_bond * dz;
			}
		}
	}
	E->bond = potential_Bond;
}
