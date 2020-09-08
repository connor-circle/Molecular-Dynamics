#include "all_def.h"
#include "function.h"
void Vanderwaalsforce_term(struct AtomType* ptc, struct Energy* E, struct NonbondingForce* NF, float(*LJ)[4],
	int** TIP3P_Water, int** nonbond_list, int** scale_1_4, float** r)

{
	///////Force routine for MD simulation, Lennard-Jones atoms.
	float dx, dy, dz, dr2, id2, id6;
	float ff;
	float rcut2 = R_CUT * R_CUT;
	float epsilon0, epsilon1, sigma0, sigma1, sigma, epsilon;
	float pair12, pair6;
	int itype, jtype;
	float potentialtemp = 0.;
	int ip, jp;
#pragma omp parallel private(ip,jp,dx, dy, dz, dr2, id2, id6,ff,epsilon0,epsilon1, sigma0,sigma1, sigma, epsilon,pair12,pair6,itype, jtype) \
	shared(nonbond_list,TIP3P_Water,scale_1_4,ptc,r,LJ,NF) reduction(+:potentialtemp) firstprivate(rcut2)
	{
#pragma omp for 
		for (ip = 1; ip <= N_atom; ip++)
		{
			for (jp = ip + 1; jp <= N_atom; jp++)
			{
				if ((nonbond_list[ip][jp] == 0) || (scale_1_4[ip][jp] == 1)) //&&(TIP3P_Water[ip][jp] == 0)) 
				{
					itype = ptc[ip].type;//give the type of atom
					jtype = ptc[jp].type;

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

					dr2 = dx * dx + dy * dy + dz * dz;

					if (dr2 < rcut2)
					{
						sigma1 = LJ[itype][0];
						sigma0 = LJ[jtype][0];
						epsilon0 = LJ[itype][1];
						epsilon1 = LJ[jtype][1];

						if (scale_1_4[ip][jp] == 1)
						{
							epsilon0 = LJ[itype][3];
							sigma0 = LJ[itype][2];

							epsilon1 = LJ[jtype][3];
							sigma1 = LJ[jtype][2];
						}

						sigma = sigma0 + sigma1;
						epsilon = sqrtf(epsilon0 * epsilon1);

						pair12 = powf(sigma, 12.) * epsilon;
						pair6 = powf(sigma, 6.) * epsilon;

						id2 = 1.0 / dr2;
						id6 = id2 * id2 * id2;

						ff = 12.0 * id6 * ((double)pair12 * id6 - pair6) * id2;

						potentialtemp = potentialtemp + id6 * ((double)pair12 * id6 - 2. * pair6);
#pragma omp critical
						{
							NF[ip].F_vander[0] += ff * dx;
							NF[ip].F_vander[1] += ff * dy;
							NF[ip].F_vander[2] += ff * dz;

							NF[jp].F_vander[0] -= ff * dx;
							NF[jp].F_vander[1] -= ff * dy;
							NF[jp].F_vander[2] -= ff * dz;
						}
					}
				}
			}
		}
	}
	E->vander = potentialtemp;
}