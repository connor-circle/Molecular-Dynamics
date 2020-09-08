#include "all_def.h"
#include "function.h"
void Verlet_Velocity(struct BAD_Control* Control, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF, struct NonbondingForce* NF,
	float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst,
	int* CoulList, float* G_k, int** nonbond_list, int** Improperlist, float** Improperconst, int*** shakelist, float* shakeconst,
	int** TIP3P_Water, int** scale_1_4, float** r, float** rold, float** v, float** F_all, float** ksin, float** kcos, float** rk)
{
	float tem;
	float hdt = 0.5 * step_t;
	int ip;

	for (ip = 1; ip <= N_atom; ip++)
	{
		rold[ip][0] = r[ip][0];
		rold[ip][1] = r[ip][1];
		rold[ip][2] = r[ip][2];

		tem = hdt / ptc[ip].m;

		v[ip][0] = v[ip][0] + tem * F_all[ip][0];
		v[ip][1] = v[ip][1] + tem * F_all[ip][1];
		v[ip][2] = v[ip][2] + tem * F_all[ip][2];

		r[ip][0] = r[ip][0] + v[ip][0] * step_t;
		r[ip][1] = r[ip][1] + v[ip][1] * step_t;
		r[ip][2] = r[ip][2] + v[ip][2] * step_t;
	}

	force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
		Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all, ksin, kcos, rk);
	//new velocity in next step

	for (ip = 1; ip <= N_atom; ip++)
	{
		tem = hdt / ptc[ip].m;

		v[ip][0] = v[ip][0] + tem * F_all[ip][0];
		v[ip][1] = v[ip][1] + tem * F_all[ip][1];
		v[ip][2] = v[ip][2] + tem * F_all[ip][2];
	}


}