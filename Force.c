#include "all_def.h"
#include "function.h"
void force(struct BAD_Control* Control, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF, struct NonbondingForce* NF,
	float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst,
	int* CoulList, float* G_k, int** nonbond_list, int** scale_1_4, int** Improperlist, float** Improperconst, int** TIP3P_Water,
	float** r, float** F_all, float** ksin, float** kcos, float** rk)
{
	int ip;
	/*#pragma omp parallel private(ip) shared(BF,NF,F_all)
	{
	#pragma omp for schedule(static)*/
	for (ip = 0; ip <= N_nonwater_atom; ip++)
	{
		BF[ip].F_tor[0] = 0;
		BF[ip].F_tor[1] = 0;
		BF[ip].F_tor[2] = 0;

		BF[ip].F_impro[0] = 0;
		BF[ip].F_impro[1] = 0;
		BF[ip].F_impro[2] = 0;

		BF[ip].F_ub[0] = 0;
		BF[ip].F_ub[1] = 0;
		BF[ip].F_ub[2] = 0;

		BF[ip].F_angle[0] = 0;
		BF[ip].F_angle[1] = 0;
		BF[ip].F_angle[2] = 0;

		BF[ip].F_bond[0] = 0;
		BF[ip].F_bond[1] = 0;
		BF[ip].F_bond[2] = 0;

	}

	for (ip = 0; ip <= N_atom; ip++)
	{
		NF[ip].F_vander[0] = 0;
		NF[ip].F_vander[1] = 0;
		NF[ip].F_vander[2] = 0;

		NF[ip].F_elec[0] = 0;
		NF[ip].F_elec[1] = 0;
		NF[ip].F_elec[2] = 0;

		F_all[ip][0] = 0;
		F_all[ip][1] = 0;
		F_all[ip][2] = 0;
	}

	E->vander = 0;
	E->real = 0;
	E->fourier = 0;
	E->all = 0;
	E->bond = 0;
	E->angle = 0;
	E->torsion = 0;
	E->kineticE = 0;
	E->improper = 0;
	E->ureyBradley = 0;

	torsion_angle_term(Control, E, BF, DihedralList, Dihedralconst, r);

	bond_stretching_term(Control, E, BF, Blist, Bconst, r);

	angle_bending_term(Control, E, BF, Alist, Aconst, r);

	Improper_term(Control, E, BF, Improperlist, Improperconst, r);

	Vanderwaalsforce_term(ptc, E, NF, LJ, TIP3P_Water, nonbond_list, scale_1_4, r);

	Coulombforce_term(ptc, E, NF, CoulList, G_k, nonbond_list, scale_1_4, r, ksin, kcos, rk);

	for (ip = 1; ip <= N_nonwater_atom; ip++)
	{
		F_all[ip][0] = BF[ip].F_angle[0] + BF[ip].F_bond[0] + BF[ip].F_tor[0] + BF[ip].F_impro[0] + BF[ip].F_ub[0];

		F_all[ip][1] = BF[ip].F_angle[1] + BF[ip].F_bond[1] + BF[ip].F_tor[1] + BF[ip].F_impro[1] + BF[ip].F_ub[1];

		F_all[ip][2] = BF[ip].F_angle[2] + BF[ip].F_bond[2] + BF[ip].F_tor[2] + BF[ip].F_impro[2] + BF[ip].F_ub[2];
	}

	for (ip = 1; ip <= N_atom; ip++)
	{
		F_all[ip][0] = F_all[ip][0] + NF[ip].F_vander[0] + NF[ip].F_elec[0];
		F_all[ip][1] = F_all[ip][1] + NF[ip].F_vander[1] + NF[ip].F_elec[1];
		F_all[ip][2] = F_all[ip][2] + NF[ip].F_vander[2] + NF[ip].F_elec[2];
	}

}