#include "all_def.h"
#include "function.h"
void steepest_descent(struct BAD_Control* Control, struct PDB_format* PDB, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF,
	struct NonbondingForce* NF, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList,
	float** Dihedralconst, int* CoulList, float* G_k, int*** shakelist, float* shakeconst, int** nonbond_list,
	int** Improperlist, float** Improperconst, int** TIP3P_Water, int** scale_1_4, float** r, float** v, float** rold, float** F_all,
	float** ksin, float** kcos, float** rk)
{
	int  k = 0;
	int count_model = 1;
	float epsilon = 1e-5;
	float beta = 0;
	float ¦Á;
	float E_K1 = 0, E_K = 0;

	FILE* steepestfile;
	if ((steepestfile = fopen("../../Output/SteepestCordinate-A.pdb", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
		Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all, ksin, kcos, rk);

	E->all = E->real + E->fourier + E->vander + E->bond + E->angle + E->torsion + E->ureyBradley + E->improper - E->self_K_SPACE_SUM;
	E_K = 0;
	E_K1 = E->all;
	printf("%f\n", E_K1);

	while (k <= STEEPEST_STEP)
	{
		if (fabs((double)E_K - E_K1) < epsilon)
		{
			break;
		}
		if (k == 0)
		{
			for (int im = 1; im <= N_atom; im++)
			{
				rold[im][0] = r[im][0];
				rold[im][1] = r[im][1];
				rold[im][2] = r[im][2];

				beta = sqrt((double)F_all[im][0] * F_all[im][0] + (double)F_all[im][1] * F_all[im][1] + (double)F_all[im][2] * F_all[im][2]);

				F_all[im][0] = F_all[im][0] / beta;
				F_all[im][1] = F_all[im][1] / beta;
				F_all[im][2] = F_all[im][2] / beta;
			}

			¦Á = 0.008;

		}
		else
		{

			if ((E_K1 - E_K) > 0)
			{
				for (int im = 1; im <= N_atom; im++)
				{
					r[im][0] = rold[im][0];
					r[im][1] = rold[im][1];
					r[im][2] = rold[im][2];
				}

				force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
					Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all,
					ksin, kcos, rk);

				for (int im = 1; im <= N_atom; im++)
				{
					beta = sqrt((double)F_all[im][0] * F_all[im][0] + (double)F_all[im][1] * F_all[im][1] + (double)F_all[im][2] * F_all[im][2]);

					F_all[im][0] = F_all[im][0] / beta;
					F_all[im][1] = F_all[im][1] / beta;
					F_all[im][2] = F_all[im][2] / beta;
				}

				¦Á = ¦Á * 0.5;
			}
			else if ((E_K1 - E_K) < 0)
			{
				for (int im = 1; im <= N_atom; im++)
				{
					rold[im][0] = r[im][0];
					rold[im][1] = r[im][1];
					rold[im][2] = r[im][2];
				}

				for (int im = 1; im <= N_atom; im++)
				{
					beta = sqrt((double)F_all[im][0] * F_all[im][0] + (double)F_all[im][1] * F_all[im][1] + (double)F_all[im][2] * F_all[im][2]);

					F_all[im][0] = F_all[im][0] / beta;
					F_all[im][1] = F_all[im][1] / beta;
					F_all[im][2] = F_all[im][2] / beta;
				}

				¦Á = ¦Á * 1.2;
			}
		}
		for (int im = 1; im <= N_atom; im++)
		{
			r[im][0] = rold[im][0] + ¦Á * F_all[im][0];
			r[im][1] = rold[im][1] + ¦Á * F_all[im][1];
			r[im][2] = rold[im][2] + ¦Á * F_all[im][2];

			//printf("%f,%f,%f\n", ptc[im].r[0], ptc[im].r[1], ptc[im].r[2]);
		}

		printf("k and  ¦Á  %d  %f\n", k, ¦Á);
		k++;
		E_K = E_K1;

		//RattleA(ptc,r, rold,v,shakelist, shakeconst);

		VWrapALL(r);

		force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
			Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all,
			ksin, kcos, rk);
		E->all = E->real + E->fourier + E->vander + E->bond + E->angle + E->torsion + E->ureyBradley + E->improper - E->self_K_SPACE_SUM;
		E_K1 = E->all;

		for (int ip = 1; ip <= N_atom; ip++)
		{
			r[ip][0] = r[ip][0] * (SIGMA);
			r[ip][1] = r[ip][1] * (SIGMA);
			r[ip][2] = r[ip][2] * (SIGMA);
		}

		fprintf(steepestfile, "%s%8s%d\n", "MODEL", " ", count_model);
		for (int ip = 1; ip <= N_atom; ip++)
		{
			fprintf(steepestfile, "%-6s%5d%1s%-4s%1s%3s%1s%1s%4d%4s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%-4s%2s\n", "ATOM", PDB[ip].serial, " ",
				PDB[ip].name, " ", PDB[ip].resName, " ", PDB[ip].chainID, PDB[ip].resSeq, " ", r[ip][0], r[ip][1],
				r[ip][2], 1.00, 0.00, " ", PDB[ip].segID, PDB[ip].element);
		}
		fprintf(steepestfile, "%s%6s%d%6s%s %s   %d\n", "TER", " ", N_atom + 1, " ", "MET", PDB[N_atom].chainID, PDB[N_atom].resSeq);
		fprintf(steepestfile, "%s\n", "ENDMDL");
		count_model++;

		for (int ip = 1; ip <= N_atom; ip++)
		{
			r[ip][0] = r[ip][0] / (SIGMA);
			r[ip][1] = r[ip][1] / (SIGMA);
			r[ip][2] = r[ip][2] / (SIGMA);
		}

		printf("P_vander is %f\n", E->vander);
		printf("P_fourier is %f\n", E->fourier);
		printf("P_real is %f\n", E->real);
		printf("P_angle is %f\n", E->angle);
		printf("P_bond is %f\n", E->bond);
		printf("P_all is %f\n", E->all);
		printf("\n");
	}
	if (fclose(steepestfile))
	{
		printf("can not close the file\n");
		exit(0);
	}
	printf("the end of steepest\n");
	printf("\n");

}

void steepest_descentB(struct BAD_Control* Control, struct PDB_format* PDB, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF,
	struct NonbondingForce* NF, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList,
	float** Dihedralconst, int* CoulList, float* G_k, int*** shakelist, float* shakeconst, int** nonbond_list,
	int** Improperlist, float** Improperconst, int** TIP3P_Water, int** scale_1_4, float** r, float** v, float** rold, float** F_all,
	float** ksin, float** kcos, float** rk)
{
	int  k = 0;
	int count_model = 1;
	float epsilon = 1e-5;
	float beta = 0;
	float ¦Á;
	float E_K1 = 0, E_K = 0;

	FILE* steepestfile;
	if ((steepestfile = fopen("../../Output/SteepestCordinate-B.pdb", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
		Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all, ksin, kcos, rk);

	E->all = E->real + E->fourier + E->vander + E->bond + E->angle + E->torsion + E->ureyBradley + E->improper - E->self_K_SPACE_SUM;
	E_K = 0;
	E_K1 = E->all;
	printf("%f\n", E_K1);

	while (k <= STEEPEST_STEP)
	{
		if (fabs((double)E_K - E_K1) < epsilon)
		{
			break;
		}
		if (k == 0)
		{
			for (int im = 1; im <= N_nonwater_atom; im++)
			{
				rold[im][0] = r[im][0];
				rold[im][1] = r[im][1];
				rold[im][2] = r[im][2];

				beta = sqrt((double)F_all[im][0] * F_all[im][0] + (double)F_all[im][1] * F_all[im][1] + (double)F_all[im][2] * F_all[im][2]);

				F_all[im][0] = F_all[im][0] / beta;
				F_all[im][1] = F_all[im][1] / beta;
				F_all[im][2] = F_all[im][2] / beta;
			}

			¦Á = 0.008;

		}
		else
		{

			if ((E_K1 - E_K) > 0)
			{
				for (int im = 1; im <= N_nonwater_atom; im++)
				{
					r[im][0] = rold[im][0];
					r[im][1] = rold[im][1];
					r[im][2] = rold[im][2];
				}

				force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
					Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r,
					F_all, ksin, kcos, rk);

				for (int im = 1; im <= N_nonwater_atom; im++)
				{
					beta = sqrt((double)F_all[im][0] * F_all[im][0] + (double)F_all[im][1] * F_all[im][1] + (double)F_all[im][2] * F_all[im][2]);

					F_all[im][0] = F_all[im][0] / beta;
					F_all[im][1] = F_all[im][1] / beta;
					F_all[im][2] = F_all[im][2] / beta;
				}

				¦Á = ¦Á * 0.5;
			}
			else if ((E_K1 - E_K) < 0)
			{
				for (int im = 1; im <= N_nonwater_atom; im++)
				{
					rold[im][0] = r[im][0];
					rold[im][1] = r[im][1];
					rold[im][2] = r[im][2];
				}

				for (int im = 1; im <= N_nonwater_atom; im++)
				{
					beta = sqrtf((double)F_all[im][0] * F_all[im][0] + (double)F_all[im][1] * F_all[im][1] + (double)F_all[im][2] * F_all[im][2]);

					F_all[im][0] = F_all[im][0] / beta;
					F_all[im][1] = F_all[im][1] / beta;
					F_all[im][2] = F_all[im][2] / beta;
				}

				¦Á = ¦Á * 1.2;
			}
		}
		for (int im = 1; im <= N_nonwater_atom; im++)
		{
			r[im][0] = rold[im][0] + ¦Á * F_all[im][0];
			r[im][1] = rold[im][1] + ¦Á * F_all[im][1];
			r[im][2] = rold[im][2] + ¦Á * F_all[im][2];

			//printf("%f,%f,%f\n", ptc[im].r[0], ptc[im].r[1], ptc[im].r[2]);
		}

		printf("k and  ¦Á  %d  %f\n", k, ¦Á);
		k++;
		E_K = E_K1;

		//RattleA(ptc,r, rold,v,shakelist, shakeconst);

		VWrapALL(r);

		force(Control, ptc, E, BF, NF, LJ, Blist, Bconst, Alist, Aconst, DihedralList,
			Dihedralconst, CoulList, G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all,
			ksin, kcos, rk);

		E->all = E->real + E->fourier + E->vander + E->bond + E->angle + E->torsion + E->ureyBradley + E->improper - E->self_K_SPACE_SUM;

		E_K1 = E->all;

		for (int ip = 1; ip <= N_nonwater_atom; ip++)
		{
			r[ip][0] = r[ip][0] * (SIGMA);
			r[ip][1] = r[ip][1] * (SIGMA);
			r[ip][2] = r[ip][2] * (SIGMA);
		}

		fprintf(steepestfile, "%s%8s%d\n", "MODEL", " ", count_model);
		for (int ip = 1; ip <= N_atom; ip++)
		{
			fprintf(steepestfile, "%-6s%5d%1s%-4s%1s%3s%1s%1s%4d%4s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%-4s%2s\n", "ATOM", PDB[ip].serial, " ",
				PDB[ip].name, " ", PDB[ip].resName, " ", PDB[ip].chainID, PDB[ip].resSeq, " ", r[ip][0], r[ip][1],
				r[ip][2], 1.00, 0.00, " ", PDB[ip].segID, PDB[ip].element);
		}
		fprintf(steepestfile, "%s%6s%d%6s%s %s   %d\n", "TER", " ", N_atom + 1, " ", "MET", PDB[N_atom].chainID, PDB[N_atom].resSeq);
		fprintf(steepestfile, "%s\n", "ENDMDL");
		count_model++;

		for (int ip = 1; ip <= N_nonwater_atom; ip++)
		{
			r[ip][0] = r[ip][0] / (SIGMA);
			r[ip][1] = r[ip][1] / (SIGMA);
			r[ip][2] = r[ip][2] / (SIGMA);
		}

		printf("P_vander is %f\n", E->vander);
		printf("P_fourier is %f\n", E->fourier);
		printf("P_real is %f\n", E->real);
		printf("P_angle is %f\n", E->angle);
		printf("P_bond is %f\n", E->bond);
		printf("P_all is %f\n", E->all);
		printf("\n");
	}
	if (fclose(steepestfile))
	{
		printf("can not close the file\n");
		exit(0);
	}
	printf("the end of steepest\n");
	printf("\n");
}