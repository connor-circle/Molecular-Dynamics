#include "all_def.h"
#include "function.h"
int main()
{

	/*                                      Normol Array                                           */

	float shakeconst[3];/*the 2 are energy const and distance const*/

	float LJ_set[LJ_list][4];/*resent different types selecting from 56 atoms ; the 2 are energy const and distance const*/

	int count_model = 1;//control the output step

	/*                                 Velocity. Position. Sum of Force                              */

	float** r;

	if ((r = (float**)malloc((N_atom + 1) * sizeof(float*))) == NULL) return NULL;

	if ((r[0] = (float*)malloc((N_atom + 1) * 3 * sizeof(float))) == NULL) return NULL;

	printf("able to allocate memeory for position\n");

	for (int i = 1; i < N_atom + 1; i++)
	{
		r[i] = r[i - 1] + 3;
	}

	float** rold;

	if ((rold = (float**)malloc((N_atom + 1) * sizeof(float*))) == NULL) return NULL;

	if ((rold[0] = (float*)malloc((N_atom + 1) * 3 * sizeof(float))) == NULL) return NULL;

	printf("able to allocate memeory for old_position\n");

	for (int i = 1; i < N_atom + 1; i++)
	{
		rold[i] = rold[i - 1] + 3;
	}

	float** F_all;

	if ((F_all = (float**)malloc((N_atom + 1) * sizeof(float*))) == NULL) return NULL;

	if ((F_all[0] = (float*)malloc((N_atom + 1) * 3 * sizeof(float))) == NULL) return NULL;

	printf("able to allocate memeory for F_all\n");

	for (int i = 1; i < N_atom + 1; i++)
	{
		F_all[i] = F_all[i - 1] + 3;
	}

	float** v;

	if ((v = (float**)malloc((N_atom + 1) * sizeof(float*))) == NULL) return NULL;

	if ((v[0] = (float*)malloc((N_atom + 1) * 3 * sizeof(float))) == NULL) return NULL;

	printf("able to allocate memeory for F_all\n");

	for (int i = 1; i < N_atom + 1; i++)
	{
		v[i] = v[i - 1] + 3;
	}

	float** rk;
	if ((rk = (float**)malloc(KSQ_num * sizeof(float*))) == NULL) return NULL;

	if ((rk[0] = (float*)malloc(KSQ_num * 4 * sizeof(float))) == NULL) return NULL;

	printf("able to allocate memeory for rk\n");

	for (int i = 1; i < KSQ_num; i++)
	{
		rk[i] = rk[i - 1] + 4;
	}

	/*                                  Structure Dynamic Array Initialization                           */

	struct AtomType* part = (struct AtomType*)calloc((N_atom + 1), sizeof(struct AtomType));

	struct NonbondingForce* NF = (struct  NonbondingForce*)calloc((N_atom + 1), sizeof(struct  NonbondingForce));

	struct BondingForce* BF = (struct  BondingForce*)calloc((N_nonwater_atom + 1), sizeof(struct  BondingForce));

	struct PDB_format* PDB = (struct PDB_format*)calloc((N_atom + 1), sizeof(struct PDB_format));

	struct PSFINPUT* PSF = (struct PSFINPUT*)calloc((N_atom + 1), sizeof(struct PSFINPUT));

	struct BAD_Control* Control = (struct BAD_Control*)calloc(1, sizeof(struct PSFINPUT));



	/*                                      Array Initialization                                   */
	float* G_k;

	G_k = (float*)malloc(125 * sizeof(float));

	int* CoulList;

	CoulList = (int*)malloc(Coul_List * sizeof(int));

	/*==============================================================================================*/

	int** nonbond_list;/* nonbond_list[N_atom + 1][N_atom + 1]*/

	if ((nonbond_list = (int**)malloc((N_atom + 1) * sizeof(int*))) == NULL) return NULL;

	if ((nonbond_list[0] = (int*)malloc((N_atom + 1) * (N_atom + 1) * sizeof(int))) == NULL) return NULL;

	printf("able to allocate memeory for nonbond_list\n");

	for (int i = 1; i < N_atom + 1; i++)
	{
		nonbond_list[i] = nonbond_list[i - 1] + (N_atom + 1);
	}

	int** scale_1_4;/*int scale_1_4[N_atom + 1][N_atom + 1]*/

	if ((scale_1_4 = (int**)malloc((N_atom + 1) * sizeof(int*))) == NULL) return NULL;

	if ((scale_1_4[0] = (int*)malloc((N_atom + 1) * (N_atom + 1) * sizeof(int))) == NULL) return NULL;

	printf("able to allocate memeory for scale_1_4\n");

	for (int i = 1; i < N_atom + 1; i++)
	{
		scale_1_4[i] = scale_1_4[i - 1] + (N_atom + 1);
	}

	int** waterlist;/*waterlist[N_water_Molecular][3]*/

	if ((waterlist = (int**)malloc(N_water_Molecular * sizeof(int*))) == NULL) return NULL;

	if ((waterlist[0] = (int*)malloc(N_water_Molecular * 3 * sizeof(int))) == NULL) return NULL;

	printf("able to allocate memeory for waterlist\n");

	for (size_t i = 1; i < N_water_Molecular; i++)
	{
		waterlist[i] = waterlist[i - 1] + 3;
	}

	int** TIP3P_Water;/*int TIP3P_Water[N_atom + 1][N_atom + 1]*/

	if ((TIP3P_Water = (int**)malloc((N_atom + 1) * sizeof(int*))) == NULL) return NULL;

	if ((TIP3P_Water[0] = (int*)malloc((N_atom + 1) * (N_atom + 1) * sizeof(int))) == NULL) return NULL;

	printf("able to allocate memeory for TIP3P_Water\n");

	for (size_t i = 1; i < N_atom + 1; i++)
	{
		TIP3P_Water[i] = TIP3P_Water[i - 1] + (N_atom + 1);
	}

	int** Blist;/*Blist[N_Bond][2]*/

	if ((Blist = (int**)malloc(N_Bond * sizeof(int*))) == NULL) return NULL;

	if ((Blist[0] = (int*)malloc(N_Bond * 2 * sizeof(int))) == NULL) return NULL;

	printf("able to allocate memeory for Blist\n");

	for (size_t i = 1; i < N_Bond; i++)
	{
		Blist[i] = Blist[i - 1] + 2;
	}

	float** Bconstraintlist;/* Bconstraintlist[N_Bond][2] 0 energy const and bond 1const*/

	if ((Bconstraintlist = (float**)malloc(N_Bond * sizeof(float*))) == NULL) return NULL;

	if ((Bconstraintlist[0] = (float*)malloc(N_Bond * 2 * sizeof(float))) == NULL) return NULL;

	for (size_t i = 1; i < N_Bond; i++)
	{
		Bconstraintlist[i] = Bconstraintlist[i - 1] + 2;
	}

	int** Alist;/*Alist[N_Angle][3]*/

	if ((Alist = (int**)malloc(N_Angle * sizeof(int*))) == NULL) return NULL;

	if ((Alist[0] = (int*)malloc(N_Angle * 3 * sizeof(int))) == NULL) return NULL;

	for (size_t i = 1; i < N_Angle; i++)
	{
		Alist[i] = Alist[i - 1] + 3;
	}

	float** Aconstraintlist;/*Aconstraintlist[N_Angle][4] || 0 1 present  energy const and angle const
							   2 3 present UB-energy const and UB-distance const of 1-3 atoms */

	if ((Aconstraintlist = (float**)malloc(N_Angle * sizeof(float*))) == NULL) return NULL;

	if ((Aconstraintlist[0] = (float*)malloc(N_Angle * 4 * sizeof(float))) == NULL) return NULL;

	for (size_t i = 1; i < N_Angle; i++)
	{
		Aconstraintlist[i] = Aconstraintlist[i - 1] + 4;
	}

	int** DihedralList;/*DihedralList[N_Dihedral][4]*/

	if ((DihedralList = (int**)malloc(N_Dihedral * sizeof(int*))) == NULL) return NULL;

	if ((DihedralList[0] = (int*)malloc(N_Dihedral * 4 * sizeof(int))) == NULL) return NULL;

	for (size_t i = 1; i < N_Dihedral; i++)
	{
		DihedralList[i] = DihedralList[i - 1] + 4;
	}

	float** Dihedralconst;/*Dihedralconst[N_Dihedral][3] || 0 energy const ,1 N const and 2 angle const*/

	if ((Dihedralconst = (float**)malloc(N_Dihedral * sizeof(float*))) == NULL) return NULL;

	if ((Dihedralconst[0] = (float*)malloc(N_Dihedral * 3 * sizeof(float))) == NULL) return NULL;

	for (size_t i = 1; i < N_Dihedral; i++)
	{
		Dihedralconst[i] = Dihedralconst[i - 1] + 3;
	}

	int** Improperlist;/*Improperlist[N_Improper][4] || the numbers of four atoms*/

	if ((Improperlist = (int**)malloc(N_Improper * sizeof(int*))) == NULL) return NULL;

	if ((Improperlist[0] = (int*)malloc(N_Improper * 4 * sizeof(int))) == NULL) return NULL;

	for (size_t i = 1; i < N_Improper; i++)
	{
		Improperlist[i] = Improperlist[i - 1] + 4;
	}

	float** Improperconst;/*Improperconst[N_Improper][2] || 0 energy const and 1 angle const(usually zero)*/

	if ((Improperconst = (float**)malloc(N_Improper * sizeof(float*))) == NULL) return NULL;

	if ((Improperconst[0] = (float*)malloc(N_Improper * 2 * sizeof(float))) == NULL) return NULL;

	for (size_t i = 1; i < N_Improper; i++)
	{
		Improperconst[i] = Improperconst[i - 1] + 2;
	}

	int*** shakelist;/*shakelist[N_water_Molecular][3][2]
					 Control the three distance of one water molecular and 2 are the number of two atoms */

	shakelist = malloc(N_water_Molecular * sizeof(int**) + /* level1 pointer */
		N_water_Molecular * 3 * sizeof(int*) + /* level2 pointer */
		N_water_Molecular * 3 * 2 * sizeof(int)); /* data pointer */
	for (size_t i = 0; i < N_water_Molecular; ++i)
	{
		shakelist[i] = (int**)(shakelist + N_water_Molecular) + i * 3;
		for (size_t j = 0; j < 3; ++j)
			shakelist[i][j] = (int*)(shakelist + N_water_Molecular + N_water_Molecular * 3) + i * 3 * 2 + j * 2;
	}
	/*                                             ksin   kocs                             */
	float** ksin;

	if ((ksin = (float**)malloc((N_atom + 1) * sizeof(float*))) == NULL) return NULL;

	if ((ksin[0] = (float*)malloc((N_atom + 1) * MAXK * sizeof(float))) == NULL) return NULL;

	for (int i = 1; i < N_atom + 1; i++)
	{
		ksin[i] = ksin[i - 1] + MAXK;
	}
	float** kcos;

	if ((kcos = (float**)malloc((N_atom + 1) * sizeof(float*))) == NULL) return NULL;

	if ((kcos[0] = (float*)malloc((N_atom + 1) * MAXK * sizeof(float))) == NULL) return NULL;

	for (int i = 1; i < N_atom + 1; i++)
	{
		kcos[i] = kcos[i - 1] + MAXK;
	}
	//////////////////////////build files//////////////////////////////
	FILE* fpv;
	if ((fpv = fopen("../../Output/velocity.txt", "w")) == NULL)
	{
		printf("File open error----\n");
		exit(0);
	}

	FILE* coordinate;
	if ((coordinate = fopen("../../Output/coordinate.pdb", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpforce;
	if ((fpforce = fopen("../../Output/force.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	/*"Vanderwaals", "Coulomb", "Bond", "Angle", "Dihedral","Improper","UB"*/
	FILE* fpVanderwaals;
	if ((fpVanderwaals = fopen("../../Output/VanderwaalsEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpCoulomb;
	if ((fpCoulomb = fopen("../../Output/CoulombEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpBond;
	if ((fpBond = fopen("../../Output/BondEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpAngle;
	if ((fpAngle = fopen("../../Output/AngleEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpDihedral;
	if ((fpDihedral = fopen("../../Output/DihedralEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpImproper;
	if ((fpImproper = fopen("../../Output/ImproperEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fpUB;
	if ((fpUB = fopen("../../Output/UBEnergy.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* Outpot;
	if ((Outpot = fopen("../../Output/SumPotential.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* Outkinet;
	if ((Outkinet = fopen("../../Output/kinetic.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	FILE* fwdcd;
	if ((fwdcd = fopen("../../Output/Trajectory.dcd", "wb+")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	//////////////////////////initial the system state //////////////////////////
	srand((unsigned)time(NULL));

	parameters_input(PSF, PDB, part, LJ_set, Blist, Bconstraintlist, Alist, Aconstraintlist, DihedralList,
		Dihedralconst, CoulList, shakelist, Improperlist, Improperconst, waterlist, r);//Input kinds of parameters

	Self_part_of_k_space_sum(part, &U);

	ExpiKspace(G_k, rk);

	//////////////////the numbers of nonbond atoms ///////////////////
	int temporary1 = 0, temporary2 = 0, temporary3 = 0;

	for (int i = 0; i <= N_atom; i++)
	{
		for (int it = 0; it <= N_atom; it++)
		{
			nonbond_list[i][it] = 0;

			TIP3P_Water[i][it] = 0;

		}
	}

	for (int i = 0; i <= N_atom; i++)
	{
		for (int it = 0; it <= N_atom; it++)
		{
			scale_1_4[i][it] = 0;
		}
	}

	for (int ip = 1; ip <= N_nonwater_atom; ip++)
	{
		for (int ifour = 0; ifour < N_Dihedral; ifour++)
		{
			if (DihedralList[ifour][0] == ip)
			{
				temporary1 = DihedralList[ifour][1];
				nonbond_list[ip][temporary1] = 1;
				temporary2 = DihedralList[ifour][2];
				nonbond_list[ip][temporary2] = 1;
				temporary3 = DihedralList[ifour][3];
				nonbond_list[ip][temporary3] = 1;
			}
			if (DihedralList[ifour][1] == ip)
			{
				temporary1 = DihedralList[ifour][0];
				nonbond_list[ip][temporary1] = 1;
				temporary2 = DihedralList[ifour][2];
				nonbond_list[ip][temporary2] = 1;
				temporary3 = DihedralList[ifour][3];
				nonbond_list[ip][temporary3] = 1;
			}
			if (DihedralList[ifour][2] == ip)
			{
				temporary1 = DihedralList[ifour][0];
				nonbond_list[ip][temporary1] = 1;
				temporary2 = DihedralList[ifour][1];
				nonbond_list[ip][temporary2] = 1;
				temporary3 = DihedralList[ifour][3];
				nonbond_list[ip][temporary3] = 1;
			}
			if (DihedralList[ifour][3] == ip)
			{
				temporary1 = DihedralList[ifour][0];
				nonbond_list[ip][temporary1] = 1;
				temporary2 = DihedralList[ifour][1];
				nonbond_list[ip][temporary2] = 1;
				temporary3 = DihedralList[ifour][2];
				nonbond_list[ip][temporary3] = 1;
			}
		}

	}


	/*1-4 scale*/
	int sc14_tem1 = 0, sc14_tem2 = 0;

	for (int ifour = 0; ifour < N_Dihedral; ifour++)
	{
		sc14_tem1 = DihedralList[ifour][0];
		sc14_tem2 = DihedralList[ifour][3];
		scale_1_4[sc14_tem1][sc14_tem2] = 1;
		scale_1_4[sc14_tem2][sc14_tem1] = 1;
	}

	printf("%d", CoulList[1]);
	/*check*/
	for (int i = 1; i <= N_nonwater_atom; i++)
	{
		for (int it = 1; it <= N_nonwater_atom; it++)
		{
			printf("%d", scale_1_4[i][it]);

		}
		printf("\n");
	}

	//int num_procs, max_threads;
	//num_procs = omp_get_num_procs();
	//printf("num_procs = %d\n", num_procs);

	//max_threads = omp_get_max_threads();
	//printf("max_threads = %d\n", max_threads);

	if (STEEPEST)
	{
		printf("============================================================================\n");

		steepest_descent(Control, PDB, part, &U, BF, NF, LJ_set, Blist, Bconstraintlist, Alist, Aconstraintlist, DihedralList,
			Dihedralconst, CoulList, G_k, shakelist, shakeconst, nonbond_list, Improperlist, Improperconst, TIP3P_Water, scale_1_4,
			r, v, rold, F_all, ksin, kcos, rk);

		initial_velocity(part, &U, v);/*the initialization of velocity  procedure must execute after energy minimization because in the "steep"
								   ,the constraint subroutine has been run ,the velocity has changed*/
	}
	else {
		initial_velocity(part, &U, v);

		force(Control, part, &U, BF, NF, LJ_set, Blist, Bconstraintlist, Alist, Aconstraintlist, DihedralList, Dihedralconst, CoulList,
			G_k, nonbond_list, scale_1_4, Improperlist, Improperconst, TIP3P_Water, r, F_all, ksin, kcos, rk);

	}
	/*Write the header of DCD file */
	Out_Put(0, true, r, fwdcd, coordinate, PDB);
	///////////////////////////ensemble/////////////////////////////

	for (int step = 0; step < loop; step++)
	{
		/*calculate the consuming time  of every step*/
		time_t start, finish;
		float cost = 0.;
		time(&start);

		printf("the step %6d\n", step);

		///////////////////////////////integration/////////////////////////////
		Verlet_Velocity(Control, part, &U, BF, NF, LJ_set, Blist, Bconstraintlist, Alist, Aconstraintlist, DihedralList,
			Dihedralconst, CoulList, G_k, nonbond_list, Improperlist, Improperconst, shakelist, shakeconst, TIP3P_Water, scale_1_4,
			r, rold, v, F_all, ksin, kcos, rk);

		Berendsen_Thermostat(part, &U, v);

		VWrapALL(r);


		///////////////////////////////Output energy////////////////////////////
		U.all = U.real + U.fourier + U.vander + U.bond + U.angle + U.torsion - U.self_K_SPACE_SUM + U.ureyBradley + U.improper;
		///////////////////////////////Output//////////////////////////////////////
		fprintf(Outpot, "%5d  %f\n", count_model, U.all);
		fprintf(Outkinet, "%5d  %f\n", count_model, U.kineticE);
		//count_model++;
		Out_Put(step, false, r, fwdcd, coordinate, PDB);
		///////////////////Output kinds of energy on the screnn/////////////////
		printf("the steptime  %f\n", step_t);
		printf("P_vander      %f\n", U.vander);
		printf("P_fourier     %f\n", U.fourier);
		printf("P_real        %f\n", U.real);
		printf("P_angle       %f\n", U.angle);
		printf("P_bond        %f\n", U.bond);
		printf("P_torsion     %f\n", U.torsion);
		printf("P_improper    %f\n", U.improper);
		printf("P_ureyBradley %f\n", U.ureyBradley);
		printf("P_all         %f\n", U.all);
		printf("P_kineticE    %f\n", U.kineticE);

		time(&finish);
		cost = difftime(finish, start);
		printf("step spends %.5f second\n", cost);
	}
	return 0;
}