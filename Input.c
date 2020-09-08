#include "all_def.h"
#include "function.h"
int parameters_input(struct PSFINPUT* PSF, struct PDB_format* PDB, struct AtomType* ptc, float(*LJ)[4], int** Blist, float** Bconst,
	int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst, int* CoulList, int*** shakelist, int** Improperlist,
	float** Improperconst, int** waterlist, float** r)
{
	/////////////////////////导入Blist//////////////////////////////////////

	FILE* Bl = fopen("../../Input/Blist.txt", "r");
	if (Bl == NULL)
	{
		printf("Blist open error\n");
		exit(0);
	}
	for (int ip = 0; ip < N_Bond; ip++)
	{
		if (!(fscanf(Bl, "%d %d", &Blist[ip][0], &Blist[ip][1])))
			printf("scanf bl wrong\n");
	}

	if (fclose(Bl))
	{
		printf("can not close Blist \n");
		exit(0);
	}

	FILE* Bl_review = fopen("../../Review/review_Blist.txt", "w");
	if (Bl_review == NULL)
	{
		printf("Bl_review open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Bond; ip++)
	{
		fprintf(Bl_review, "%d   %d\n", Blist[ip][0], Blist[ip][1]);
	}


	if (fclose(Bl_review))
	{
		printf("can not close Bl_review \n");
		exit(0);
	}
	//Real-time output 

	for (int ip = 0; ip < N_Bond; ip++)
	{
		printf("B list  %d  %d  \n", Blist[ip][0], Blist[ip][1]);

	}

	//////////////////////////////Bconstraintlist///////////////////////////////////
	float kb = 0, r_bond = 0;
	FILE* Bc = fopen("../../Input/Bondconst.txt", "r");
	if (Bc == NULL)
	{
		printf("Bondconst open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Bond; ip++)
	{
		fscanf(Bc, "%lf  %lf", &kb, &r_bond);//0 presents r 1 presents k
		Bconst[ip][0] = r_bond / SIGMA;
		Bconst[ip][1] = (kb / EPSLON) * SIGMA * SIGMA;
	}

	if (fclose(Bc))
	{
		printf("can not close Bondconst \n");
		exit(0);
	}

	FILE* Bc_review = fopen("../../Review/review_Bconstraintlist.txt", "w");
	if (Bc_review == NULL)
	{
		printf("can not open Bondconst_review \n");
		exit(0);
	}

	for (int ip = 0; ip < N_Bond; ip++)
	{
		fprintf(Bc_review, "%f, %f\n", Bconst[ip][0], Bconst[ip][1]);
	}

	if (fclose(Bc_review))
	{
		printf("can not close the Bondconst_review \n");
		exit(0);
	}

	//Real-time output
	for (int ip = 0; ip < N_Bond; ip++)
	{
		printf("Bconst  %f, %f\n", Bconst[ip][0], Bconst[ip][1]);
	}

	////////////////////////Anglelist//////////////////////////////

	FILE* Al = fopen("../../Input/Alist.txt", "r");
	if (Al == NULL)
	{
		printf("Alist open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Angle; ip++)
	{
		fscanf(Al, "%d %d  %d", &Alist[ip][0], &Alist[ip][1], &Alist[ip][2]);

	}

	if (fclose(Al))
	{
		printf("can not close Al \n");
		exit(0);
	}

	FILE* Al_review = fopen("../../Review/review_Alist.txt", "w");
	if (Al_review == NULL)
	{
		printf("Al_re open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Angle; ip++)
	{
		fprintf(Al_review, "%d   %d    %d\n", Alist[ip][0], Alist[ip][1], Alist[ip][2]);
	}

	if (fclose(Al_review))
	{
		printf("can not close Al_re \n");
		exit(0);
	}
	//Real-time output
	for (int ip = 0; ip < N_Angle; ip++)
	{
		printf("Alist%d  %d  %d\n", Alist[ip][0], Alist[ip][1], Alist[ip][2]);

	}

	////////////////////////Angleconstraint//////////////////////////////////
	float theta0 = 0, ktheta = 0;
	float theta01 = 0, ktheta1 = 0;
	FILE* Ac = fopen("../../Input/Angleconst.txt", "r");
	if (Ac == NULL)
	{
		printf("Ac open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Angle; ip++)
	{
		fscanf(Ac, "%lf  %lf %lf %lf", &ktheta, &theta0, &ktheta1, &theta01);// 0 k  1 angle
		Aconst[ip][0] = ktheta / EPSLON;
		Aconst[ip][1] = (180. - theta0) * PAI / 180.;
		Aconst[ip][2] = ktheta1 / EPSLON * SIGMA * SIGMA;
		Aconst[ip][3] = theta01 / SIGMA;
	}

	if (fclose(Ac))
	{
		printf("can not close Ac \n");
		exit(0);
	}

	FILE* Ac_review = fopen("../../Review/review_Angleconstraintlist.txt", "w");
	if (Ac_review == NULL)
	{
		printf("can not open Ac_re \n");
		exit(0);
	}

	for (int ip = 0; ip < N_Angle; ip++)
	{
		fprintf(Ac_review, "%6.4f  %6.4f  %6.4f  %6.4f\n", Aconst[ip][0], Aconst[ip][1], Aconst[ip][2], Aconst[ip][3]);
	}

	if (fclose(Ac_review))
	{
		printf("can not close Ac_re \n");
		exit(0);
	}

	for (int ip = 0; ip < N_Angle; ip++)
	{
		printf("Aconst %6.4f  %6.4f  %6.4f  %6.4f\n", Aconst[ip][0], Aconst[ip][1], Aconst[ip][2], Aconst[ip][3]);
	}

	////////////////////////////Dihedrallist///////////////////////////

	FILE* Dihedrallist = fopen("../../Input/Dihedrallist.txt", "r");

	if (Dihedrallist == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Dihedral; ip++)
	{
		fscanf(Dihedrallist, "%d %d %d %d", &DihedralList[ip][0], &DihedralList[ip][1], &DihedralList[ip][2], &DihedralList[ip][3]);
	}

	if (fclose(Dihedrallist))
	{
		printf("can not close the file \n");
		exit(0);
	}

	FILE* Dihedrallist_review = fopen("../../Review/review_Dihedrallist.txt", "w");
	if (Dihedrallist_review == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Dihedral; ip++)
	{
		fprintf(Dihedrallist_review, "%d   %d    %d   %d\n", DihedralList[ip][0], DihedralList[ip][1], DihedralList[ip][2],
			DihedralList[ip][3]);
	}

	if (fclose(Dihedrallist_review))
	{
		printf("can not close the file \n");
		exit(0);
	}
	//Real-time output
	for (int ip = 0; ip < N_Dihedral; ip++)
	{
		printf("Dihedrallist %d  %d  %d  %d\n", DihedralList[ip][0], DihedralList[ip][1], DihedralList[ip][2], DihedralList[ip][3]);
	}

	/////////////////////Dihedralconst/////////////////////////////

	float kd = 0, thetaD;
	FILE* Dconst = fopen("../../Input/Dihedralconst.txt", "r");
	if (Dconst == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Dihedral; ip++)
	{//0表示Vn/2  1 为 n  2为 γ (角度)
		fscanf(Dconst, "%lf %lf %lf ", &kd, &Dihedralconst[ip][1], &thetaD);
		Dihedralconst[ip][0] = kd / EPSLON;
		Dihedralconst[ip][2] = thetaD * PAI / 180.;
	}

	if (fclose(Dconst))
	{
		printf("can not close the file \n");
		exit(0);
	}

	FILE* Dihedralconst_review = fopen("../../Review/review_Dihedralconst.txt", "w");
	if (Dihedralconst_review == NULL)
	{
		printf("File open error\n");
		exit(0);
	}

	for (int ip = 0; ip < N_Dihedral; ip++)
	{
		fprintf(Dihedralconst_review, "%f   %f    %f  \n", Dihedralconst[ip][0], Dihedralconst[ip][1], Dihedralconst[ip][2]);
	}

	if (fclose(Dihedralconst_review))
	{
		printf("can not close the file \n");
		exit(0);
	}
	//实时输出

	for (int ip = 0; ip < N_Dihedral; ip++)
	{
		printf("Dihedralconst   %f  %f  %f  \n", Dihedralconst[ip][0], Dihedralconst[ip][1], Dihedralconst[ip][2]);
	}

	////////////////////////////////////PDB文件格式输入////////////////////////////////////////////
	float linshi1 = 0, linshi2 = 0;
	FILE* fppdb;
	if ((fppdb = fopen("../../Input/PDB.txt", "r")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	for (int ip = 1; ip <= N_atom; ip++)
	{
		fscanf(fppdb, "%s %d %s %s %s %d %lf %lf %lf %lf %lf %s %s", PDB[ip].atom, &PDB[ip].serial, PDB[ip].name, PDB[ip].resName,
			PDB[ip].chainID, &PDB[ip].resSeq, &r[ip][0], &r[ip][1], &r[ip][2], &linshi1, &linshi2, PDB[ip].segID, PDB[ip].element);
	}
	if (fclose(fppdb))
	{
		printf("the fppdb file can not close\n");
		exit(0);
	}
	//transition 
	for (int ip = 1; ip <= N_atom; ip++)
	{
		r[ip][0] = r[ip][0] / (SIGMA);
		r[ip][1] = r[ip][1] / (SIGMA);
		r[ip][2] = r[ip][2] / (SIGMA);
	}

	FILE* fp_review_initial_location;
	if ((fp_review_initial_location = fopen("../../Review/review_initial_location.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	for (int ip = 1; ip <= N_atom; ip++)
	{
		fprintf(fp_review_initial_location, " %d,  %f,  %f,  %f\n", ip, r[ip][0], r[ip][1], r[ip][2]);
	}
	if (fclose(fp_review_initial_location))
	{
		printf("can not close the file \n");
		exit(0);
	}
	//坐标检查
	for (int ip = 1; ip <= N_atom; ip++)
	{
		printf("%s %d %s %s %s %d %f %f %f %s\n", PDB[ip].atom, PDB[ip].serial, PDB[ip].name, PDB[ip].resName, PDB[ip].chainID, PDB[ip].resSeq,
			r[ip][0], r[ip][1], r[ip][2], PDB[ip].element);
	}

	//////////////////////////////PSF文件输入质量和电荷///////////
	FILE* psf_input;
	if ((psf_input = fopen("../../Input/PSF.txt", "r")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	for (int ip = 1; ip <= N_atom; ip++)
	{
		fscanf(psf_input, "%d %s %d %s %s %s %lf %lf %d", &PSF[ip].atomID, PSF[ip].segmentname, &PSF[ip].residueID, PSF[ip].residuename,
			PSF[ip].atomname, PSF[ip].atomtype, &ptc[ip].Q, &ptc[ip].m, &PSF[ip].unused);
	}

	if (fclose(psf_input))
	{
		printf("the psf_input file can not close\n");
		exit(0);
	}

	for (int ip = 1; ip <= N_atom; ip++)
	{
		printf("lxm%d  %s  %d  %s  %s  %s  %f  %f  %d\n", PSF[ip].atomID, PSF[ip].segmentname, PSF[ip].residueID, PSF[ip].residuename,
			PSF[ip].atomname, PSF[ip].atomtype, ptc[ip].Q, ptc[ip].m, PSF[ip].unused);
	}
	//////////////////////////////LJlist///////////////////////  
	char* LJlist[] = { "HCA1","HCA2","HCA3","HCP1","HT","CC30A","CC31A","CC32A","CC33A","CC326A","HCA25A",
						"CC325A","CC325B" };
	for (int ip = 1; ip <= N_atom; ip++)
	{
		for (int i = 0; i < 56; i++)
		{
			if (strcmp(PSF[ip].atomtype, LJlist[i]) == 0)
			{
				ptc[ip].type = i;

				break;
			}
			if (i == 55)
			{
				printf("did not find the type\n");
				printf("%d\n", ip);
				exit(0);
			}


		}

	}

	for (int ip = 1; ip <= N_atom; ip++)
	{
		printf("%d  ", ptc[ip].type);
	}
	printf("\n");
	FILE* fpljlist = fopen("../../Input/LJ interaction set.txt", "r");
	float e = 0, r0 = 0, e1 = 0, r1 = 0;
	if (fpljlist == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	for (int i = 0; i < LJ_list; i++)
	{
		fscanf(fpljlist, "%lf  %lf %lf  %lf", &e, &r0, &e1, &r1);
		LJ[i][0] = r0 / SIGMA;
		LJ[i][1] = e / EPSLON;
		LJ[i][2] = r1 / SIGMA;
		LJ[i][3] = e1 / EPSLON;
	}
	if (fclose(fpljlist))
	{
		printf("the fpljlist file can not close\n");
		exit(0);
	}

	//将电荷原子编号导入CoulList
	int count = 0;
	for (int ip = 1; ip <= N_atom; ip++)
	{
		if (fabs(ptc[ip].Q) > TOL)
		{
			CoulList[count] = ip;
			printf("%d  ", count);
			count++;
		}
	}
	printf("\n");
	printf("==================================%d\n", count);

	for (int i = 0; i < count; i++)
	{
		printf("%d\n", CoulList[i]);
	}
	/////////////////////////review_mass_and_charge///////////////////////////////
	FILE* fpreview = fopen("../../Review/review_mass_and_charge.txt", "w");
	if (fpreview == NULL)
	{
		printf("File open error\n");
		exit(0);
	};
	for (int ip = 1; ip <= N_atom; ip++)
	{
		fprintf(fpreview, "%f,  %f\n", ptc[ip].m, ptc[ip].Q);
	}

	if (fclose(fpreview))
	{
		printf("can not close the file \n");
		exit(0);
	}
	//Real-time output
	for (int ip = 1; ip <= N_atom; ip++)
	{
		printf("mass charge,  %d,  %f,   %f\n", ip, ptc[ip].m, ptc[ip].Q);
	}

	return 0;
}