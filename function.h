#include "all_def.h"
struct AtomType
{
	int  type;//LJ type
	float m;
	float Q;
};//part[N_atom +1]

struct NonbondingForce
{
	float F_vander[3];
	float F_elec[3];
};//NF[N_atom + 1]

struct BondingForce
{
	float F_tor[3];
	float F_impro[3];
	float F_ub[3];
	float F_angle[3];
	float F_bond[3];
};//BF[N_nonwater_atom + 1]

struct Energy
{
	float vander;
	float real;
	float fourier;
	float all;
	float bond;
	float angle;
	float torsion;
	float improper;
	float ureyBradley;
	float self_K_SPACE_SUM;
	float kineticE;
}U;

struct BAD_Control/*------------noting:N_BondN_AngleN_Dihedral--------------*/
{
	float Bond_L[N_Bond];//real bond length
	float Angle_D[N_Angle];//real degree
	float Dihedral_D[N_Dihedral];//real dihedral degree

};//Control[1]

struct PDB_format
{
	int resSeq;//residue number in a protein or peptides

	int serial;//atom number 

	char atom[5];//ATOM

	char name[5];//atom name

	char resName[4];//residue name

	char chainID[2];//chain number eg.A B New

	char element[3];//element name

	char segID[5];//VMD used
};//PDB[N_atom +1]

struct PSFINPUT
{
	int atomID;
	int unused;
	int residueID;

	char segmentname[3];
	char residuename[4];

	char atomname[5];
	char atomtype[5];
};//PSF[N_atom +1]

float Rondom(void);

void ExpiKspace(float* G_k, float** rK);

void Self_part_of_k_space_sum(struct AtomType* ptc, struct Energy* E);

void initial_velocity(struct AtomType* ptc, struct Energy* E, float** v);

void bond_stretching_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, int** Bl, float** Bc, float** r);

void angle_bending_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, int** Al, float** Ac, float** r);

void torsion_angle_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, int** DihedralList, float** Dihedralconst, float** r);

void Improper_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, int** ImproperList, float** Improperconst, float** r);

void Vanderwaalsforce_term(struct AtomType* ptc, struct Energy* E, struct NonbondingForce* NF, float(*LJ)[4], int** TIP3P_Water, int** nonbond_list, int** scale_1_4, float** r);

void Coulombforce_term(struct AtomType* ptc, struct Energy* E, struct NonbondingForce* NF, int* CoulList, float* G_k, int** nonbond_list, int** scale_1_4, float** r, float** ksin, float** kcos, float** rk);

void Verlet_Velocity(struct BAD_Control* Control, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF, struct NonbondingForce* NF, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst, int* CoulList, float* G_k, int** nonbond_list, int** Improperlist, float** Improperconst, int*** shakelist, float* shakeconst, int** TIP3P_Water, int** scale_1_4, float** r, float** rold, float** v, float** F_all, float** ksin, float** kcos, float** rk);

void steepest_descent(struct BAD_Control* Control, struct PDB_format* PDB, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF, struct NonbondingForce* NF, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst, int* CoulList, float* G_k, int*** shakelist, float* shakeconst, int** nonbond_list, int** Improperlist, float** Improperconst, int** TIP3P_Water, int** scale_1_4, float** r, float** v, float** rold, float** F_all, float** ksin, float** kcos, float** rk);

void steepest_descentB(struct BAD_Control* Control, struct PDB_format* PDB, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF, struct NonbondingForce* NF, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst, int* CoulList, float* G_k, int*** shakelist, float* shakeconst, int** nonbond_list, int** Improperlist, float** Improperconst, int** TIP3P_Water, int** scale_1_4, float** r, float** v, float** rold, float** F_all, float** ksin, float** kcos, float** rk);

int parameters_input(struct PSFINPUT* PSF, struct PDB_format* PDB, struct AtomType* ptc, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst, int* CoulList, int*** shakelist, int** Improperlist, float** Improperconst, int** waterlist, float** r);

void force(struct BAD_Control* Control, struct AtomType* ptc, struct Energy* E, struct BondingForce* BF, struct NonbondingForce* NF, float(*LJ)[4], int** Blist, float** Bconst, int** Alist, float** Aconst, int** DihedralList, float** Dihedralconst, int* CoulList, float* G_k, int** nonbond_list, int** scale_1_4, int** Improperlist, float** Improperconst, int** TIP3P_Water, float** r, float** F_all, float** ksin, float** kcos, float** rk);

void VWrapALL(float** r);

void Berendsen_Thermostat(struct AtomType* ptc, struct Energy* E, float** v);

void Out_Put(int step, bool option, float** r, FILE * fwdcd, FILE * coordinate, struct PDB_format* PDB);
