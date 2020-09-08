#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>
#include<string.h>

/*------------------------------use charmm22 force field----------------------------*/
/*---------------------------------temperature 310.65K------------------------------*/

#define loop 100000  /*loop times*/
#define SIGMA 1.7682 /*Oxygen's length const from water - unit A*/
#define EPSLON 0.1521/*Oxygen's energy const from water - unit kcal/mol*/
#define kT 4.0586    /*multiply Boltzmann constant by temperature -1.38064852 *pow(10,-23) J/k*/

#define N_atom 26 /*all numbers of atoms*/
#define N_nonwater_atom 26/*all numbers of atoms which do not include water*/
#define first_water_number 76  /*number of the first water*/
#define N_water_Molecular 148  /*numbers of all water moleculars*/
#define Coul_List 26     /*the numbers of elec atoms*/


#define step_t 0.00112679      /*step time*/
#define xboxL 35.00 /*should be reduced unit*/
#define yboxL 35.00 /*should be reduced unit*/
#define zboxL 35.00 /*should be reduced unit*/
#define L 10.
#define alpha 6./L


#define R_CUT 12./SIGMA                /*reduced unit*/
#define Elec_Cut 12./SIGMA             /*reduced unit*/


//bond numbers
#define N_Bond 25
#define N_Dihedral 63
#define N_Angle 48
#define N_Improper 3
#define LJ_list 13            /*categories of atoms*/

// Constraint temperature
#define taoT 0.5

//cell list 

#define Ncell 1331 /*the numbers of all cell*/
#define ncellx 11
#define ncelly 11
#define ncellz 11
#define SIZE_LIST 750
//water const
#define NATOM_H2O 3
#define HO_bond 0.9572/SIGMA
#define HH_bond 1.5139/SIGMA

/*#define EPSLON_0 8.854187817*pow(10, -12)       #define NA 6.02214076*pow(10, 23)*/
#define elec_constant 1230.082436/*(1.6*pow(10, -19)*1.6*pow(10, -19)) / ((EPSLON * 4184. / NA) * SIGMA*pow(10, -10) * 4*¦°*EPSLON_0)*/ 
#define PAI 3.141592653/*¦Ð*/
#define KMAX 4// reciprocal space  
#define KSQMAX 64
#define MAXK 405
#define KSQ_num 124
//Constraint parameter
#define Iter_MAX 200
#define TOL 1E-9

//Steepest const
#define STEEPEST 0
#define STEEPEST_STEP 10000
#define CONSTRAIN false

#define Gup_time 10

