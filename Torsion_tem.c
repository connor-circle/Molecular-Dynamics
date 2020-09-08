#include "all_def.h"
#include "function.h"
void torsion_angle_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, \
	 int** DihedralList,float** Dihedralconst, float** r)
{
	int i1, i2, i3, i4, dih;
	float dr1x, dr1y, dr1z, dr2x, dr2y, dr2z, dr3x, dr3y, dr3z;
	float c11, c12, c13, c22, c23, c33;
	float ca, cb1, cb2, cd, c;
	float t1, t2, t3, t4, t5, t6, cr1, cr2;
	float w1x, w1y, w1z, w2x, w2y, w2z;
	float f, fia, f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z;
	float potential_torsion = 0;
#pragma omp parallel private(i1, i2, i3, i4,dih,dr1x, dr1y, dr1z, dr2x, dr2y, dr2z, dr3x, dr3y, dr3z,c11, c12, c13, c22, c23, c33,\
   ca, cb1, cb2, cd, c,t1, t2, t3, t4, t5, t6, cr1, cr2,w1x, w1y, w1z, w2x, w2y, w2z,f,fia,f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z)\
    shared(DihedralList,r,BF) reduction(+:potential_torsion)
	{
#pragma omp for schedule(static)
		for (dih = 0; dih < N_Dihedral; dih++)
		{
			i1 = DihedralList[dih][0];
			i2 = DihedralList[dih][1];
			i3 = DihedralList[dih][2];
			i4 = DihedralList[dih][3];

			// 1st bond
			dr1x = r[i2][0] - r[i1][0];
			dr1y = r[i2][1] - r[i1][1];
			dr1z = r[i2][2] - r[i1][2];

			if (dr1x > 0.5 * xboxL)
			{
				dr1x = dr1x - xboxL;
			}
			else if (dr1x < -0.5 * xboxL)
			{
				dr1x = dr1x + xboxL;
			}
			if (dr1y > 0.5 * yboxL)
			{
				dr1y = dr1y - yboxL;
			}
			else if (dr1y < -0.5 * yboxL)
			{
				dr1y = dr1y + yboxL;
			}
			if (dr1z > 0.5 * zboxL)
			{
				dr1z = dr1z - zboxL;
			}
			else if (dr1z < -0.5 * zboxL)
			{
				dr1z = dr1z + zboxL;
			}

			// 2nd bond
			dr2x = r[i3][0] - r[i2][0];
			dr2y = r[i3][1] - r[i2][1];
			dr2z = r[i3][2] - r[i2][2];

			if (dr2x > 0.5 * xboxL)
			{
				dr2x = dr1x - xboxL;
			}
			else if (dr2x < -0.5 * xboxL)
			{
				dr2x = dr2x + xboxL;
			}
			if (dr2y > 0.5 * yboxL)
			{
				dr2y = dr2y - yboxL;
			}
			else if (dr1y < -0.5 * yboxL)
			{
				dr2y = dr2y + yboxL;
			}
			if (dr2z > 0.5 * zboxL)
			{
				dr2z = dr2z - zboxL;
			}
			else if (dr2z < -0.5 * zboxL)
			{
				dr2z = dr2z + zboxL;
			}

			// 3rd bond
			dr3x = r[i4][0] - r[i3][0];
			dr3y = r[i4][1] - r[i3][1];
			dr3z = r[i4][2] - r[i3][2];

			if (dr3x > 0.5 * xboxL)
			{
				dr3x = dr3x - xboxL;
			}
			else if (dr3x < -0.5 * xboxL)
			{
				dr3x = dr3x + xboxL;
			}
			if (dr3y > 0.5 * yboxL)
			{
				dr3y = dr3y - yboxL;
			}
			else if (dr3y < -0.5 * yboxL)
			{
				dr3y = dr3y + yboxL;
			}
			if (dr3z > 0.5 * zboxL)
			{
				dr3z = dr3z - zboxL;
			}
			else if (dr3z < -0.5 * zboxL)
			{
				dr3z = dr3z + zboxL;
			}

			c11 = dr1x * dr1x + dr1y * dr1y + dr1z * dr1z;
			c12 = dr1x * dr2x + dr1y * dr2y + dr1z * dr2z;
			c13 = dr1x * dr3x + dr1y * dr3y + dr1z * dr3z;
			c22 = dr2x * dr2x + dr2y * dr2y + dr2z * dr2z;
			c23 = dr2x * dr3x + dr2y * dr3y + dr2z * dr3z;
			c33 = dr3x * dr3x + dr3y * dr3y + dr3z * dr3z;

			ca = c13 * c22 - c12 * c23;
			cb1 = c11 * c22 - c12 * c12;
			cb2 = c22 * c33 - c23 * c23;

			cd = sqrt((double)cb1 * cb2);
			c = ca / cd;
			if (c < -1)
			{
				c = -1;
			}
			if (c > 1)
			{
				c = 1;
			}
			fia = acos(c);

			//Control->Dihedral_D[dih] = fia * 180 / PAI;

			f = -(double)Dihedralconst[dih][0] * Dihedralconst[dih][1] * sin((double)Dihedralconst[dih][1] * fia - Dihedralconst[dih][2]) / sin(fia);
			t1 = ca;
			t2 = c11 * c23 - c12 * c13;
			t3 = -cb1;
			t4 = cb2;
			t5 = c13 * c23 - c12 * c33;
			t6 = -ca;
			cr1 = c12 / c22;
			cr2 = c23 / c22;

			w1x = t1 * dr1x + t2 * dr2x + t3 * dr3x;
			w1y = t1 * dr1y + t2 * dr2y + t3 * dr3y;
			w1z = t1 * dr1z + t2 * dr2z + t3 * dr3z;
			w1x = w1x * f * c22 / (cd * cb1);
			w1y = w1y * f * c22 / (cd * cb1);
			w1z = w1z * f * c22 / (cd * cb1);

			w2x = t4 * dr1x + t5 * dr2x + t6 * dr3x;
			w2y = t4 * dr1y + t5 * dr2y + t6 * dr3y;
			w2z = t4 * dr1z + t5 * dr2z + t6 * dr3z;
			w2x = w2x * f * c22 / (cd * cb2);
			w2y = w2y * f * c22 / (cd * cb2);
			w2z = w2z * f * c22 / (cd * cb2);
			//计算四个原子上的力

			//1st
			f1x = w1x;
			f1y = w1y;
			f1z = w1z;
			//2nd
			f2x = cr2 * w2x - (1. + cr1) * w1x;
			f2y = cr2 * w2y - (1. + cr1) * w1y;
			f2z = cr2 * w2z - (1. + cr1) * w1z;
			//3rd
			f3x = cr1 * w1x - (1. + cr2) * w2x;
			f3y = cr1 * w1y - (1. + cr2) * w2y;
			f3z = cr1 * w1z - (1. + cr2) * w2z;
			//4th
			f4x = w2x;
			f4y = w2y;
			f4z = w2z;
#pragma omp critical
			{

				BF[i1].F_tor[0] += f1x;
				BF[i1].F_tor[1] += f1y;
				BF[i1].F_tor[2] += f1z;

				BF[i2].F_tor[0] += f2x;
				BF[i2].F_tor[1] += f2y;
				BF[i2].F_tor[2] += f2z;

				BF[i3].F_tor[0] += f3x;
				BF[i3].F_tor[1] += f3y;
				BF[i3].F_tor[2] += f3z;

				BF[i4].F_tor[0] += f4x;
				BF[i4].F_tor[1] += f4y;
				BF[i4].F_tor[2] += f4z;
			}
			potential_torsion += Dihedralconst[dih][0] * (1 + cos((double)Dihedralconst[dih][1] * fia - Dihedralconst[dih][2]));

		}
	}
	E->torsion = potential_torsion;

}