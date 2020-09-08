#include "all_def.h"
#include "function.h"
void angle_bending_term(struct BAD_Control* Control, struct Energy* E, struct BondingForce* BF, int** Al, float** Ac, float** r)
{
	float dijx, dijy, dijz, dkjx, dkjy, dkjz, dikx, diky, dikz, dik2, dik, c11, c12, c22, cd, c, dtheta;
	float fix, fiy, fiz, fkx, fky, fkz, f, thetareal;
	int ip, jp, kp, ian;
	float potential_Angle = 0;
	float potential_UB = 0;
#pragma omp parallel private(dijx, dijy, dijz, dkjx, dkjy, dkjz, dikx, diky, dikz, dik2,dik, c11, c12, c22,cd,c, dtheta,fix, fiy, \
    fiz, fkx, fky, fkz,f,thetareal,ip, jp, kp, ian) shared(Al,r,Ac,BF)   reduction(+:potential_Angle,potential_UB)
	{
#pragma omp for schedule(static)
		for (ian = 0; ian < N_Angle; ian++)
		{
			ip = Al[ian][0];
			jp = Al[ian][1];
			kp = Al[ian][2];

			dijx = r[jp][0] - r[ip][0];
			dijy = r[jp][1] - r[ip][1];
			dijz = r[jp][2] - r[ip][2];

			if (dijx > 0.5 * xboxL)
			{
				dijx = dijx - xboxL;
			}
			else if (dijx < -0.5 * xboxL)
			{
				dijx = dijx + xboxL;
			}
			if (dijy > 0.5 * yboxL)
			{
				dijy = dijy - yboxL;
			}
			else if (dijy < -0.5 * yboxL)
			{
				dijy = dijy + yboxL;
			}
			if (dijz > 0.5 * zboxL)
			{
				dijz = dijz - zboxL;
			}
			else if (dijz < -0.5 * zboxL)
			{
				dijz = dijz + zboxL;
			}

			c11 = dijx * dijx + dijy * dijy + dijz * dijz;

			dkjx = r[kp][0] - r[jp][0];
			dkjy = r[kp][1] - r[jp][1];
			dkjz = r[kp][2] - r[jp][2];

			if (dkjx > 0.5 * xboxL)
			{
				dkjx = dkjx - xboxL;
			}
			else if (dkjx < -0.5 * xboxL)
			{
				dkjx = dkjx + xboxL;
			}
			if (dkjy > 0.5 * xboxL)
			{
				dkjy = dkjy - xboxL;
			}
			else if (dkjy < -0.5 * xboxL)
			{
				dkjy = dkjy + xboxL;
			}
			if (dkjz > 0.5 * xboxL)
			{
				dkjz = dkjz - xboxL;
			}
			else if (dkjz < -0.5 * xboxL)
			{
				dkjz = dkjz + xboxL;
			}

			c22 = dkjx * dkjx + dkjy * dkjy + dkjz * dkjz;

			c12 = dijx * dkjx + dijy * dkjy + dijz * dkjz;

			cd = sqrtf(c11*c22);
			
			c = c12 / cd;

			thetareal = acosf(c);

			//Control->Angle_D[ian] = 180. - thetareal * 180. / PAI;

			dtheta = (thetareal - Ac[ian][1]);

			potential_Angle += Ac[ian][0] * dtheta * dtheta;//U=k(theta-theta0)^2

			f = 2. * Ac[ian][0] * dtheta / sinf(thetareal);//-partialU/partialcos¦È=-2*k*(theta-theta0)/-sin¦È

			fix = (c12 / c11 * dijx - dkjx) * (f / cd);
			fiy = (c12 / c11 * dijy - dkjy) * (f / cd);
			fiz = (c12 / c11 * dijz - dkjz) * (f / cd);

			fkx = -(c12 / c22 * dkjx - dijx) * (f / cd);
			fky = -(c12 / c22 * dkjy - dijy) * (f / cd);
			fkz = -(c12 / c22 * dkjz - dijz) * (f / cd);
#pragma omp critical
			{
				BF[ip].F_angle[0] += fix;
				BF[ip].F_angle[1] += fiy;
				BF[ip].F_angle[2] += fiz;

				BF[jp].F_angle[0] -= (fix + fkx);
				BF[jp].F_angle[1] -= (fiy + fky);
				BF[jp].F_angle[2] -= (fiz + fkz);

				BF[kp].F_angle[0] += fkx;
				BF[kp].F_angle[1] += fky;
				BF[kp].F_angle[2] += fkz;
			}
			if (Ac[ian][2] != 0)
			{
				dikx = r[ip][0] - r[kp][0];
				diky = r[ip][1] - r[kp][1];
				dikz = r[ip][2] - r[kp][2];

				if (dikx > 0.5 * xboxL)
				{
					dikx = dikx - xboxL;
				}
				else if (dikx < -0.5 * xboxL)
				{
					dikx = dikx + xboxL;
				}
				if (diky > 0.5 * yboxL)
				{
					diky = diky - yboxL;
				}
				else if (diky < -0.5 * yboxL)
				{
					diky = diky + yboxL;
				}
				if (dikz > 0.5 * zboxL)
				{
					dikz = dikz - zboxL;
				}
				else if (dikz < -0.5 * zboxL)
				{
					dikz = dikz + zboxL;
				}

				dik2 = dikx * dikx + diky * diky + dikz * dikz;
				dik = sqrtf(dik2);

				potential_UB += Ac[ian][2] * (dik - Ac[ian][3]) * (dik - Ac[ian][3]);

				f = 2* Ac[ian][2] * (dik-Ac[ian][3]) / dik;
#pragma omp critical
				{
					BF[ip].F_ub[0] += f * dikx;
					BF[ip].F_ub[1] += f * diky;
					BF[ip].F_ub[2] += f * dikz;

					BF[kp].F_ub[0] -= f * dikx;
					BF[kp].F_ub[1] -= f * diky;
					BF[kp].F_ub[2] -= f * dikz;
				}
			}

		}

	}

	E->angle = potential_Angle;
	E->ureyBradley = potential_UB;

}