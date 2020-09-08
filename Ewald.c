#include "all_def.h"
#include "function.h"
void Coulombforce_term(struct AtomType* ptc, struct Energy* E, struct NonbondingForce* NF, int* CoulList, float* G_k,
	int** nonbond_list, int** scale_1_4, float** r, float** ksin, float** kcos, float** rk)
{
	/*                                  realspace parameters                                  */
	float dx, dy, dz, r_ab;
	float charge_2;
	float f_recip, f_real;
	float potential_reciprocal = 0, potential_real = 0;
	float temp2, tem;
	int a, b, i, h, temp3;
	/*                                 Reciprocal    parameters                               */

	size_t ip = 0, kcount = 0, count = 0;
	float V = xboxL * yboxL * zboxL;
	float factor;

	float* SumQSin;

	if ((SumQSin = (float*)calloc(MAXK, sizeof(float))) == NULL) exit(1);

	float* SumQCos;

	if ((SumQCos = (float*)calloc(MAXK, sizeof(float))) == NULL) exit(1);

	int nx, ny, nz;

	/*                                      Real Space                                 */
#pragma omp parallel private(dx, dy, dz, r_ab,charge_2,f_real,i,h,a,b,temp2, tem) shared(V,CoulList,nonbond_list,scale_1_4,\
	ptc,NF,ksin,kcos,SumQSin,SumQCos,rk) reduction(+:potential_real) 
	{
#pragma omp for schedule(static)
		for (i = 0; i < Coul_List; i++)
		{
			for (h = i + 1; h < Coul_List; h++)
			{
				a = CoulList[i];//atom number
				b = CoulList[h];//atom number

				if ((nonbond_list[a][b] == 0) || (scale_1_4[a][b] == 1))
				{
					dx = r[a][0] - r[b][0];
					dy = r[a][1] - r[b][1];
					dz = r[a][2] - r[b][2];

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

					r_ab = sqrtf(dx * dx + dy * dy + dz * dz);

					if (r_ab < Elec_Cut)
					{
						charge_2 = ptc[a].Q * ptc[b].Q;

						temp2 = (r_ab * 2. * alpha * expf(-alpha * alpha * (double)r_ab * r_ab) / sqrtf(PAI) + erfcf(alpha * r_ab)) / ((double)r_ab * r_ab);

						f_real = elec_constant * temp2;

#pragma omp critical 
						{
							NF[a].F_elec[0] = NF[a].F_elec[0] + charge_2 * f_real * dx / r_ab;
							NF[a].F_elec[1] = NF[a].F_elec[1] + charge_2 * f_real * dy / r_ab;
							NF[a].F_elec[2] = NF[a].F_elec[2] + charge_2 * f_real * dz / r_ab;

							NF[b].F_elec[0] = NF[b].F_elec[0] + charge_2 * f_real * dx / -r_ab;
							NF[b].F_elec[1] = NF[b].F_elec[1] + charge_2 * f_real * dy / -r_ab;
							NF[b].F_elec[2] = NF[b].F_elec[2] + charge_2 * f_real * dz / -r_ab;
						}
						//Potential		
						potential_real = potential_real + elec_constant * charge_2 * erfcf(alpha * r_ab) / r_ab;
					}
				}
			}

		}

#pragma omp for private(ip,tem) schedule(static) 
		for (count = 0; count < KSQ_num; count++)
		{
			for (a = 0; a < Coul_List; a++)
			{
				ip = CoulList[a];
				tem = rk[count][0] * r[ip][0] + rk[count][1] * r[ip][1] + rk[count][2] * r[ip][2];
				ksin[ip][count] = sin(tem);
				kcos[ip][count] = cos(tem);
				SumQSin[count] += ptc[ip].Q * ksin[ip][count];
				SumQCos[count] += ptc[ip].Q * kcos[ip][count];
			}
		}

#pragma omp for private(ip,kcount,count,factor,f_recip,temp3) schedule(static) reduction(+:potential_reciprocal)
		for (a = 0; a < Coul_List; a++)
		{
			ip = CoulList[a];
			kcount = 0;

			for (size_t count = 0; count < KSQ_num; count++)
			{
				if (rk[count][0] < 0.1)//this formula is equivalent to rk[count][0] ==0
				{
					factor = 1.0;
				}
				else
				{
					factor = 2.0;
				}
				temp3 = (int)rk[count][3];

				f_recip = 2.0 * factor * sqrtf(elec_constant) * ptc[ip].Q * G_k[temp3] * ((double)SumQCos[kcount] * ksin[ip][kcount]
					+ (double)SumQSin[kcount] * kcos[ip][kcount]) / V;

				NF[ip].F_elec[0] -= f_recip * rk[count][0];
				NF[ip].F_elec[1] -= f_recip * rk[count][1];
				NF[ip].F_elec[2] -= f_recip * rk[count][2];

				if (ip == 1)
				{//Conjugate complex,exp(-i*k.ri) and exp(i*k.ri)
					potential_reciprocal += sqrt(elec_constant) * G_k[temp3] * ((double)SumQCos[kcount] * SumQCos[kcount] + (double)SumQSin[kcount] * SumQSin[kcount]) / V;

				}
				kcount++;
			}
		}
	}
	E->fourier = potential_reciprocal;
	E->real = potential_real;
	free(SumQSin);
	free(SumQCos);
}