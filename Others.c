#include "all_def.h"
#include "function.h"
float Rondom(void)
{
	//Box-Muller method 
	float x1;
	float rand1;
	float rand2;
	rand1 = (rand() % 10000) * 0.0001;//rand function
	rand2 = (rand() % 10000) * 0.0001;
	x1 = (cos(2 * PAI * rand1)) * sqrt(-2 * log(rand2));
	return x1;
}

void ExpiKspace(float* G_k, float** rk)
{
	int KSQ;
	float rkx, rky, rkz, rksq;
	float tpiLx = 2. * PAI / L, tpiLy = 2. * PAI / L, tpiLz = 2. * PAI / L;
	int kcount = 0;
	for (int nx = 0; nx <= KMAX; nx++)
	{
		rkx = nx * tpiLx;
		for (int ny = 0; ny <= KMAX; ny++)
		{
			rky = ny * tpiLy;
			for (int nz = 0; nz <= KMAX; nz++)
			{
				rkz = nz * tpiLz;
				KSQ = nx * nx + ny * ny + nz * nz;

				if ((KSQ < KSQMAX) && (KSQ != 0))
				{
					rk[kcount][0] = rkx;
					rk[kcount][1] = rky;
					rk[kcount][2] = rkz;
					rk[kcount][3] = KSQ;

					rksq = rkx * rkx + rky * rky + rkz * rkz;//K=2*¦°*n/L the square of vector
					G_k[KSQ] = 2. * PAI * exp(-0.25 / (alpha * alpha) * rksq) / rksq;
					printf("%d         %f,%f,%f,%f\n", kcount, rk[kcount][0], rk[kcount][1], rk[kcount][2], rk[kcount][3]);
					kcount++;
				}
			}
		}
	}
	printf("the kcount is %d\n", kcount);
}

void Self_part_of_k_space_sum(struct AtomType* ptc, struct Energy* E)
{
	float VS = 0.0;
	int i;

	for (i = 1; i <= N_atom; i++)
	{
		VS = VS + ptc[i].Q * ptc[i].Q;
	}

	E->self_K_SPACE_SUM = (1. / sqrt(PAI)) * alpha * VS * elec_constant;

	printf("the self energy %f\n", E->self_K_SPACE_SUM);
}

void initial_velocity(struct AtomType* ptc, struct Energy* E, float** v)
{
	FILE* fp_initial_V = NULL;
	if ((fp_initial_V = fopen("../../Review/Initial_Velocity.txt", "w")) == NULL)
	{
		printf("File open error\n");
		exit(0);
	}
	float v2 = 0;
	float sumvx = 0, sumvy = 0, sumvz = 0, sumv2 = 0;
	float tx = 0, ty = 0, tz = 0;
	float Mass_all = 0;
	E->kineticE = 0;
	for (int ip = 1; ip <= N_atom; ip++)
	{
		v[ip][0] = kT * Rondom();
		v[ip][1] = kT * Rondom();
		v[ip][2] = kT * Rondom();

		printf("velocity is %f, %f, %f\n", v[ip][0], v[ip][1], v[ip][2]);

		sumvx += v[ip][0] * ptc[ip].m;
		sumvy += v[ip][1] * ptc[ip].m;
		sumvz += v[ip][2] * ptc[ip].m;

		Mass_all += ptc[ip].m;
	}
	//keep sum of momentum is zero

	printf("sumvx   %5.3f\n", sumvx);
	printf("sumvy   %5.3f\n", sumvy);
	printf("sumvz   %5.3f\n", sumvz);

	sumvx = sumvx / Mass_all;
	sumvy = sumvy / Mass_all;
	sumvz = sumvz / Mass_all;

	printf("sumvx1   %5.3f\n", sumvx);
	printf("sumvy1   %5.3f\n", sumvy);
	printf("sumvz1   %5.3f\n", sumvz);

	for (int ip = 1; ip <= N_atom; ip++)
	{
		tx = v[ip][0] - sumvx;
		ty = v[ip][1] - sumvy;
		tz = v[ip][2] - sumvz;

		v[ip][0] = tx;
		v[ip][1] = ty;
		v[ip][2] = tz;

		printf("%f %f %f\n", v[ip][0], v[ip][1], v[ip][2]);
		v2 += ptc[ip].m * (tx * tx + ty * ty + tz * tz);
	}

	float fs;
	fs = sqrt(3. * N_atom * kT / v2);
	v2 = 0;
	for (int ip = 1; ip <= N_atom; ip++)
	{
		tx = v[ip][0] * fs;
		ty = v[ip][1] * fs;
		tz = v[ip][2] * fs;

		v[ip][0] = tx;
		v[ip][1] = ty;
		v[ip][2] = tz;

		v2 += ptc[ip].m * (tx * tx + ty * ty + tz * tz);

	}
	E->kineticE = 0.5 * v2;

	printf("The initial kinetic energy  %f\n", E->kineticE);

	//Output the initial velocity
	for (int ip = 1; ip <= N_atom; ip++)
	{
		fprintf(fp_initial_V, "%-5d,  %-8.6f,  %-8.6f,  %-8.6f\n", ip, v[ip][0], v[ip][1], v[ip][2]);
	}

	if (fclose(fp_initial_V))
	{
		printf("can not close the fp_initial_V file \n");
		exit(0);
	}
}
