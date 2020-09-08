#include "all_def.h"
#include "function.h"
void Berendsen_Thermostat(struct AtomType* ptc, struct Energy* E, float** v)
{
	float Te = 0., Xi = 0.;
	int ip;

	for (ip = 1; ip <= N_atom; ip++)
	{
		Te += ptc[ip].m * (v[ip][0] * v[ip][0] + v[ip][1] * v[ip][1] + v[ip][2] * v[ip][2]);
	}

	E->kineticE = Te;
	Te = Te / (3. * N_atom);
	Xi = sqrt(1 + step_t * (kT / Te - 1) / taoT);
	E->kineticE *= (0.5 * Xi * Xi);

	for (ip = 1; ip <= N_atom; ip++)
	{
		v[ip][0] *= Xi;
		v[ip][1] *= Xi;
		v[ip][2] *= Xi;
	}

}