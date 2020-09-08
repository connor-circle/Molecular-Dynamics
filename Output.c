#include "all_def.h"
#include "function.h"
void Out_Put(int step, bool option, float** r, FILE* fwdcd, FILE* coordinate, struct PDB_format* PDB)
{
	int count_model = 1;
	if (step)
	{
		fprintf(coordinate, "%s%8s%d\n", "MODEL", " ", count_model);
		for (int ip = 1; ip <= N_atom; ip++)
		{
			fprintf(coordinate, "%-6s%5d%1s%-4s%1s%3s%1s%1s%4d%4s%8.3f%8.3f%8.3f%18s%-4s%2s\n", "ATOM", PDB[ip].serial, " ",
				PDB[ip].name, " ", PDB[ip].resName, " ", PDB[ip].chainID, PDB[ip].resSeq, " ", r[ip][0] * SIGMA, r[ip][1] * SIGMA,
				r[ip][2] * SIGMA, " ", PDB[ip].segID, PDB[ip].element);
		}
		fprintf(coordinate, "%s%6s%d%6s%s %s   %d\n", "TER", " ", N_atom + 1, " ", "MET", PDB[N_atom].chainID, PDB[N_atom].resSeq);
		fprintf(coordinate, "%s\n", "ENDMDL");
		count_model++;

	}
	if (option)
	{
		int out_integer;
		float out_float;
		char title_string[200];
		out_integer = 84;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);
		strcpy(title_string, "CORD");
		fwrite((char*)title_string, sizeof(char), 4, fwdcd);
		out_integer = 0;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = Gup_time;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = 0;
		for (int i = 0; i < 5; i++)
		{
			fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);
		}

		out_float = (float)step_t;
		fwrite((char*)&out_float, sizeof(float), 1, fwdcd);

		out_integer = 1;

		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = 0;
		for (int i = 0; i < 8; i++)
		{
			fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);
		}
		out_integer = 22;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd); /* CHARMM version 22 */

		out_integer = 84;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = 164;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = 2;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		char title1;
		title1 = "Coordination of Atoms";
		char title2;
		title2 = "Created by YuanRuiKang in 20/5/2019";

		fwrite((char*)title1, 80, 1, fwdcd);
		fwrite((char*)title2, 80, 1, fwdcd);

		out_integer = 164;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = 4;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = N_atom;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

		out_integer = 4;
		fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);
	}
	/*
	* Write a timestep to a dcd file
	* Input: fd - a file struct for which a dcd header has already been written
	*       curframe: Count of frames written to this file, starting with 1.
	*        curstep: Count of timesteps elapsed = istart + curframe * nsavc.
	*         natoms: number of elements in x, y, z arrays
	*        x, y, z: pointers to atom coordinates
	* Output: 0 on success, negative error code on failure.
	* Side effects: coordinates are written to the dcd file.
	*/

	int out_integer = 48;
	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);
	/*unitcell - space for six floats to hold the unit cell data.
	  Not set if no unit cell data is present.*/

	float dtemp;
	dtemp = xboxL;
	fwrite((char*)&dtemp, sizeof(float), 1, fwdcd);
	out_integer = 0;
	fwrite((char*)&out_integer, sizeof(float), 1, fwdcd);

	dtemp = yboxL;
	fwrite((char*)&dtemp, sizeof(float), 1, fwdcd);
	out_integer = 0;
	fwrite((char*)&out_integer, sizeof(float), 1, fwdcd);


	dtemp = zboxL;
	fwrite((char*)&dtemp, sizeof(float), 1, fwdcd);
	out_integer = 0;
	fwrite((char*)&out_integer, sizeof(float), 1, fwdcd);

	out_integer = 48;
	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	out_integer = 4 * N_atom;
	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	float ftemp;
	for (int i = 1; i <= N_atom; i++)
	{
		ftemp = (float)(SIGMA * r[i][0]);
		fwrite((char*)&ftemp, sizeof(float), 1, fwdcd);
	}

	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	for (int i = 1; i <= N_atom; i++)
	{
		ftemp = (float)(SIGMA * r[i][1]);
		fwrite((char*)&ftemp, sizeof(float), 1, fwdcd);
	}

	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	for (int i = 1; i <= N_atom; i++)
	{
		ftemp = (float)(SIGMA * r[i][2]);
		fwrite((char*)&ftemp, sizeof(float), 1, fwdcd);
	}

	fwrite((char*)&out_integer, sizeof(int), 1, fwdcd);

	/* update the DCD header information */
	int curframe = 0, curstep = 0, nsav = 0;

	fseek(fwdcd, 8L, SEEK_SET);
	fread(&curframe, sizeof(int), 1, fwdcd);
	curframe++;
	fseek(fwdcd, 20L, SEEK_SET);
	fwrite(&nsav, sizeof(int), 1, fwdcd);

	fseek(fwdcd, 16L, SEEK_SET);
	fwrite(&curstep, sizeof(int), 1, fwdcd);

	curstep += nsav;

	fseek(fwdcd, 8L, SEEK_SET);
	fwrite(&curframe, sizeof(int), 1, fwdcd);
	fseek(fwdcd, 20L, SEEK_SET);
	fwrite(&curstep, sizeof(int), 1, fwdcd);

	fseek(fwdcd, 0, SEEK_END);
}

