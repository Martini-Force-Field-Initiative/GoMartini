#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include "contact_map.h"

// GLOBAL VARIABLES
pdb_str        *pdb      = NULL;
atomaux_str    *atomauxs = NULL;
residue_str    *residues = NULL;
surface_str    *surface  = NULL;

#define string_zeroing(T, a, n) do { \
	T *b = (a);                  \
	size_t m = (n);              \
	for (;m>0;--m,++b) *b=(T){0};\
} while (0)

bool string_ncopy(char * destination, const char *source, size_t num) {
	if (destination == NULL) return false;
	char *ptr = destination;
	while(*source && num--) { *ptr = *source; ptr++; source++; }
	return true;
}

void make_surface(surface_str *surf, float3 c,int fiba,int fibb,float vrad) {
	float x=c.x,y=c.y,z=c.z;
	int phi_aux = 0;
	for(int k=0; k<fibb; k++){
		surface_str * s = surf + k;

		// INITIALIZE
		s->i = -1;
	      	s->j = -1;
		s->d = FLT_MAX;

		phi_aux += fiba;
		if(phi_aux > fibb) phi_aux = phi_aux - fibb;

		// POINT
		float theta = acos(1.0-2.0*k/fibb);
		float phi   = 2.0*M_PI*phi_aux/fibb;
		s->x = x + vrad*sin(theta)*cos(phi);
		s->y = y + vrad*sin(theta)*sin(phi);
		s->z = z + vrad*cos(theta);
	}
}

float dist(float3 c1,float3 c2) {
	float x1=c1.x,y1=c1.y,z1=c1.z;
	float x2=c2.x,y2=c2.y,z2=c2.z;
	float dx=(x2-x1),dy=(y2-y1),dz=(z2-z1);
	return sqrtf(dx*dx + dy*dy + dz*dz);
}


// USEFUL GET FUNCTIONS
float3 SURFACE_COORDINATES(int k) {
	surface_str *s = surface + k;
	float3 r = {s->x,s->y,s->z};
	return r;
}
float3 RESIDUE_CENTER_OF_MASS(int i) {
	// RESIDUE CENTER OF MASS
	residue_str *res = residues+i;
	float3 rc = {res->x,res->y,res->z};
	return rc;
}
float RESIDUE_RADIUS(int i) {
	residue_str *res = residues+i;
	return res->rad;
}
float3 ATOM_COORDINATES(int i,int j) {
	residue_str *res = residues+i;
	atom_pdb_str *a = res->a;
	assert(j>=0 && j<res->n);
	a = a+j;
	float3 c = {a->x,a->y,a->z};
	return c;
}
int ATOM_INDEX(int i1,int j1) {
	residue_str * res1 = residues+i1;
	atom_pdb_str *a1 = res1->a;
	atom_pdb_str * atom1 = a1+j1;
	ptrdiff_t k1 = atom1 - pdb->atoms;
	return k1;
}
float vdW_RADIUS(int i1,int j1) {
	residue_str * res1 = residues+i1;
	atom_pdb_str *a1 = res1->a;
	atom_pdb_str * atom1 = a1+j1;
	ptrdiff_t k1 = atom1 - pdb->atoms;
	atomaux_str *aux1 = atomauxs + k1;
	return aux1->vrad;
}

int BASES_IN_RESIDUE(int i1) {
	residue_str * res1=residues + i1;
	int n1 = res1->n;
	return n1;
}

void init_atomauxs_and_residues() {
	int resCount = -1;
	int resSeq = -1;
	for (int k=0;k<pdb->natoms;k++) {
		atom_pdb_str * atom = pdb->atoms + k;
		// new residue
		if (resSeq != atom->resSeq) {
			resSeq = atom->resSeq;
			resCount++;
			residue_str * res = residues + resCount;
			res->a = atom;
			strncpy(res->name,atom->resName,4);
			res->seq = atom->resSeq;
			strncpy(res->chain,atom->chainID,1);
			res->model = atom->model;
		}

		// init aux
		atomaux_str * aux = atomauxs + k;
		protein_map(atom,aux);

		// print info
		//printf("%4s %4d | ",atom->resName,aux->keyresidue);
		//printf("%4s %4d | ",atom->name,aux->keyatom);
		//printf("%2d %4.2f %2d | ",aux->nb,aux->vrad,aux->atype);
		//printf("%5d | ",k+1); print_atom(atom);

		// init residue
		residue_str * res = residues + resCount;
		res->n++;
		res->x += atom->x;
		res->y += atom->y;
		res->z += atom->z;
	}
	assert(resCount+1 == pdb->nresidues);
	// CENTROID
	for (int k=0; k<pdb->nresidues; k++) {
		residue_str * res = residues + k;
		res->x /= res->n;
		res->y /= res->n;
		res->z /= res->n;
	}
}

float DISTANCE_C_ALPHA(int i1,int i2) {
	// C-ALPHA DISTANCES
	int j1=1,j2=1; // C-ALPHA IS GIVEN AS THE SECOND ATOM OF A RESIDUES IN PDB
	return dist(ATOM_COORDINATES(i1,j1),ATOM_COORDINATES(i2,j2));
}

/*
void printA(int i1, int j1, int i2, int j2) {
	int k1 = ATOM_INDEX(i1,j1);
	int k2 = ATOM_INDEX(i2,j2);
	float d = dist(ATOM_COORDINATES(i1,j1),ATOM_COORDINATES(i2,j2));
	int t1 = atomauxs[k1].atype, t2 = atomauxs[k2].atype;
	printf("A%5d  %3s   %1s %-3s %2d"
               " %5d  %3s   %1s %-3s %2d"
               "    %2d   %8.4f\n"
        ,pdb->atoms[k1].resSeq,pdb->atoms[k1].resName,pdb->atoms[k1].chainID,pdb->atoms[k1].name,t1,
	 pdb->atoms[k2].resSeq,pdb->atoms[k2].resName,pdb->atoms[k2].chainID,pdb->atoms[k2].name,t2,
	BONDTYPE(t1,t2), d
);
}
*/
/*
printf("\n\n"
"Atom  - Atom contacts\n"
"I1,I2 - id of residues in contact\n"
"A1,A2 - names of atoms in contact\n"
"C     - chain\n"
"T     - id of atom class\n"
"I     - id of interaction type\n"
"Surf  - surface of interaction in CSU algorithm\n"
"S0    - whole surface of overlap in CSU algorithm\n"
"Cont  - number of contacts between atoms\n"
"    I1  AA    C A1   T    I2  AA    C A2   T     I   DISTANCE       Surf      S0    Cont\n"
"========================================================================================\n"
);
*/
int main (int argc,char **argv) {
	printf(
"                         CONTACT MAPS FROM PDB FILES                          \n"
"                                                                              \n"
" This software is written by:                                                 \n"
"       Rodrigo Azevedo Moreira da Silva                                       \n"
"                                                                              \n"
" Copyright (c) 2020 - IPPT-PAN                                                \n"
"       Institute of Fundamental Techonological Research                       \n"
"       Polish Academy of Sciences                                             \n"
" MIT LICENSE, check out LICENSE for more informations.                        \n"
"                                                                              \n"
        );


	size_t size_pdb = sizeof(pdb_str) + 100000*sizeof(atom_pdb_str);

	printf("Reading file:    %s\n",argv[1]);

	// ALOCATING MEMORY TO PDB FILE

	pdb = (pdb_str *)malloc(size_pdb);
	assert(pdb != NULL);
	memset(pdb,0,size_pdb);

	// READ PDB "ATOM" INFO
	read_pdb(argv[1],pdb);
	int natoms = pdb->natoms;
	int nresidues = pdb->nresidues;
	printf("pdb natoms:      %d\n",natoms);
	printf("pdb nresidues:   %d\n",nresidues);

	printf("Memory usage:  %7.2f MB\n",(
		size_pdb +
		natoms*sizeof(atomaux_str) +
		nresidues*sizeof(residue_str) +
		nresidues*nresidues*sizeof(int) * 4
		)/(1024.0*1024.0));

	// AUXILIARY STRUCTURES
	atomauxs = (atomaux_str*)calloc(natoms,sizeof(atomaux_str));    assert(atomauxs != NULL);
	residues = (residue_str*)calloc(nresidues,sizeof(residue_str)); assert(residues != NULL);
	init_atomauxs_and_residues();

	// FIBONACCI
	int fib = 14, fiba = 0, fibb = 1;
	for (int f=0;f<fib;f++) { int fibc=fiba+fibb; fiba=fibb; fibb=fibc; }

	surface = (surface_str*)calloc(fibb,sizeof(surface_str));
	printf("Fibonacci grid:  %d\n",fibb);

	int * overlapcounter      = (int *)calloc(nresidues*nresidues,sizeof(int)); assert(overlapcounter);
	int * contactcounter      = (int *)calloc(nresidues*nresidues,sizeof(int)); assert(contactcounter);
	int * stabilizercounter   = (int *)calloc(nresidues*nresidues,sizeof(int)); assert(stabilizercounter);
	int * destabilizercounter = (int *)calloc(nresidues*nresidues,sizeof(int)); assert(destabilizercounter);
	//float scmap[pdb->nresidues][pdb->nresidues][MAXBONDTYPE];

	// MODEL VARIABLES
	float alpha        = 1.24; // ENLARGMENT FACTOR TO ACCCOUNT ATTRACTION EFFECTS
	float water_radius = 2.80; // RADIUS OF WATER MOLECULE
	printf("ALPHA:        %7.2f\n",alpha);
	printf("WATER_RADIUS: %7.2f\n",water_radius);

	// MAKE CONTACT LIST
	for(unsigned int i1=0; i1 < nresidues; i1++) {
	for(unsigned int j1=0; j1 < BASES_IN_RESIDUE(i1); j1++) {
		// INITIALIZE AND MAKE SURFACE
		make_surface(surface,ATOM_COORDINATES(i1,j1),fiba,fibb,vdW_RADIUS(i1,j1)+water_radius);

		// FIND NEIGHBORS OF i1 j1
		for(unsigned int i2=0; i2 < nresidues; i2++) {
		// CHECK RESIDUES CONTACT AND ON THE SAME MODEL
		if (
			(residues[i1].model == residues[i2].model) &&
			(dist(RESIDUE_CENTER_OF_MASS(i1),RESIDUE_CENTER_OF_MASS(i2)) <= 3.5*4) // roughly one aminoacid = 3.5A
		) {
		for(unsigned int j2=0; j2 < BASES_IN_RESIDUE(i2); j2++) {
		// CHECK NOT THE SAME CENTER
		if (!(i1==i2 && j1==j2)) {

			float distance = dist(ATOM_COORDINATES(i1,j1),ATOM_COORDINATES(i2,j2));

			// Enlarged overlap (OV) contact
      			if(distance <= ((vdW_RADIUS(i1,j1)+vdW_RADIUS(i2,j2))*alpha)) {
				overlapcounter[i1 + nresidues*i2] = 1;
			}

			// Contacts of Structural Units (CSU)
			if (distance <= vdW_RADIUS(i1,j1)+vdW_RADIUS(i2,j2)+water_radius) {
				// FIND I2 J2 CLOSEST TO A GIVEN POINT ON THE SURFACE
				for (int k=0; k<fibb; k++) {
				surface_str * s = surface + k;
				if( dist(SURFACE_COORDINATES(k),ATOM_COORDINATES(i2,j2)) < vdW_RADIUS(i2,j2)+water_radius
                                    && distance <= s->d) {	
      					s->d = distance;
      					s->i = i2;
      					s->j = j2;
				}}
			}

		}}}}

		// AREA OF SURFACE
		//float sarea = vdW_RADIUS(i1,j1)*vdW_RADIUS(i1,j1)*M_PI*4/fibb;

		for(int k=0; k<fibb; k++) {
			surface_str * s = surface + k;
			int i2 = s->i;
			int j2 = s->j;
			if( i2 >= 0 && j2 >= 0 ) {
				int at1 = ATOMTYPE(i1,j1);
				int at2 = ATOMTYPE(i2,j2);
				// CHECK ONLY MAPPED ATOM
				if ( at1>0 && at2>0 ) {
					contactcounter[i1 + nresidues*i2] += 1;
					//if ( i1 != i2 ) printA(i1,j1,i2,j2);
					int btype = BONDTYPE(at1,at2);
					//scmap[i1][i2][btype-1] += sarea;
					// COUNT STABILIZING   BONDS
					if(btype <= 4) stabilizercounter[i1 + nresidues*i2] += 1;
					// COUNT DESTABILIZING BONDS
					if(btype == 5) destabilizercounter[i1 + nresidues*i2] += 1;
				}
			}
		}
	}}

	printf("\n"
"Residue-Residue Contacts\n"
"\n"
"ID       - atom identification\n"
"I1,I2    - serial residue id\n"
"AA       - 3-letter code of aminoacid\n"
"C        - chain\n"
"I(PDB)   - residue number in PDB file\n"
"DCA      - distance between CA\n"
"CMs      - OV , CSU , oCSU , rCSU\n"
"           (CSU does not take into account chemical properties of atoms)\n"
"rCSU     - net contact from rCSU\n"
"Count    - number of contacts between residues\n"
"MODEL    - model number\n"
//"aSurf    - surface of attractive connections\n"
//"rSurf    - surface of repulsive connections\n"
//"nSurf    - surface of neutral connections\n"
"\n"
"      ID    I1  AA  C I(PDB)     I2  AA  C I(PDB)        DCA       CMs    rCSU   Count Model\n"
"============================================================================================\n"
//       
);
	int count = 0;
	for(int i1=0;i1 < nresidues;i1++) {
	for(int i2=0;i2 < nresidues;i2++) {
		int over = overlapcounter[i1 + nresidues*i2];
		int cont = contactcounter[i1 + nresidues*i2];
		int stab = stabilizercounter[i1 + nresidues*i2];
		int dest = destabilizercounter[i1 + nresidues*i2];
		int ocsu = stab;
		int rcsu = stab - dest;
		int mdl1 = residues[i1].model;
		int mdl2 = residues[i2].model;
		//float sum = 0.0;
		//for (int k=0;k<MAXBONDTYPE;k++) sum += scmap[i1][i2][k];
		if (mdl1 == mdl2 && i1 != i2 && (over > 0 || cont > 0)) {
			count++;
			printf("R %6d ",count);
			printf("%5d %4s %1s %4d",i1+1,residues[i1].name,residues[i1].chain,residues[i1].seq);
                        printf("    ");
			printf("%5d %4s %1s %4d",i2+1,residues[i2].name,residues[i2].chain,residues[i2].seq);
			printf("     %8.4f     ",DISTANCE_C_ALPHA(i1,i2));
			printf("%d %d %d %d",
				over,
				cont != 0 ? 1 : 0,
				ocsu != 0 ? 1 : 0,
				rcsu >  0 ? 1 : 0);
			printf("%6d  %6d %4d\n",rcsu,cont,residues[i1].model);
		}
	}}

	free(overlapcounter);
	free(contactcounter);
	free(stabilizercounter);
	free(destabilizercounter);
	free(pdb);
	free(atomauxs);
	free(residues);
	free(surface);
}













