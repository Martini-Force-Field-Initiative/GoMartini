#ifndef contact_map_h
#define contact_map_h

#include <stddef.h>
#include <stdbool.h>

typedef struct {
	float x,y,z;
} float3;

typedef struct {
	char  field[7],name[5],resName[5],altLoc[2],chainID[2],iCode[2],element[3],charge[3];
	int   serial,resSeq; 
	float x,y,z,occupancy,tempFactor;
	int model;
} atom_pdb_str;

typedef struct {
	int natoms,nresidues;
	atom_pdb_str atoms[1];
} pdb_str;

typedef struct {
	int n; 			// NUMBER OF ATOMS IN RESIDUE
	float x,y,z;		// COORDINATES OF CENTER OF MASS
	float rad;		// RESIDUE RADIUS
	char name[4];		// RESIDUE 3 LETTER NAME
	int seq;		// RESIDUE NUMBER IN PDB FILE
	char chain[1];		// CHAIN ID
	int model;		// MODEL IN PDB FILE
	atom_pdb_str * a;
} residue_str;

typedef struct {
	int atype; 		// TYPES OF ATOMS
	int nb; 		// NUMBER OF ATOMS IN RESIDUES
	float vrad; 		// VAN DER WAALS RADII
	int keyresidue;		// RESIDUE KEY FROM ENUM PROTEIN MAP
	int keyatom;		// ATOM KEY FROM ENUM PROTEIN MAP
} atomaux_str;

typedef struct {
	float x,y,z;		// COORDINATES OF SURFACE POINTS
	float d;		// DISTANCE OF CLOSEST CENTER
	int i,j;		// I and J OF CLOSEST CENTER
} surface_str;

extern bool protein_map(atom_pdb_str *atom, atomaux_str *vdw);
extern void read_pdb (char *name, pdb_str *pdb);
extern float ATOMTYPE(unsigned int i,unsigned int j);
extern int BONDTYPE(int i, int j);

#endif
