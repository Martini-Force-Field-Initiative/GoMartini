#include <assert.h>
#include "contact_map.h"

// REFERENCE: J. Chem. Phys. 143, 243105 (2015); https://doi.org/10.1063/1.4929599

extern residue_str *residues;
extern pdb_str *pdb;
extern atomaux_str *atomauxs;

const int MAXATOMTYPE = 10;
const int MAXBONDTYPE = 6;

float ATOMTYPE(unsigned int i,unsigned int j) {
	// Atom classes:
	// 1  -- hydrophilic - N and O atoms which can be both donor and acceptor of hydrogen bond
	// 2  -- hydrogen bond acceptor
	// 3  -- hydrogen bond donor,
	// 4  -- hydrophobic - all C atoms not connected to N or O and not in aromatic ring
	// 5  -- aromatic - all C atoms in aromatic rings
	// 6  -- neutral - S atoms and C atoms which
	// 7  -- neutral-donor - C atoms
	// 8  -- neutral-acceptor - C atoms
	// 9  -- positively charged - atoms with positive charge
	// 10 -- negatively charged - atoms with negative charge
	residue_str * res = residues+i;
	atom_pdb_str *a = res->a;
	atom_pdb_str * atom = a+j;
	ptrdiff_t k = atom - pdb->atoms;
	atomaux_str *aux = atomauxs + k;
	int atype = aux->atype;
	assert(atype >= 0 && atype <= MAXATOMTYPE );
	return atype;
}

int BONDTYPE(int i, int j) {
	assert(i >= 1 && i <= MAXATOMTYPE);
	assert(j >= 1 && j <= MAXATOMTYPE);
	i -= 1;
	j -= 1;
	// BOND TYPE
	//Types of contacts:
	//HB -- 1 -- hydrogen-bond
	//PH -- 2 -- hydrophobic
	//AR -- 3 -- aromatic - contacts between aromatic rings
	//IB -- 4 -- ionic bridge - contacts created by two atoms with different charges
	//DC -- 5 -- destabilizing contact - contacts which are in general repulsive
	//OT -- 6 -- denotes negligible other contacts.
	// 1-HB,2-PH,3-AR,4-IP,5-DC,6-OT
	const int BONDTYPE[10][10] = {
	{1,1,1,5,5,6,6,6,1,1},
	{1,5,1,5,5,6,6,6,1,5},
	{1,1,5,5,5,6,6,6,5,1},
	{5,5,5,2,2,6,6,6,5,5},
	{5,5,5,2,3,6,6,6,5,5},
	{6,6,6,6,6,6,6,6,6,6},
	{6,6,6,6,6,6,6,6,6,6},
	{6,6,6,6,6,6,6,6,6,6},
	{1,1,5,5,5,6,6,6,5,4},
	{1,5,1,5,5,6,6,6,4,5}
	};
	return BONDTYPE[i][j];
};

