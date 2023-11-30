#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "contact_map.h"

static int subStr(const char *src, char *dest, int begin, int end) {
	begin--; end--;
	memset(dest,'\0',end-begin+2);
	assert(begin >= 0 && end >= begin);
	while (src[begin] == ' ' && begin < end   ) begin++;
	while (src[end]   == ' ' && end   > begin ) end--;
	int length = end - begin + 1;
	strncpy(dest, src + begin, length);
	//printf("len %4d | %4d -> %4d | %s\n",length,begin+1,end+1,dest);
	//dest[length] = 0;
	return length;
}

static int getInt (const char *line, int begin, int end) {
	char str[81];
	subStr(line,str,begin,end);
	return atoi(str);
}

static double getdouble(const char *line, int begin, int end) {
	char str[81];
	subStr(line,str,begin,end);
	return atof(str);
}

static int getStr(const char *line, char str[], int begin, int end) {
	return subStr(line,str,begin,end);
}

extern void print_atom (atom_pdb_str *a) {
	printf("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n", a->field,a->serial,a->name,a->altLoc,a->resName,a->chainID,a->resSeq,a->iCode,a->x,a->y,a->z,a->occupancy,a->tempFactor,a->element,a->charge);
}

static atom_pdb_str read_atom (char *line) {
	atom_pdb_str a;
	getStr(line,a.field,1,6);
	a.serial = getInt(line,7,11);
	getStr(line,a.name,13,16);
	getStr(line,a.altLoc,17,17);
	getStr(line,a.resName,18,20);
	getStr(line,a.chainID,22,22);
	if (strncmp(a.chainID," ",1) == 0) {
		//printf("ChainID not found, assigning it to chain Z\n");
		strncpy(a.chainID,"Z\0",1);
	}
	a.resSeq = getInt(line,23,26);
	getStr(line,a.iCode,27,27);
	a.x = getdouble(line,31,38);
	a.y = getdouble(line,39,46);
	a.z = getdouble(line,47,54);
	a.occupancy = getdouble(line,55,60);
	a.tempFactor = getdouble(line,61,66);
	getStr(line,a.element,77,78);
	getStr(line,a.charge,79,80);
	//printf("%s",line);
	//print_atom(&a);
	return a;
}

extern void read_pdb (char *name, pdb_str *pdb) {
	FILE *pFile = fopen(name,"r");
	if(pFile == NULL) {
		fprintf(stderr,"Error while opening PDB file.");
		assert(false);
	}
	char line[120];
	int resSeq = -1;
	int model = 0;
	while (fgets(line,120,pFile) != NULL) {
		char tag[6];
		sscanf(line,"%6s",tag);
		if (strncmp(tag,"MODEL",5) == 0) {
			model = getInt(line,11,14);
			printf("\tReading PDB model %4d\n",model); fflush(stdout);
		}
		if (strncmp(tag,"ATOM",4) == 0) {
			atom_pdb_str atom = read_atom(line);
			if (resSeq != atom.resSeq) {
				resSeq = atom.resSeq;
				pdb->nresidues++;
			}
			atom.model = model;
			//print_atom(res);
			//printf("natoms: %d\n",pdb->natoms);
			atom_pdb_str *atm = pdb->atoms + pdb->natoms;
			memcpy(atm,&atom,sizeof(atom_pdb_str));
			pdb->natoms++;
		}
	}
	fclose(pFile);
}
