#!/usr/bin/env python3

import sys
import argparse
import subprocess
import numpy as np


def user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', help='File containing the coarse-grained structure of the protein in pdb format.')
    parser.add_argument('-f', help='File containing the contact analysis of the (atomistic) protein structure obtained from the webserver http://info.ifpan.edu.pl/~rcsu/rcsu/index.html.')
    parser.add_argument('--moltype', default='molecule_0', help='Molecule name used as prefix in your output file names and the virtual bead names (default: molecule_0). If you will combine your Go-like model with a coarse-grained protein generated with martinize2, you must use the same name as specified with the --govs-moltype flag of martinize2!')
    parser.add_argument('--go_eps', type=float, default=0.0, help='Dissociation energy [kJ/mol] of the Lennard-Jones potential used in the Go-like model (default: 0.0).')
    parser.add_argument('--cutoff_short', type=float, default=0.3, help='Lower cutoff distance [nm]: contacts with a shorter distance than cutoff_short are not included in the Go-like interactions (default: 0.3).')
    parser.add_argument('--cutoff_long', type=float, default=1.1, help='Upper cutoff distance [nm]: contacts with a longer distance than cutoff_long are not included in the Go-like interactions (default: 1.1).')
    parser.add_argument('--Natoms', type=int, help='Number of coarse-grained beads in the protein excluding the virtual Go beads.')
    parser.add_argument('--missres', type=int, default=0, help='Number of missing residues at the beginning of the atomistic pdb structure which is needed if the numbering of the coarse-grained structure starts at 1 (default: 0).')
    parser.add_argument('--idp_eps', type=float, default=0.0, help='Dissociation energy [kJ/mol] of the Lennard-Jones potential used to increase the water interaction of the protein backbone of IDPs via the virtual backbone beads (default: 0.0).')
    parser.add_argument('--idp_start', type=int, default=0, help='First resid of the IDP fragment of the protein (default: first residue of the protein).')
    parser.add_argument('--idp_end', type=int, default=0, help='Last resid of the IDP fragment of the protein (default: last residue of the protein).')
    parser.add_argument('--idp_auto', type=bool, default=False, help='Automatically apply water-BB interaction adjustments based on secondary structure.')
    parser.add_argument('--itp', help='File containing the Martini 3 protein itp.')
    args = parser.parse_args()
    return args


def get_settings():
    # names for temporary files:
    file_BB = 'BB.pdb'
    file_OV = 'OV.map'
    file_rCSU = 'rCSU.map'

    # some other rudimentary variables
    header_lines = 0
    seqDist = 4         # minimal distance in the sequence to add a elastic bond (ElNedyn=3 [Perriole2009]; Go=4 [Poma2017])
    cols = [5, 9, 10]   # colums of interest in the OV and rCSU contact map file
    missAt = 0          # number of missing atoms at the beginning of pdb structure
                        # (this has to result in the correct atom number when added to "k_at" compared to the .itp file)
    c6c12 = 0           # if set to 1, the C6 and C12 term are expected in the .itp file; if set to 0, sigma and go_eps are used
    #idp_sig = 0.465	# sigma for the additional LJ interaction between virtual BB beads and W for IDPs
    #deprecated in favor of having sig per residue
    return (file_BB, file_OV, file_rCSU, header_lines,
             seqDist, cols, missAt, c6c12)#, idp_sig)


def read_pdb(struct_pdb, file_BB, header_lines):
    # preparation of temporary files for reading
    subprocess.call("grep 'BB' " + struct_pdb + " > " + file_BB, shell=True)
    subprocess.call("echo '' >> " + file_BB, shell=True)

    # read coarse-grained BB bead positions
    with open(file_BB,'r') as fid:
        dat = fid.readlines()
    dat = dat[header_lines:-1]
    print('Number of coarse-grained BB beads in your protein: ' 
          + str(len(dat)))

    indBB = []
    nameAA = []
    for k in range(0, len(dat)):
        tmp = dat[k]
        tmp = tmp.split()
        indBB.append([ int(tmp[1]), float(tmp[6]), 
                      float(tmp[7]), float(tmp[8]) ])
        nameAA.append(tmp[3])
    indBB = np.array(indBB)
    print(nameAA)
    return indBB, nameAA


def read_contmap(file_contacts, file_OV, file_rCSU, header_lines, cols):
    # preparation of temporary files for reading
    subprocess.call("grep '1 [01] [01] [01]' " 
                    + file_contacts + " > " + file_OV, shell=True)
    subprocess.call("echo '' >> " + file_OV, shell=True)
    subprocess.call("grep '0 [01] [01] 1' " + file_contacts 
                    + " > " + file_rCSU, shell=True)
    subprocess.call("echo '' >> " + file_rCSU, shell=True)

    # read OV contact map
    with open(file_OV,'r') as fid:
        dat = fid.readlines()
    dat = dat[header_lines:-1]
    print('Number of contacts read from your OV contact map file: ' 
          + str(len(dat)))

    map_OVrCSU = []
    row = []
    for k in range(0, len(dat)):
        tmp = dat[k]
        tmp = tmp.replace('\t',' ')
        tmp = tmp.split()
        for l in cols:
            row.append(float(tmp[l]))
        map_OVrCSU.append(row)
        row = []

    # read rCSU contact map
    with open(file_rCSU,'r') as fid:
        dat = fid.readlines()
    dat = dat[header_lines:-1]
    print('Number of contacts read from your rCSU contact map file: ' 
          + str(len(dat)))

    for k in range(0, len(dat)):
        tmp = dat[k]
        tmp = tmp.replace('\t',' ')
        tmp = tmp.split()
        for l in cols:
            row.append(float(tmp[l]))
        map_OVrCSU.append(row)
        row = []

    return map_OVrCSU


def get_go(indBB, nameAA, map_OVrCSU, cutoff_short, 
           cutoff_long, go_eps, seqDist, missRes):
    # calculate the distances based on the coordinates of the CG BB bead
    for k in range(0, len(map_OVrCSU)):
        dist_vec = (indBB[ int(map_OVrCSU[k][1])-missRes-1 ,1:4] - 
                    indBB[ int(map_OVrCSU[k][0])-missRes-1 ,1:4])
        map_OVrCSU[k][2] = np.linalg.norm(dist_vec) / 10     # [Ang] to [nm]

    pairs = []
    for k in range(0, len(map_OVrCSU)):
        if ((map_OVrCSU[k][2] > cutoff_short) and 
           (map_OVrCSU[k][2] < cutoff_long) and 
           ( abs(map_OVrCSU[k][1]-map_OVrCSU[k][0]) >= seqDist )):
            # parameters for LJ potential
            # calc sigma for the LJ potential in [nm]
            sigma = map_OVrCSU[k][2] / 1.12246204830   
            Vii = 4.0 * pow(sigma,6) * go_eps
            Wii = 4.0 * pow(sigma,12) * go_eps
            pairs.append([indBB[ int(map_OVrCSU[k][0])-missRes-1 ,0],
                          indBB[ int(map_OVrCSU[k][1])-missRes-1 ,0], Vii,
                          Wii, map_OVrCSU[k][0], map_OVrCSU[k][1],
                          map_OVrCSU[k][2], sigma])
            Vii = []
            Wii = []
        elif map_OVrCSU[k][2] > cutoff_long:
            print('This contact is excluded due to distance > cutoff_long: '
                   + str(map_OVrCSU[k]))
        elif map_OVrCSU[k][2] < cutoff_short:
            print('This contact is excluded due to distance < cutoff_short: '
                   + str(map_OVrCSU[k]))
        elif abs(map_OVrCSU[k][1]-map_OVrCSU[k][0]) < seqDist:
            print('This contact is excluded because the AA have less than '
                   + str(seqDist-1) + ' other AA between each other: '
                   + str(map_OVrCSU[k]))

    sym_pairs = []
    # count contacts only once; 
    # exclude asymmetric rCSU contacts (cf. doi 10.1063/1.4929599)
    for k in range(0, len(pairs)):
        if pairs[k][0] < pairs[k][1]:
            for l in range(k+1, len(pairs)):
                if (pairs[l][0] == pairs[k][1]) and (pairs[l][1] == pairs[k][0]):
                    sym_pairs.append(pairs[k])

    print ('- - - -')
    print ('These results exclude the contacts with distances higher than the cutoff_long (' + str(cutoff_long) + ' nm), shorter than the cutoff_short (' + 
            str(cutoff_short) + ' nm), or where the AA have less than ' + str(seqDist-1) + ' other AA between each other:')
    print ('Sum of symmetric (doubly counted) and asymmetric OV + rCSU contacts: ' + str(len(pairs)))
    print ('Only symmetric OV + rCSU contacts (singly counted):' + str(len(sym_pairs)))

    return sym_pairs


def write_genfiles(file_pref, sym_pairs, missAt, indBB, missRes):
    # writes the files required for both cases - Go-like model and IDP solubility
    # write supplementary file: BB virtual particle definitions for martini.itp
    with open(file_pref + '_BB-part-def_VirtGoSites.itp','w') as f:
        f.write('; protein BB virtual particles \n')
        for k in range(0, len(indBB)):
            s2print = "%s_%s 0.0 0.000 A 0.0 0.0 \n" % (file_pref, str(k+1 +missRes)) # residue index adapted due to missing residues
            f.write(s2print)
    subprocess.call("echo '#include \"" + file_pref + "_BB-part-def_VirtGoSites.itp\"' >> BB-part-def_VirtGoSites.itp ", shell=True)

    # write supplementary file: exclusions for protein.itp
    with open(file_pref + '_exclusions_VirtGoSites.itp','w') as f:
        f.write(';[ exclusions ] \n')
        f.write('; OV + symmetric rCSU contacts \n')
        for k in range(0, len(sym_pairs)):
            s2print = " %s  %s  \t ;  %s  %s \n" % (str(int(sym_pairs[k][0]) +missAt), str(int(sym_pairs[k][1]) +missAt), 
                                                    str(int(sym_pairs[k][4] -missRes)), str(int(sym_pairs[k][5] -missRes))) 
                                        # atom index and residue index adapted due to missing residues
            f.write(s2print)


def write_gofiles(file_pref, sym_pairs, missAt, indBB, missRes, Natoms, go_eps, c6c12):
    # writes the files required solely for the Go-like model
    # write the interaction table for the Go-like bonds
    with open(file_pref + '_go-table_VirtGoSites.itp','w') as f:
        f.write('; OV + symmetric rCSU contacts \n')
        if (c6c12 == 1):
            for k in range(0, len(sym_pairs)):
                # to write the LJ potential itp:
                s2print = " %s_%s  %s_%s    1  %.10f  %.10f  ;  %s  %s  %.3f \n" % (file_pref, str(int(sym_pairs[k][4])), file_pref, str(int(sym_pairs[k][5])), 
                                                                          sym_pairs[k][2], sym_pairs[k][3], str(int(sym_pairs[k][0]) +missAt), 
                                                                          str(int(sym_pairs[k][1]) +missAt), sym_pairs[k][6]) 
                                                            # atom index and residue index adapted due to missing residues
                f.write(s2print)
        else:
            for k in range(0, len(sym_pairs)):
                # to write the LJ potential itp:
                s2print = " %s_%s  %s_%s    1  %.10f  %.10f  ;  %s  %s  %.3f \n" % (file_pref, str(int(sym_pairs[k][4])), file_pref, str(int(sym_pairs[k][5])), 
                                                                          sym_pairs[k][7], go_eps, str(int(sym_pairs[k][0]) +missAt), 
                                                                          str(int(sym_pairs[k][1]) +missAt), sym_pairs[k][6]) 
                                                            # atom index and residue index adapted due to missing residues
                f.write(s2print)

    subprocess.call("echo '#include \"" + file_pref + "_go-table_VirtGoSites.itp\"' >> go-table_VirtGoSites.itp ", shell=True)

    # write supplementary file: Go-like bonds as harmonic bonds for visulization of the protein
    with open(file_pref + '_go4view_harm.itp','w') as f:
        f.write('; Go bonds as harmonic bonds between the virtual particles: \n')
        f.write('; OV + symmetric rCSU contacts \n')
        for k in range(0, len(sym_pairs)):
            # to write the harmonic bonds itp:
            s2print = " %s  %s  1  %s  1250  ; %s_%s  %s_%s \n" % (str(int(sym_pairs[k][4] -missRes +Natoms)), str(int(sym_pairs[k][5] -missRes +Natoms)), str(round(sym_pairs[k][6],3)), 
                                                                 file_pref, str(int(sym_pairs[k][4] -missRes)), file_pref, str(int(sym_pairs[k][5] -missRes)))
                                                               # the bonds are added between the virtual particles
            f.write(s2print)
        for k in range(0, len(indBB)):
            # print (np.sum(np.array(sym_pairs)[:,4]==k+1) + np.sum(np.array(sym_pairs)[:,5]==k+1))
            if (np.sum(np.array(sym_pairs)[:,4]==k+1) + np.sum(np.array(sym_pairs)[:,5]==k+1)) == 0:
                s2print = " %s  %s  1  1.  1     ; %s_%s  %s_%s --> added for vmd \n" % (str(int(k+1 +Natoms)), str(int(k +Natoms)), file_pref, str(k+1), file_pref, str(k))
                f.write(s2print)


def write_idpfiles(file_pref, indBB, missRes, idp_sig, idp_eps, idp_start, idp_end, c6c12):
    # writes the files required solely for the IDP-solubility model
    # write the interaction table for the virtual Go beads with water
    with open(file_pref + '_go-table_IDPsolubility.itp','w') as f:
        f.write('; additional Lennard-Jones interaction between virtual BB bead and W\n')

        # adapt residue index due to missing residues if no custom residue interval is specified
        if idp_start == 0:
            idp_start = missRes + 1
        if idp_end == 0:
            idp_end = missRes + len(indBB)
        if idp_start > idp_end:
            sys.exit('Error: \nYour IDP sequence has a larger ending than starting residue index. Please check your input! \nWill terminate now.')

        if (c6c12 == 1):
            for k in range(idp_start, idp_end+1):
                # to write the LJ potential itp:
                Vii = 4.0 * pow(idp_sig[k-1],6) * idp_eps
                Wii = 4.0 * pow(idp_sig[k-1],12) * idp_eps

                s2print = " %s_%s  %s    1  %.10f  %.10f ; \n" % (file_pref, str(int(k)), "W", Vii, Wii) 
                f.write(s2print)
        else:
            for k in range(idp_start, idp_end+1):
                # to write the LJ potential itp:
                s2print = " %s_%s  %s    1  %.10f  %.10f ; \n" % (file_pref, str(int(k)), "W", idp_sig[k-1], idp_eps)
                f.write(s2print)

    subprocess.call("echo '#include \"" + file_pref + "_go-table_IDPsolubility.itp\"' >> go-table_VirtGoSites.itp ", shell=True)


def write_autoidp(file_pref, indBB, missRes, idp_sig, itp):
    # writes the files required solely for the IDP-solubility model
    # write the interaction table for the virtual Go beads with water
    with open(file_pref + '_go-table_IDPsolubility.itp','w') as f:
        f.write('; additional Lennard-Jones interaction between virtual BB bead and W\n')

        # adapt residue index due to missing residues if no custom residue interval is specified
        idp_start = missRes + 1
        idp_end = missRes + len(indBB)
        
        ss = get_ss(itp)

        eps=[]
        for aa in ss:
            if aa == 'C':
                eps.append(0.5)
            elif aa == 'H':
                eps.append(-0.5)
            else:
                eps.append(0.0)

        for k in range(idp_start, idp_end+1):
            # to write the LJ potential itp:
            s2print = " %s_%s  %s    1  %.10f  %.10f ; \n" % (file_pref, str(int(k)), "W", idp_sig[k-1], eps[k-1])
            f.write(s2print)

    subprocess.call("echo '#include \"" + file_pref + "_go-table_IDPsolubility.itp\"' >> go-table_VirtGoSites.itp ", shell=True)


def get_idp_sig(nameAA):
    idp_sig=[]
    for aa in nameAA:
        if aa in ['GLY', 'ALA', 'VAL', 'PRO']:
            idp_sig.append(0.430)
        else:
            idp_sig.append(0.470)
    return idp_sig


def get_ss(itp):
    with open(itp,'r') as fid:
        dat = fid.readlines()
    dat = dat[5][2:-1]
    return dat


def main():
    args = user_input()
    #file_BB, file_OV, file_rCSU, header_lines, seqDist, cols, missAt, c6c12, idp_sig = get_settings()
    file_BB, file_OV, file_rCSU, header_lines, seqDist, cols, missAt, c6c12 = get_settings()
    indBB, nameAA = read_pdb(args.s, file_BB, header_lines)
    idp_sig = get_idp_sig(nameAA)

    if args.go_eps != 0.0:
        map_OVrCSU = read_contmap(args.f, file_OV, file_rCSU, header_lines, cols)
        sym_pairs = get_go(indBB, nameAA, map_OVrCSU, args.cutoff_short, args.cutoff_long, args.go_eps, seqDist, args.missres)
        write_gofiles(args.moltype, sym_pairs, missAt, indBB, args.missres, args.Natoms, args.go_eps, c6c12)
        print('All symmetric OV and rCSU contacts written! Have fun!')
    
    if args.idp_auto:
        sym_pairs = []
        write_autoidp(args.moltype, indBB, args.missres, idp_sig, args.itp)
        print('Additional BB-W interactions written for coils and helices! Have fun!')

    elif args.idp_eps != 0.0:
        sym_pairs = []
        write_idpfiles(args.moltype, indBB, args.missres, idp_sig, args.idp_eps, args.idp_start, args.idp_end, c6c12)
        print('Additional BB-W interactions written! Have fun!')

    write_genfiles(args.moltype, sym_pairs, missAt, indBB, args.missres)

if __name__ == '__main__':
    main()
