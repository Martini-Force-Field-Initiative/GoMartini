# GōMartini — powered by create_goVirt

`create_goVirt` powers the deployment of the new virtual-site implementation of GōMartini which can be combined with the latest iteration of the [Martini](http://cgmartini.nl/) model, [Martini 3](https://doi.org/10.1038/s41592-021-01098-3). 

## Installation

`create_goVirt` requires python 3.6 or greater. It is distributed via PyPi, and can be installed using the pip command:
``````
pip install GoMartini
``````

This installs the last released version. You can update an existing installation by running `pip install -U GoMartini`. In some cases you may want to experiment with running the latest development version. You can install this version with the following command:
``````
pip install git+https://github.com/Martini-Force-Field-Initiative/GoMartini
``````

Note that development versions, may contain bugs that cause it to produce incorrect topologies. Check the produced output carefully!

The behavior of the pip command can vary depending of the specificity of your python installation. See the documentation on installing a python package to learn more.

## Basic usage

### 1. Martinize

To set up a GōMartini model, the atomistic protein structure must first be martinized using [martinize2](https://github.com/marrink-lab/vermouth-martinize). 

The `--govs-include` martinize2 flag prepares a `protein.itp` with the `include` statements required to incorporate the Gō-like model generated by `create_goVirt`. A second flag (`--govs-moltype`) specifies the molecule name. This name is used as prefix for the files containing the Gō-like model (which will be generated in a following step by `create_goVirt`). 

A typical Martinize command for generating a Gō-ready `protein.itp` may look something like this:

````
martinize2 -f pdb.pdb -x cg.pdb -o topol.top -ff martini3001 -dssp dssp -scfix -govs-include -govs-moltype molecule_0 
````


### 2. Create the contact map

Gō-like contacts are defined by two criteria: a residue overlap criterion, and a [restricted chemical structural units (rCSU) criterion](https://doi.org/10.1021/acs.jctc.6b00986). This rCSU criterion is dependent on a contact map, which will be required in the following step. 

This rCSU contact map can be obtained from a [webserver](http://pomalab.ippt.pan.pl/GoContactMap/). An alternative older version of the [webserver](http://info.ifpan.edu.pl/~rcsu/rcsu/index.html) is also available. Input your original atomistic `protein.pdb` file and use the default values for the `radii` and `Fibonacci number`. Note that the requirements for the `pdb` file uploaded to the webserver are quite strict. Thus do carefully check if your contact map is meaningful before using it in the next step. In particular, the table listing the “Residue residue contacts” will be used.

Alternatively, a locally executable version of the webserver is also available. The source files for the `contact_map` executable can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3817447.svg)](https://doi.org/10.5281/zenodo.3817447). If you use the locally executable version of the contact map, please cite:

> Koehler, M., Ray, A., Moreira, R.A. et al. Molecular insights into receptor binding energetics and neutralization of SARS-CoV-2 variants. Nat Commun 12, 6977 (2021). https://doi.org/10.1038/s41467-021-27325-1


### 3. create_goVirt

The `cg.pdb` structure and `.itp` obtained from martinize2 (see Step 1) are required for the next step, together with a contact map (`.map`) of the atomistic protein structure (see Step 2). 

`create_goVirt` generates all additional files required for GōMartini. 

#### 3.1. Input

In its most basic application, `create_goVirt` requires the CG protein structure of the protein in pdb format (`-s`), the number of CG beads in the protein excluding virtual Gō beads (`--Natoms`), the number of missing residues in the atomistic structure (`--missres`), as well as the contact map file (`-f`). 

The prefix of the generated files can be specified (`--moltype`) to e.g. distinguish different protein chains which are not connected by Gō-like bonds. This prefix must match the one assigned in Step 1 on martinize2 with the flag `--govs-moltype`.

In addition, the dissociation energy of the Lennard-Jones potentials (`--go_eps`) as well as the lower (`--cutoff_short`) and upper (`--cutoff_long`) distance cutoffs for two connected backbone beads can also be specified.

The following settings are strongly recommended:
``````
--go_eps 9.414 (kJ/mol)
--cutoff_short 0.3 (nm)
--cutoff_long 1.1 (nm)
``````

A typical `create_goVirt` command for generating a Gō-ready `protein.itp` may look something like this:

````
create_goVirt -s cg.pdb -f contactmap.map --moltype molecule_0 --go_eps 9.414
````

#### 3.2. Output

The script generates four files with parameter details required for the GōMartini model. The file names have a prefix to specify the respective protein or chain (defined by `--moltype`). The name ending as well as the content of each file is listed below: 

| File name | Description |
| ----------- | ----------- |
| `<moltype>_BB-part-def_VirtGoSites.itp` | Bead definitions for the virtual particles. |
| `<moltype>_go-table_VirtGoSites.itp` | Gō-like bonds interaction table. | 
| `<moltype>_exclusions_VirtGoSites.itp` | Exclusions table. | 
| `<moltype>_go4view_harm.itp` | Harmonic Gō-like bonds for visualization only.


Do not change the names of the files as they are included in various other `.itp` files in a way that requires the correct file names. Additionally, all these files, as well as the main Martini force field ``.itp`` should be present in the same directory.


Because the main Martini force field ``.itp`` should not be changed, the script adds the necessary bead definitions to two files which are included in the force field (via include statements) if the variable `GO_VIRT` is defined in the ``.top`` file. This variable is automatically set in the ``.top`` file supplied by martinize2 if the `--govs-include` flag is specified.

### 4. Scaling backbone — Water (BB—W) interactions

`create_goVirt` can also be used to apply a bias to scale protein-water interactions. The virtual interaction sites are typically used to solely encode the LJ potentials between virtual site pairs that make up the Gō-like scaffold. However, non-bonded interactions between these virtual sites and other regular Martini 3 beads — such as the water bead, W — can also be defined. This can be particularly useful to correct some issues pertaining to unstable helical transmembrane behavior or over-aggregation of intrinsically disordered proteins. 

#### 4.1. Bias a sequence of consecutive amino acid residues
A bias can be easily applied to a sequence of consecutive amino acid residues by using the `--idp_start` and `--idp_end` flags to define the first and last amino acid residue number of the sequence to which the bias should be applied and  `--idp_eps` to define the dissociation energy (kJ/mol) of the Lennard-Jones potential used to modulate the BB-W interaction via the virtual backbone beads. Note that `--idp_eps` can be either a positive value (increasing BB-W interactions) or a negative value (decreasing BB-W interactions). 

A typical `create_goVirt` command for generating a Gō-ready `protein.itp` with added BB-W interaction scaling may look something like this:

````
create_goVirt -s cg.pdb -f contactmap.map --moltype molecule_0 --go_eps 9.414 --idp_start 4 --idp_end 30 --idp_eps 0.5
````

Note that you can also use `create_goVirt` to build a model that applies the BB-W interaction scaling, without applying the Gō-like scaffold by not defining `--go_eps` or setting `--go_eps 0`.

#### 4.2. Automatically apply bias depending on secondary structure

`create_goVirt` can automatically apply a particular bias depending on the secondary structure motif (alpha-helix, beta-sheet or random coil) associated with each protein residue. To do this `--bias_auto True` can be set. The dissociation energy (kJ/mol) of the Lennard-Jones potential used to modulate the BB-W interactions associated with each secondary structure motif can be set using `--bias_alfa`, `--bias_beta`, and `--bias_coil` for alpha-helix, beta sheet and random coils, respectively. If set to zero (`0`) no bias is applied to residues assigned to this motif. As the secondary structure information is obtained from the protein `.itp` file produced by martinize2, the `--itp` flag must be applied.

A typical `create_goVirt` command for generating a Gō-ready `protein.itp` with added BB-W interaction scaling may look something like this:

````
create_goVirt -s cg.pdb -f contactmap.map --moltype molecule_0 --go_eps 9.414 --bias_auto True --itp molecule_0.itp --bias_alfa -0.5 --bias_beta 0.0 --bias_coil 0.5
````
