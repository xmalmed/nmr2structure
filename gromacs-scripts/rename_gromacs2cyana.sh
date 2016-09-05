#!/bin/bash
# rename "in-place" pdb from gromacs to Cyana 
if [ "$#" -eq 0 ]; then
    echo "input structure .pdb"
    exit
fi

if [ "$#" -gt 1 ]; then
    echo "Too much args."
    exit
fi

#PDB edit Hx1 -> Hx3 for -CH2- groups
# reneme back ILE CD1, [123]HD1
sed -i '
s/1HG1 ILE/3HG1 ILE/
s/CD  ILE/CD1 ILE/
s/ HD1 ILE/1HD1 ILE/
s/ HD2 ILE/2HD1 ILE/
s/ HD3 ILE/3HD1 ILE/
s/H1  PRO/H3  PRO/
s/HB1 PRO/HB3 PRO/
s/HG1 PRO/HG3 PRO/
s/HD1 PRO/HD3 PRO/
s/HB1 LEU/HB3 LEU/
s/HB1 ASP/HB3 ASP/
s/HB1 ASN/HB3 ASN/
s/HA1 GLY/HA3 GLY/
s/HB1 GLU/HB3 GLU/
s/HG1 GLU/HG3 GLU/
s/HB1 GLN/HB3 GLN/
s/HG1 GLN/HG3 GLN/
s/HB1 LYS/HB3 LYS/
s/HG1 LYS/HG3 LYS/
s/HD1 LYS/HD3 LYS/
s/HE1 LYS/HE3 LYS/
s/HB1 SER/HB3 SER/
s/HB1 ARG/HB3 ARG/
s/HG1 ARG/HG3 ARG/
s/HD1 ARG/HD3 ARG/
s/HB1 CYS/HB3 CYS/
s/HB1 MET/HB3 MET/
s/HG1 MET/HG3 MET/
s/HB1 PHE/HB3 PHE/
s/HB1 HIS/HB3 HIS/
s/HB1 TRP/HB3 TRP/
s/HB1 TYR/HB3 TYR/
' $1


