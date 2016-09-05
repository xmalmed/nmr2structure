#!/bin/bash

show_help() {
cat <<EOF
This script rename mostly -CH2- groups from H2 and H3 to H2 and H1, (3 -> 1), e.g.:
SER HB2 -> SER HB2
SER HB3 -> SER HB1
You should rename .pdb file and noe distances .upl before converting these files to gromacs.

Usage: $0 -h [-s structures.pdb] [-o output.pdb] [-u neo.upl] [-p output-noe.upl]

    -f input .pdb file
    -o output of renamed .pdb (output.pdb)
    -u input .upl file
    -p output of renamed .upl (output.upl)
    -h help
EOF
}

output="output.pdb"
outputUPL="output.upl"
input=
inputUPL=
flagS=0
flagU=0
OPTIND=1
while getopts "s:o:u:p:h?" opt; do
    case $opt in
        s)
            echo "Structures $OPTARG" >&2
            input=$OPTARG
            flagS=1
            ;;
        o)
            output=$OPTARG
            ;;
        u)
            echo "NOE distances $OPTARG" >&2
            inputUPL=$OPTARG
            flagU=1
            ;;
        p)
            outputUPL=$OPTARG
            ;;
        h|\?)
            show_help >&2
            exit 0
            ;;
    esac
done

shift $((OPTIND-1))
[ "$1" = "--" ] && shift
# next args can be processed, $1, $2...

# script body...

if [ $flagS == 1 ]; then
    # RENAME of structures in pdb Hx3 -> Hx1 for -CH2- groups 
    # gmx pdb2gmx rename ILE CD1, HD1[123] to CD, HD[123]
    if [ ! -f "$input" ]; then
        show_help >&2
        echo "ERROR, input $input does not exist." >&2
        exit 1
    fi 

    if [ -f "$output" ]; then
        mv $output ${output}.BAK
        echo "NOTE, old $output was backuped to ${output}.BAK" >&2
    fi 
    
    sed -e '
    /^ATOM/ s/HG13 ILE/HG11 ILE/
    /^ATOM/ s/H3  PRO/H1  PRO/
    /^ATOM/ s/HB3 PRO/HB1 PRO/
    /^ATOM/ s/HG3 PRO/HG1 PRO/
    /^ATOM/ s/HD3 PRO/HD1 PRO/
    /^ATOM/ s/HB3 LEU/HB1 LEU/
    /^ATOM/ s/HB3 ASP/HB1 ASP/
    /^ATOM/ s/HB3 ASN/HB1 ASN/
    /^ATOM/ s/HA3 GLY/HA1 GLY/
    /^ATOM/ s/HB3 GLU/HB1 GLU/
    /^ATOM/ s/HG3 GLU/HG1 GLU/
    /^ATOM/ s/HB3 GLN/HB1 GLN/
    /^ATOM/ s/HG3 GLN/HG1 GLN/
    /^ATOM/ s/HB3 LYS/HB1 LYS/
    /^ATOM/ s/HG3 LYS/HG1 LYS/
    /^ATOM/ s/HD3 LYS/HD1 LYS/
    /^ATOM/ s/HE3 LYS/HE1 LYS/
    /^ATOM/ s/HB3 SER/HB1 SER/
    /^ATOM/ s/HB3 ARG/HB1 ARG/
    /^ATOM/ s/HG3 ARG/HG1 ARG/
    /^ATOM/ s/HD3 ARG/HD1 ARG/
    /^ATOM/ s/HB3 CYS/HB1 CYS/
    /^ATOM/ s/HB3 MET/HB1 MET/
    /^ATOM/ s/HG3 MET/HG1 MET/
    /^ATOM/ s/HB3 PHE/HB1 PHE/
    /^ATOM/ s/HB3 HIS/HB1 HIS/
    /^ATOM/ s/HB3 TRP/HB1 TRP/
    /^ATOM/ s/HB3 TYR/HB1 TYR/
    ' $input > $output
    echo "$input was renamed and stored in $output"
fi

if [ $flagU == 1 ]; then
    # UPL edit with a consistency to PDB...
    if [ ! -f "$inputUPL" ]; then
        show_help >&2
        echo "ERROR, input $inputUPL does not exist." >&2
        exit 1
    fi 

    if [ -f "$outputUPL" ]; then
        mv $outputUPL ${outputUPL}.BAK
        echo "NOTE, old $outputUPL was backuped to ${outputUPL}.BAK" >&2
    fi 
    
    sed -e '
    s/ILE  HG13/ILE  HG11/g
    s/PRO  H3 /PRO  H1 /g
    s/PRO  HB3/PRO  HB1/g
    s/PRO  HG3/PRO  HG1/g
    s/PRO  HD3/PRO  HD1/g
    s/LEU  HB3/LEU  HB1/g
    s/ASP  HB3/ASP  HB1/g
    s/ASN  HB3/ASN  HB1/g
    s/GLY  HA3/GLY  HA1/g
    s/GLU  HB3/GLU  HB1/g
    s/GLU  HG3/GLU  HG1/g
    s/GLN  HB3/GLN  HB1/g
    s/GLN  HG3/GLN  HG1/g
    s/LYS  HB3/LYS  HB1/g
    s/LYS  HG3/LYS  HG1/g
    s/LYS  HD3/LYS  HD1/g
    s/LYS  HE3/LYS  HE1/g
    s/SER  HB3/SER  HB1/g
    s/ARG  HB3/ARG  HB1/g
    s/ARG  HG3/ARG  HG1/g
    s/ARG  HD3/ARG  HD1/g
    s/CYS  HB3/CYS  HB1/g
    s/MET  HB3/MET  HB1/g
    s/MET  HG3/MET  HG1/g
    s/PHE  HB3/PHE  HB1/g
    s/HIS  HB3/HIS  HB1/g
    s/TRP  HB3/TRP  HB1/g
    s/TYR  HB3/TYR  HB1/g
    ' $inputUPL > $outputUPL
    echo "$inputUPL was renamed and stored in $outputUPL"
fi    
