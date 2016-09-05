#!/bin/bash

# help
show_help() {
cat <<EOF

Convert IUPAC nomenclature of atoms for CYANA.
Prepares file for TALOS (-t)

Usage: $0 [-ht] -f file.shifts -s file.seq [-o output.prot] 
-s sequence file:
    GLY -2
    SER -1
    ASP 0
    ALA 1
-f chemical shifts file in XEASY format with IUPAC nomenclature:
   1  177.681  0.000      C 1
   2   52.266  0.011     CA 1
   3   23.209  0.022     CB 1
   4    8.604  0.003      H 1
   5    5.082  0.002     HA 1
   6    1.260  0.002     MB 1
   7  119.059  0.027      N 1
-o output file (output.prot)
-t Creates output file (output.tab) which can be used by TALOS (NMRPipe) for prediction of dihedral angles.
-h help            

EOF
}

output="output.prot"
talos=0
input=
tab=
sequence=
OPTIND=1
while getopts "s:f:o:th?" opt; do
    case $opt in
        s)
            echo "Sequence $OPTARG" >&2
            sequence=$OPTARG
            ;;
        f)
            echo "Chemical shifts $OPTARG" >&2
            input=$OPTARG
            ;;
        o)
            echo "Output $OPTARG" >&2
            output=$OPTARG
            ;;
        t)
            talos=1
            ;;
        h|\?)
            show_help >&2
            exit 0
            ;;
    esac
done

# checks:
if [ ! -f "$sequence" ]; then
    echo "ERROR, Sequence does not exist." >&2
    show_help >&2
    exit 1
fi 

if [ ! -f "$input" ]; then
    echo "ERROR, Chemical shift file does not exist." >&2
    show_help >&2
    exit 1
fi 


# TALOS file:
if [ $talos ] ; then
    tab="`echo $output | grep -o '.*\.'`tab"
    echo "Creating file with resonances for TALOS:" >&2
    ./xeasy-resonances2talos.sh -f "$input" -s "$sequence" -o "$tab"
    echo "Outputting file for TALOS is $tab"
fi 


# some parameters:
firstRES=`head -n 1 $sequence | awk '{print $2}'`
lines=`cat $sequence |wc -l`
end=$((lines + firstRES))


# Conversion from IUPAC to CYANA:
cat $input | awk -v first=$firstRES -v end=$end -v seq=$sequence -v output=$output '
#read sequence file and store in array a[]
BEGIN {
    numShift = 0
    for (i=first; i<end; i++) {getline < seq; a[i] = $1}
}

#substitutions from IUPAC XEASY format to cyana
{
    changed = 0
    $1 += numShift
    comment = "# "$4
    if (a[$5] == "ALA") {
        if ("MB" == $4) {$4 = "QB"; changed = 1}
        
    }
    
    if (a[$5] == "VAL") {
        if ("MG1" == $4) {$4 = "QG1"; changed = 1}
        if ("MG2" == $4) {$4 = "QG2"; changed = 1}
        if ("QMG" == $4) {$4 = "QQG"; changed = 1}
    }
    
    if (a[$5] == "THR") { 
        if ($4 == "MG2") {$4 = "QG2"; changed = 1}
    } 
    
    if (a[$5] == "MET") { 
        if ($4 == "ME") {$4 = "QE"; changed = 1}
    } 
          
    if (a[$5] == "LEU") {
        if ("MD1" == $4) {$4 = "QD1"; changed = 1}
        if ("MD2" == $4) {$4 = "QD2"; changed = 1}
        if ("QMD" == $4) {$4 = "QQD"; changed = 1}
    }
          
    if (a[$5] == "ILE") {
        if ("MG2" == $4) {$4 = "QG2"; changed = 1}
        if ("MD1" == $4) {$4 = "QD1"; changed = 1}
        if ("QMD" == $4) {$4 = "QQD"; changed = 1}
    } 
    
    if ("CQD" == $4) {
        changed = 1
        numShift ++
        $4 = "CD1"
        printf "%4i %8.3f %6.3f %6s %-4i %s\n",$1,$2,$3,$4,$5,comment > output
        $1 ++
        $4 = "CD2"
    }
    
    if ("CQE" == $4) {
        changed = 1
        numShift ++
        $4 = "CE1"
        printf "%4i %8.3f %6.3f %6s %-4i %s\n",$1,$2,$3,$4,$5,comment > output
        $1 ++
        $4 = "CE2"
    }
    
    if ( ! changed ) {comment = ""}
            
    printf "%4i %8.3f %6.3f %6s %-4i %s\n",$1,$2,$3,$4,$5,comment > output
}'
     

