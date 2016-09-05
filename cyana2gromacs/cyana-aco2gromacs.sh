#!/bin/bash

show_help() {
cat <<EOF
Prepare dihedral angle restraints for gromacs from cyana (xeasy) input (dihedrals.aco)

Usage: $0 -h -f gromacs-topology.{itp,top,gro} -a input.aco [-o dihres.itp] [-k 350]

    -f input gromacs topology file with protein atom numbering {.itp,.top}
    -a input file name with dihedral angle constraints (.aco)
    -o output file name (dihres.itp)
    -k is constant k_dihr for forcefield, kJ mol^(-1) rad^(-2) 
    -h help
    
EOF
}

k=350
output="dihres.itp"
input=
inputACO=
OPTIND=1
while getopts "f:a:o:k:h?" opt; do
    case $opt in
        f)
            echo "Input $OPTARG" >&2
            input=$OPTARG
            ;;
        a)
            echo "Input $OPTARG" >&2
            inputACO=$OPTARG
            ;;
        o)
            output=$OPTARG
            ;;
        k)
            k=$OPTARG
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

if [ ! -f "$input" ]; then
    show_help >&2
    echo "ERROR, topology input $input does not exist." >&2
    exit 1
fi 

if [ ! -f "$inputACO" ]; then
    show_help >&2
    echo "ERROR, dihedral input $inputACO does not exist." >&2
    exit 1
fi 

if [ -f "$output" ]; then
    mv $output ${output}.BAK
    echo "NOTE, old $output was backuped to ${output}.BAK" >&2
fi 

#echo "Output $OPTARG" >&2
echo "Force constant is $k" >&2

# script body...

echo "[ dihedral_restraints ] " > $output
echo "; phi: C'(n-1) - N - CA - C' " >> $output
echo "; psi: N - CA - C' - N(n+1) " >> $output
echo "; ai   aj   ak   al   type    phi    dphi   k_dihr " >> $output

lines=`cat $input | wc -l`

# a[$3 = residue number][$5 = atom name N,CA,C] = $1 ...atom number
awk -v top="$input" -v lines=$lines -v k=$k 'BEGIN {
    for (i=1; i<lines; i++) {
        getline < top
        if ($1 != ";" && NF >= 6 ) {
            # .itp
            if ($5 == "N")  { a[$3][$5] = $1 }
            if ($5 == "CA") { a[$3][$5] = $1 }
            if ($5 == "C")  { a[$3][$5] = $1 }
            # .gro
            if ($2 == "N")  { a[int($1)][$2] = $3 }
            if ($2 == "CA") { a[int($1)][$2] = $3 }
            if ($2 == "C")  { a[int($1)][$2] = $3 }
            }
        }
    }
    $3 == "PHI" {
        printf "%4d %4d %4d %4d %4d %8.2f %8.2f %8.3f\n", a[$1 - 1]["C"], a[$1]["N"], a[$1]["CA"], a[$1]["C"], 1, ($5-$4)/2 + $4, ($5-$4)/2, k
        }
    $3 == "PSI" {
        printf "%4d %4d %4d %4d %4d %8.2f %8.2f %8.3f\n", a[$1]["N"], a[$1]["CA"], a[$1]["C"], a[$1 + 1]["N"], "1", ($5-$4)/2 + $4, ($5-$4)/2, k
        }
' $inputACO >> $output
        
echo "Dihedral constraints were stored in $output" >&2



