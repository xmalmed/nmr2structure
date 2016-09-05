#!/bin/bash

show_help() {
cat <<EOF
Prepare distance restraints (NOE) for gromacs from cyana (xeasy) input (distance.upl)
The conversion is done by (g)awk script nmr2gmx.awk written by EIZO AB.

Prepare distance restraints (NOE) for gromacs from cyana (xeasy) input (distance.upl)

Usage: $0 -h -f gromacs-topology.{itp,top,gro} -u input.upl [-o disres.itp] [-k 1] [-l 0]

    -f FILE : input gromacs topology file with protein atom numbering {.itp,.top}
    -u FILE : input file name with distance upper limits constraints (.upl)
    -o FILE : output distance restraints file name (disres.itp)
    -k NUM : is factor for multiply forcefield constant of distance restraints
    -l NUM : in angstrom, lengthening of NOE distances, default 0A
    -h help 
EOF
}

k=1
l=0
output="disres.itp"
input=
inputUPL=
OPTIND=1
while getopts "f:u:o:k:l:h?" opt; do
    case $opt in
        f)
            echo "Input $OPTARG" >&2
            input=$OPTARG
            ;;
        u)
            echo "Input $OPTARG" >&2
            inputUPL=$OPTARG
            ;;
        o)
            echo "Output $OPTARG" >&2
            output=$OPTARG
            ;;
        k)
            k=$OPTARG
            ;;
        l)
            l=$OPTARG
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

if [ ! -f "$inputUPL" ]; then
    show_help >&2
    echo "ERROR, dihedral input $inputUPL does not exist." >&2
    exit 1
fi 

if [ -f "$output" ]; then
    mv $output ${output}.BAK
    echo "NOTE, old $output was backuped to ${output}.BAK" >&2
fi 

echo "Scaling constant (factor) is $k" >&2
echo "NOEs extend by $l A" >&2

# script body...
# convert QQD, QQG pseudo-atoms
awk '{print "constraint", $1, $2, $3, $4, $5, $6, $7, 0.0, 1}' $inputUPL |sed '
/QQD/ { s/ QQD / HD2# /
    s/1$/2/p
    s/ HD2#/ HD1#/
    s/2$/1/}
/QQG/ { s/ QQG / HG2# /
    s/1$/2/p
    s/ HG2#/ HG1#/
    s/2$/1/}' | sed '
/QQD/ {s/ QQD / HD2# /
    s/2$/4/p
    s/1$/2/p
    s/\(.*\) HD2#/\1 HD1#/
    s/4$/3/
    s/2$/1/}
/QQG/ {s/ QQG / HG2# /
    s/2$/4/p
    s/1$/2/p
    s/\(.*\) HG2#/\1 HG1#/
    s/4$/3/
    s/2$/1/}
    
s/ QA / HA# /g
s/ QB / HB# /g
s/ QG / HG# /g
s/ QD / HD# /g
s/ QE / HE# /g
s/ QE2 / HE2# /g
s/ QG1 / HG1# /g
s/ QG2 / HG2# /g
s/ QD1 / HD1# /g
s/ QD2 / HD2# /g ' > tmp.upl


# nmr2gmx.awk has many useful parameters
DIR=`dirname $0`
awk -f ${DIR}/nmr2gmx.awk -v weight=$k -v lengthening=$l $input tmp.upl > $output

rm tmp.upl

echo "Distance constraints were stored in $output"



