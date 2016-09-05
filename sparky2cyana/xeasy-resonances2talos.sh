#!/bin/bash

show_help() {
cat <<EOF
Transform XEASY resonance list to NMRPipe format for TALOS(+/N).
Not well tested...

Usage: $0 [-h] -f input.shifts -s sequence.seq [-o output.tab]

    -f input file name
    -s sequence file:
        GLY -2
        SER -1
        ASP 0
        ALA 1
    -o output file name (output.tab)
    -h help
EOF
}

output="output.tab"
input=
sequence=
OPTIND=1
while getopts "s:f:o:h?" opt; do
    case $opt in
        f)
            echo "Input $OPTARG" >&2
            input=$OPTARG
            ;;
        s)
            echo "Sequence $OPTARG" >&2
            sequence=$OPTARG
            ;;
        o)
            echo "Output $OPTARG" >&2
            output=$OPTARG
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
    echo "ERROR, input $input does not exist." >&2
    show_help >&2
    exit 1
fi 

if [ ! -f "$sequence" ]; then
    echo "ERROR, Sequence does not exist." >&2
    show_help >&2
    exit 1
fi 

if [ -f "$output" ]; then
    mv $output ${output}.BAK
    echo "NOTE, old $output was backuped to ${output}.BAK" >&2
fi 

# script body...

firstRES=`head -n 1 $sequence | awk '{print $2}'`
echo "First residue: $firstRES"
lines=`cat $sequence |wc -l`
end=$((lines + firstRES))

# header:
echo "DATA FIRST_RESID $firstRES" > $output


cat $input | sed 's/ H /HN /
/ QA / {s/ QA /HA2 /p ; s/HA2 /HA3 / }' | \
awk -v first=$firstRES -v end=$end -v seq=$sequence -v output=$output '
#read sequence file and store in array a[]
#continue with header
BEGIN {
sequence = "" 
for (i=first; i<end; i++) {
    getline < seq; 
    if       ($1=="ALA")  $1="A"
    else if  ($1~/^CYS/)  $1="C"
    else if  ($1=="ASP")  $1="D"
    else if  ($1=="GLU")  $1="E"
    else if  ($1=="PHE")  $1="F"
    else if  ($1=="GLY")  $1="G"
    else if  ($1=="HIS")  $1="H"
    else if  ($1=="ILE")  $1="I"
    else if  ($1=="LYS")  $1="K"
    else if  ($1=="LEU")  $1="L"
    else if  ($1=="MET")  $1="M"
    else if  ($1=="ASN")  $1="N"
    else if  ($1=="PRO")  $1="P"
    else if  ($1=="GLN")  $1="Q"
    else if  ($1=="ARG")  $1="R"
    else if  ($1=="SER")  $1="S"
    else if  ($1=="THR")  $1="T"
    else if  ($1=="VAL")  $1="V"
    else if  ($1=="TRP")  $1="W"
    else if  ($1=="TYR")  $1="Y"
    else if  ($1=="LYS+") $1="K"
    else if  ($1=="ARG+") $1="R"
    else if  ($1=="HIS+") $1="H"
    else if  ($1=="ASP-") $1="D"
    else if  ($1=="GLU-") $1="E"

    a[i] = $1
    sequence = "" sequence $1
    }
    
printf "DATA SEQUENCE %s\n", sequence >> output
print "" >> output 
print "VARS   RESID RESNAME ATOMNAME SHIFT" >> output
print "FORMAT %4d   %1s     %4s      %8.3f" >> output
print "" >> output 
system("sleep 1") # delay to write header into the file
}

# resonances
{if($4=="C" || $4=="CA" || $4=="CB" || $4=="N" || $4=="HN" || $4=="HA" || $4=="HA2" || $4=="HA3" )
    {printf "%5d %s %4s %8.3f\n", $5, a[$5], $4, $2}}' >> $output


