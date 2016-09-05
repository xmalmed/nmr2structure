#!/bin/bash

show_help() {
cat <<EOF
Transform sparky resonance list 'sparky shortcut rl' to NMRPipe format for TALOS(+/N).

Usage: $0 [-h] -f input.shifts [-o output.tab]

    -f input file name
    -o output file name (output.tab)
    -h help
EOF
}

output="output.tab"
input=
OPTIND=1
while getopts "f:o:h?" opt; do
    case $opt in
        f)
            echo "Input $OPTARG" >&2
            input=$OPTARG
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

if [ -f "$output" ]; then
    mv $output ${output}.BAK
    echo "NOTE, old $output was backuped to ${output}.BAK" >&2
fi 

# script body...
# header
echo "DATA FIRST_RESID 1" > $output
echo "DATA SEQUENCE ..." >> $output
echo "" >> $output 
echo "VARS   RESID RESNAME ATOMNAME SHIFT" >> $output
echo 'FORMAT %4d   %1s     %4s      %8.3f' >> $output
echo "" >> $output 

# resonances
cat $input | sed 's/ H /HN /
/ QA / {s/ QA /HA2 /p ; s/HA2 /HA3 / }' | \
awk '/[A-Z][0-9]+/ {if($2=="C" || $2=="CA" || $2=="CB" || $2=="N" || $2=="HN" || $2=="HA" || $2=="HA2" || $2=="HA3" )
    {printf "%4d %1s %4s %8.3f\n", substr($1,2,length($1)), substr($1,1,1), $2, $4}}' >> $output

echo "WARN, output file $output needs adjust header!"


