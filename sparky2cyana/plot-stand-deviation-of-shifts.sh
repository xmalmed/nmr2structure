#!/bin/bash

show_help() {
cat <<EOF
Create histogram of standart deviations for atoms 1H, 13C, and 15N. This script use OCTAVE.
Output are [H,C,N]-std-hist.svg

Usage: script -ahf input.file -o output.file 
    -f input file name of XEASY format of chemical shifts with 3. row of stand. dev.
    -h help
EOF
}

input=
OPTIND=1
while getopts "f:h?" opt; do
    case $opt in
        f)
            echo "Input $OPTARG" >&2
            input=$OPTARG
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
    echo "ERROR, $input does not exist." >&2
    exit 1
fi 

# script body...
#check octave
type octave >/dev/null 2>&1 || { echo >&2 "Octave is not installed.  Aborting."; exit 1; }


grep " [HM]" $input | awk '{print $3}' > tmpH1.txt
grep  " [C]" $input | awk '{print $3}' > tmpC1.txt
grep  " [N]" $input | awk '{print $3}' > tmpN1.txt
#mean and std deviation
#mean=`awk '{sum += $3} END {print sum / NR}' tmpH.txt`   
#std=`awk -v mean=$mean '{sumsq += ($3 - mean)^2} END {print sqrt( sumsq / NR )}' tmpH.txt`

echo "Creating histograms [H,C,N]-std-hist.svg"

octave << EOF > /dev/null 2>&1 
Hdata = load ('tmpH1.txt');
hist(Hdata,30);
print("-dsvg", "H-std-hist.svg");
EOF


octave << EOF > /dev/null 2>&1 
Cdata = load ('tmpC1.txt');
hist(Cdata,30);
print("-dsvg", "C-std-hist.svg");
EOF

octave << EOF > /dev/null 2>&1 
Ndata = load ('tmpN1.txt');
hist(Ndata,30);
print("-dsvg", "N-std-hist.svg");
EOF

rm tmpH1.txt tmpC1.txt tmpN1.txt
