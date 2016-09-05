#!/bin/bash
# This script runs just once for all prepared files
#
# NOTE: due to a bug in distance restraints, only save option is to run MD on 1 CPU or 1 GPU... 
# this script prepares file to be run in metacentrum with infinity on GPUs:


show_help() {
cat <<EOF
Prepare tar archive (cyana format) to upload in iCING https://nmr.le.ac.uk/icing with all final structures (models in one .pdb file) for analyses.

Usage: `basename $0` [-p parameters.txt] [-h] 
    
    -p parameters file name, default 'parameters.txt'
    -h help
EOF
}

trajExport=0
number=
input="parameters.txt"
OPTIND=1
while getopts "p:h?" opt; do
    case $opt in
        p)
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

echo "Reading parameters file $input " 

# Reading parameters parameters:
#------------------------------------------------------------------------------

name="`awk '$1 == "name" {print $2}' $input `"

dirName="`awk '$1 == "cooling_directory" {print $2}' $input `" 

dihres="`awk '$1 == "dihedral_restraints_cyana" {print $2}' $input `"
disres="`awk '$1 == "distance_restraints_cyana" {print $2}' $input `"
seq="`awk '$1 == "sequence" {print $2}' $input `"
#-----------------------------------------------------------
dirICING="${dirName}-iCING"
mkdir ${dirICING}

if [ -s pdb.tmp ]; then
    rm join.temp
fi

if [ -s remark.tmp ]; then
    rm remark.tmp
fi

num=0
for i in ${dirName}-[0-9]* ; do
    k=`echo $i | tr -d '/' | grep -oE '[0-9]+$'`  # number of run
 
    if [ -s ${i}/${name}.pdb ]; then
        echo "Concatenate ${i}/${name}.pdb"
        cat ${i}/${name}.pdb >> pdb.tmp
        echo "REMARK structure ${i}/${name}.pdb " >> remark.tmp
        grep -E "REMARK|CRYST1|TITLE" ${i}/${name}.pdb >> remark.tmp
        echo "REMARK" >> remark.tmp
        if [ "`tail -1 ${i}/${name}.pdb | cut -c1-3`" != "END" ] ; then
            echo "END" >> pdb.tmp
        fi
        num=$((num+1))    
    else
        echo "Skipping directory $i "
    fi
done

awk '
BEGIN {
    ModelNumber = 1
    NewModel    = "yes"
}

# first ATOM card signals start of model
/^ATOM/ {
    if ( NewModel == "yes" ) {
	printf( "MODEL     %4i\n", ModelNumber )
	NewModel = "no"
    }
}

/^CONECT/ { next }

/^MASTER/ { next }

# END card signals end of model
/^END/ {
    ModelNumber++
    NewModel = "yes"
    print "ENDMDL"
    next
}

# print only models...
NewModel == "no" { print $0 }

END { print "END" }
' pdb.tmp | cat remark.tmp - > ${name}-${num}mdl.pdb

rm pdb.tmp
rm remark.tmp


#ions are align wrong by 1 space
sed -i "s/ \([A-Z][A-Z]\)   \1 /\1   \1  /" ${name}-${num}mdl.pdb

cp ${name}-${num}mdl.pdb ${dirICING}/${name}.pdb

#rename to cyana, 1->3 in -CH2- groups:
./rename_gromacs2cyana.sh ${dirICING}/${name}.pdb

cp $dihres ${dirICING}/${name}.aco
cp $disres ${dirICING}/${name}.upl
cp $seq ${dirICING}/${name}.seq

tar -czf ${dirICING}.cyana.tgz ${dirICING}
echo "Archive ${dirICING}.cyana.tgz for iCING was prepared."






