#!/bin/bash
# Idea is to run this script in a loop, e.g.:  
# for i in {1..20}; do 
#   ./topology.sh $i
# done 
#
# where $i is number of starting structure from cyana.
#
# NOTE: due to a bug in distance restraints, only save option is to run MD on 1 CPU or 1 GPU... 


show_help() {
cat <<EOF
Prepare topology with solvent in periodic box.

Usage: `basename $0` -n number [-p parameters.txt] [-h]
or:    `basename $0` number

    -p parameters file name, default 'parameters.txt'
    -n number of the run and starting structure index
    -h help
EOF
}

number=
input="parameters.txt"
OPTIND=1
while getopts "n:p:h?" opt; do
    case $opt in
        p)
            input=$OPTARG
            ;;
        n)
            number=$OPTARG
            echo "Run number $number " >&2
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

if [ "$#" -eq 1 ]; then
    number=$1
fi

if [ -z "$number" ]; then
    show_help >&2
    echo "ERROR, number of the run was not given." >&2
    exit 1
fi 

echo "Run number $number "
echo "Reading parameters file $input " 
echo ""

# Reading parameters parameters:
#------------------------------------------------------------------------------

name="`awk '$1 == "name" {print $2}' $input `"

dirName="`awk '$1 == "topology_directory" {print $2}' $input `" 
pdbName=$(eval echo "`awk '$1 == "starting_structures" {print $2}' $input `")     # takes pdb name with ${number} variable from parameters.txt and eval ${number} to -n number.
saltConcentration="`awk '$1 == "salt_concentration" {print $2}' $input `"   # mol/liter 

#------------------------------------------------------------------------------
# The most common setting for protein in water was used below. You might want to change some. 

mkdir ${dirName}-$number
cp $pdbName ${dirName}-$number/
cd ${dirName}-$number


#convert pdb to gmx
# -ignh ... not used to control names of protons, due to NOE distances.
gmx pdb2gmx -ff amber99sb-ildn -f `basename $pdbName` -o ${name}.gro -water tip3p -merge all
if [ $? -ne 0 ]; then exit; fi

#periodic boundary conditions
#cube... -bt (box type) cubic, dodecahedron
#center protein ... -c 
#distance from the edge, -d 1.0 [nm]
#d=0.6 # nm; space around protein in periodic box, 1.0 is safe
gmx editconf -f ${name}.gro -o ${name}-box.gro -c -d 0.6 -bt dodecahedron
if [ $? -ne 0 ]; then exit; fi

#solvate
# -cp ... configuration of the protein
# -cs ... configuration of the solvent, spc216.gro for 3 point model (tip3p)
gmx solvate -cp ${name}-box.gro -cs spc216.gro -o ${name}-box-sol.gro -p topol.top
if [ $? -ne 0 ]; then exit; fi


# energy minimization

cat <<< "
;Preprocessing
;define = -DFLEXIBLE ;flexible waters
;define = -DPOSRES ; include posre.itp, if it is define in topology.

;Run control
integrator	= steep		; Algorithm (steep = steepest descent minimization)
nsteps		= 5000	  	; Maximum number of (minimization) steps to perform

;Energy minimization 
emtol		= 500.0  	; Stop minimization when the maximum force < 500.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size

;Output control
; nstxout     = 50        ; number of steps between writing coordinates

;Neighbor searching
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 20		; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
pbc		        = xyz 		; Periodic Boundary Conditions (xyz/no/xy)

;Electrostatics
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; Short-range electrostatic cut-off [nm]

;VdW
rvdw		    = 1.0		; Short-range Van der Waals cut-off [nm]
" > em.mdp 

#ions.tpr
gmx grompp -f em.mdp -c ${name}-box-sol.gro -p topol.top -o ions.tpr
if [ $? -ne 0 ]; then exit; fi

#ions adding
gmx make_ndx -f ${name}-box-sol.gro -o tmp.ndx <<< q | grep " SOL " > tmp
sol=`awk '$2=="SOL" {print $1}' tmp`
rm tmp tmp.ndx

gmx genion -s ions.tpr -p topol.top -o ${name}-box-sol-ions.gro -conc $saltConcentration -neutral <<< $sol 
if [ $? -ne 0 ]; then exit; fi

# EM
gmx grompp -v -f em.mdp -c ${name}-box-sol-ions.gro -p topol.top -o ${name}-box-sol-ions-em.tpr
if [ $? -ne 0 ]; then exit; fi
# run EM
gmx mdrun -deffnm ${name}-box-sol-ions-em -ntmpi 1 -ntomp 1
if [ $? -ne 0 ]; then exit; fi

#check2.gro  
gmx trjconv -f ${name}-box-sol-ions-em.gro -s ${name}-box-sol-ions-em.tpr -o check.pdb -pbc mol -ur compact <<< 0
if [ $? -ne 0 ]; then exit; fi

gmx trjconv -f ${name}-box-sol-ions-em.gro -s ${name}-box-sol-ions-em.tpr -o start.gro -pbc mol -ur compact <<< 0

echo ""
echo "Check minimized structure, check.pdb!"
echo ""

cd ..

