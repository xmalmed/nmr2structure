#!/bin/bash
# Idea is to run this script in a loop, e.g.:  
# for i in {1..20}; do 
#   ./script.sh $i
# done 
#
# where $i is number of starting structure from cyana.
#
# NOTE: due to a bug in distance restraints, only save option is to run MD on 1 CPU or 1 GPU... 

show_help() {
cat <<EOF
Quickly warm the system in a high temperature.

Usage: `basename $0` -n number [-p parameters.txt] [-h]
or:    `basename $0` number

    -p parameters file name, default 'parameters.txt'
    -n number of run and/or starting structure
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
    echo "ERROR, number of run was not given." >&2
    exit 1
fi 

echo "Run number $number "
echo "Reading parameters file $input " 

# Reading parameters parameters:
#------------------------------------------------------------------------------

name="`awk '$1 == "name" {print $2}' $input `"

dirName="`awk '$1 == "warming_directory" {print $2}' $input `" 
prevDirName="`awk '$1 == "equilibration_directory" {print $2}' $input `" 

dihres="`awk '$1 == "dihedral_restraints" {print $2}' $input `"
disres="`awk '$1 == "distance_restraints" {print $2}' $input `"

temp="`awk '$1 == "normal_temperature" {print $2}' $input `"
tempHigh="`awk '$1 == "high_temperature" {print $2}' $input `"
#-----------------------------------------------------------

mkdir ${dirName}-$number
#topology and restraints
cp ${prevDirName}-$number/*.top ${dirName}-$number/   # topop.top
cp ${prevDirName}-$number/*.itp ${dirName}-$number/   # dihres.itp, disres.itp, posres.itp, ....
#structure
cp ${prevDirName}-$number/npt.gro  ${prevDirName}-$number/npt.cpt ${dirName}-$number/
cd ${dirName}-$number/


#insert DIHRES restraints in the topology after #ifdef POSRES
for top in *.{top,itp}; do
    grep -q ' *#ifdef DIHRES *' $top || sed -i "/#ifdef POSRES$/ {N;N;N; a; Include Dihedral restraint file\n#ifdef DIHRES\n  #include \"`basename $dihres`\"\n#endif\n\n
    }" $top && grep -q ' *#ifdef DIHRES *' $top && echo "DIHRES Restraints were add in topology $top"
done

#insert DISRES restraints in the topology after #ifdef POSRES
for top in *.{top,itp}; do
    grep -q ' *#ifdef DISRES *' $top || sed -i "/#ifdef POSRES$/ {N;N;N; a; Include Distance restraint file\n#ifdef DISRES\n  #include \"`basename $disres`\"\n#endif\n\n
    }" $top && grep -q ' *#ifdef DISRES *' $top && echo "DISRES Restraints were add in topology $top"
done

echo ""

#warming
cat <<< "
; Preprocessing
define        = -DDISRES -DDIHRES   ; distance and dihedral restraints

; Run control
integrator    = md            ; leap-frog integrator
nsteps        = 25000         ; 50 ps
dt            = 0.002         ; 2 fs

; Output control
nstxout        = 1000          ; save coordinates every 2.0 ps
nstvout        = 1000          ; save velocities every  2.0 ps
nstenergy      = 1000          ; save energies every    2.0 ps
nstlog         = 1000          ; update log file every  2.0 ps
;energygrps    = Protein Non-Protein

; Neighbor searching
cutoff-scheme      = Verlet
ns_type            = grid     ; search neighboring grid cells
nstlist            = 20       ; 20 fs, largely irrelevant with Verlet
pbc                = xyz      ; Periodic boundary conditions, 3-D PBC 

; Electrostatics
coulombtype     = PME         ; Particle Mesh Ewald for long-range electrostatics
rcoulomb        = 1.0         ; short-range electrostatic cutoff (in nm)

; VdW
rvdw            = 1.0         ; short-range van der Waals cutoff (in nm)
DispCorr        = EnerPres    ; Dispersion correction, account for cut-off vdW scheme

; Ewald
fourierspacing    = 0.16      ; grid spacing for FFT
pme_order         = 4         ; cubic interpolation

; Temperature coupling ON
tcoupl        = V-rescale              ; modified Berendsen thermostat
tc-grps       = Protein Non-Protein    ; two coupling groups - more accurate (separate baths for Protein, Non-Protein)
tau_t         = 0.1       0.1          ; Coupling time constant, one for each group, in ps
ref_t         = $temp     $temp        ; reference temperature, one for each group, in K

; Pressure coupling ON
pcoupl              = no             ; Pressure coupling OFF

; Simulated annealing
annealing           = single single     ; for each tc-grps 
annealing-npoints   = 2 2               ; # of points (temp) for each group ... 2 = starting, end temp.
annealing-time      = 0 50 0 50         ; ps, times when the given temp should be reach
annealing-temp      = $temp $tempHigh $temp $tempHigh ; K, temps for each ann.-npoints

; Velocity generation
gen_vel     = no         ; assign velocities from Maxwell distribution; 

; Bonds
constraints             = all-bonds   ; all bonds (even heavy atom-H bonds) constrained
constraint_algorithm    = lincs       ; holonomic constraints 
continuation            = yes         ; SECOND dynamics run; ADD -t nvt.cpt!
lincs_order             = 4           ; also related to accuracy
lincs_iter              = 1           ; accuracy of LINCS

; NMR refinement
disre           = simple        ; distance restraints on
disre-weighting = equal
disre-mixed     = no
disre-fc        = 1000          ; kJ / mol / nm2 * factor (weight)
disre-tau       = 0             ; ps
nstdisreout     = 100          ; steps, write to energy file
orire           = no            ; orientation restraints
" > warm.mdp

#mvt.tpr   ... continue with nvt.gro, -t nvt.cpt (velocity etc)
gmx grompp -f warm.mdp -c npt.gro -t npt.cpt -p topol.top -o warm.tpr
if [ $? -ne 0 ]; then exit; fi

#MD npt equil.
gmx mdrun -deffnm warm -ntmpi 1 -ntomp 1
if [ $? -ne 0 ]; then exit; fi


#    #temperature
#    gmx energy -f warm.edr -o temperature.xvg  << EOF
#    17
#    0
#    EOF
#
#    #presure
#    gmx energy -f warm.edr -o pressure.xvg  << EOF
#    19
#    0
#    EOF
#
#    #distres
#    gmx energy -f warm.edr -o disres.xvg << EOF
#    10
#    0
#    EOF
#
#    #distres viol
#    gmx energy -f warm.edr -o disresViol.xvg << EOF
#    11
#    0
#    EOF
#
#    #dihres
#    gmx energy -f warm.edr -o dihres.xvg << EOF
#    12
#    0
#    EOF

#check.pdb
gmx trjconv -f warm.gro -s warm.tpr -o check.pdb -pbc mol -ur compact <<< 0

cd ..



