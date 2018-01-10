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
Continue with free MD from equlibrate state.

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

dirName="`awk '$1 == "stability_directory" {print $2}' $input `" 
prevDirName="`awk '$1 == "equilibration_directory" {print $2}' $input `" 
firstStructureDir="`awk '$1 == "topology_directory" {print $2}' $input `"

dihres="`awk '$1 == "dihedral_restraints" {print $2}' $input `"
disres="`awk '$1 == "distance_restraints" {print $2}' $input `"

temp="`awk '$1 == "normal_temperature" {print $2}' $input `"
tempStart="`awk '$1 == "normal_temperature" {print $2}' $input `"
#-----------------------------------------------------------
j=50 # ps, delay to warm to 300 K
MDsteps=18000000  # 36ns  ;    #100000000  # 200ns
MDps=$((MDsteps/500)) # dt = 0.002
#-----------------------------------------------------------

mkdir ${dirName}-$number
#topology and restraints
cp ${prevDirName}-$number/*.top ${dirName}-$number/   # topop.top
cp ${prevDirName}-$number/*.itp ${dirName}-$number/   # dihres.itp, disres.itp, posres.itp, ....
cp export.ndx ${dirName}-$number/   # export.ndx
cp $dihres ${dirName}-$number/
cp $disres ${dirName}-$number/


#structure
cp $firstStructureDir-$number/start*.gro ${dirName}-$number/start.gro
cp ${prevDirName}-$number/npt.gro ${prevDirName}-$number/npt.cpt ${dirName}-$number/ 
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

#stability
    cat <<< "
; Preprocessing
define        = -DDISRES -DDIHRES ; distance and dihedral restraints

; Run control
integrator    = md            ; leap-frog integrator
nsteps        = $MDsteps      ; 200 ns ; 1500 frames
dt            = 0.002         ; 2 fs

; Output control
nstxout        = 10000          ; save coordinates every 10.0 ps
nstvout        = 10000          ; save velocities every  10.0 ps
nstenergy      = 10000          ; save energies every    10.0 ps
nstlog         = 10000          ; update log file every  10.0 ps
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
ref_t         = $temp $temp    ; reference temperature, one for each group, in K

; Pressure coupling ON
pcoupl              = Parrinello-Rahman ; Pressure coupling on in NPT
pcoupltype          = isotropic         ; uniform scaling of box vectors
tau_p               = 1.0               ; time constant, in ps
ref_p               = 1.0               ; reference pressure, in bar
compressibility     = 4.5e-5            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com

; Simulated annealing
annealing           = single single     ; for each tc-grps 
annealing-npoints   = 3 3               ; # of points (temp) for each group ... 2 = starting, end temp.
annealing-time      = 0 ${j} $MDps 0 ${j} $MDps   ; ps, times when the given temp should be reach
annealing-temp      = $tempStart $temp $temp $tempStart $temp $temp  ; K, temps for each ann.-npoints

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
disre-fc        = 0          ; kJ / mol / nm2 * factor (weight)
disre-tau       = 0             ; ps
nstdisreout     = 10000          ; steps, write to energy file
orire           = no            ; orientation restraints
" > stability.mdp

# 1. submit
cat <<< '#!/bin/bash
# run as:    
# psubmit gpu_long md-run.sh ngpus=1,ncpus=1,mem=4gb,worksize=20gb,walltime=144h,props=^cl_konos 

metamodule add gromacs-5.0.5-gpu 

gmx grompp -f stability.mdp -c npt.gro -t npt.cpt -p topol.top -o stability.tpr
gmx mdrun -deffnm stability 

./convertTRR.sh
' > md-run.sh

# convert trr  
cat <<< "#!/bin/bash
#TRR postprocessing

metamodule add gromacs-5.0.5-gpu

gmx convert-tpr -s stability.tpr -o final.tpr -n export.ndx

#transform trajectory ... strip away solvent and remove jumps (based on first structure)
gmx trjconv -f stability.trr -n export.ndx -s start.gro -o tmp-noJump.trr -pbc nojump


gmx trjconv -f tmp-noJump.trr -s final.tpr -o last1.gro -b $MDps -boxcenter zero -center << EOF
0
0
EOF
# fit it to ... last frame

gmx trjconv -f tmp-noJump.trr -s last1.gro -o final.trr -fit rot+trans << EOF
0
0
EOF

#pdb
gmx editconf -f last1.gro -o ${name}.pdb

" > convertTRR.sh

chmod a+x md-run.sh
chmod a+x convertTRR.sh

cd ..



