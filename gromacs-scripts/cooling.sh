#!/bin/bash
# This script runs just once for all prepared files
#
# NOTE: due to a bug in distance restraints, only save option is to run MD on 1 CPU or 1 GPU... 
# NOTE2: The bug should be repaired in new Gromacs versions since year 2018.
# this script prepares file to be run in metacentrum with infinity on GPUs:


show_help() {
cat <<EOF
Prepare directories for cooling, 15 ns long MD, which is supposed to run on GPU (due to the bug  with distance restraints). Concretely, it works in Metacentrum with Infinity and runs on 1 GPU.

Usage: `basename $0` [-p parameters.txt] [-h] [-i]
    
    -i index file 'export.ndx' with the atoms of interest will be prepared for export from 15 ns trajectory (e.g. probably without solvent). Can reduce a data size a lot.
    -p parameters file name, default 'parameters.txt'
    -h help
EOF
}

trajExport=0
number=
input="parameters.txt"
OPTIND=1
while getopts "p:hi?" opt; do
    case $opt in
        p)
            input=$OPTARG
            ;;
        i)
            trajExport=1
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
prevDirName="`awk '$1 == "warming_directory" {print $2}' $input `" 
firstStructureDir="`awk '$1 == "equilibration_directory" {print $2}' $input `"

dihres="`awk '$1 == "dihedral_restraints" {print $2}' $input `"
disres="`awk '$1 == "distance_restraints" {print $2}' $input `"

tempCold="`awk '$1 == "cold_temperature" {print $2}' $input `"
tempHigh="`awk '$1 == "high_temperature" {print $2}' $input `"
#-----------------------------------------------------------
j=3000 # ps, delay at 1000K befor cooling to 0K
MDsteps=15000000
MDps=$((MDsteps/1000)) # dt = 0.001
#-----------------------------------------------------------


#prepare index file of what atoms should be exported from trajectory (e.g. non-water)
if [ $trajExport -eq 1 ]; then
    echo ""
    echo "Make sense to export only atoms common to all runs..."
    #echo "Otherwise, you have to create your index files with name 'export.ndx' for each topology."
    echo ""
    gmx make_ndx -f ${prevDirName}-1/warm.gro -o export.ndx  
fi


for i in ${prevDirName}-[0-9]* ; do
    # i is previous dir name
    k=`echo $i | tr -d '/' | grep -oE '[0-9]+$'`  # number of run
    mkdir ${dirName}-${k}

    cp $i/*.itp $i/*.top ${dirName}-${k}/
    cp $i/warm.gro  $i/warm.cpt ${dirName}-${k}/
    cp $dihres $disres ${dirName}-${k}/
    cp $firstStructureDir-${k}/${name}.gro ${dirName}-${k}/start.gro
    if [ -s export.ndx ]; then
        cp export.ndx ${dirName}-${k}/
    fi
    cd ${dirName}-${k}/

#cold.mdp
    cat <<< "
; Preprocessing
define        = -DDISRES -DDIHRES ; distance and dihedral restraints

; Run control
integrator    = md            ; leap-frog integrator
nsteps        = $MDsteps      ; 15 ns ; 1500 frames
dt            = 0.001         ; 1 fs

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
ref_t         = $tempHigh $tempHigh    ; reference temperature, one for each group, in K

; Pressure coupling ON
pcoupl              = no             ; Pressure coupling OFF

; Simulated annealing
annealing           = single single     ; for each tc-grps 
annealing-npoints   = 3 3               ; # of points (temp) for each group ... 2 = starting, end temp.
annealing-time      = 0 ${j} $MDps 0 ${j} $MDps   ; ps, times when the given temp should be reach
annealing-temp      = $tempHigh $tempHigh $tempCold $tempHigh $tempHigh $tempCold  ; K, temps for each ann.-npoints

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
nstdisreout     = 10000          ; steps, write to energy file
orire           = no            ; orientation restraints
" > cold.mdp

    # 1. submit
    cat <<< '#!/bin/bash
# run as:    
# psubmit gpu md-run.sh ngpus=1,ncpus=1,mem=2gb,scratch=6gb,walltime=24h,props=^cl_konos 

metamodule add gromacs-5.0.5-gpu 

gmx grompp -f cold.mdp -c warm.gro -t warm.cpt -p topol.top -o cold.tpr
gmx mdrun -deffnm cold 

./convertTRR.sh
' > md-run.sh


    # 2. submit, continue
    cat <<< "#!/bin/bash
# continue
# psubmit gpu md-run-append.sh ngpus=1,ncpus=1,mem=2gb,scratch=6gb,walltime=24h,props=^cl_konos

metamodule add gromacs-5.0.5-gpu 

#gmx grompp -f cold.mdp -c warm.gro -t warm.cpt -p topol.top -o cold.tpr

gmx mdrun -deffnm cold -cpi cold.cpt -append  

./convertTRR.sh
" > md-run-append.sh
  
    # convert trr  
    cat <<< "#!/bin/bash
#TRR postprocessing

metamodule add gromacs-5.0.5-gpu 

gmx convert-tpr -s cold.tpr -o final.tpr -n export.ndx

#transform trajectory ... strip away solvent and remove jumps (based on first structure)
gmx trjconv -f cold.trr -n export.ndx -s start.gro -o tmp-noJump.trr -pbc nojump


gmx trjconv -f tmp-noJump.trr -s final.tpr -o last.gro -b $MDps -boxcenter zero -center << EOF
0
0
EOF
# fit it to ... last frame

gmx trjconv -f tmp-noJump.trr -s last.gro -o final.trr -fit rot+trans << EOF
0
0
EOF

#pdb
gmx editconf -f last.gro -o ${name}.pdb

" > convertTRR.sh
    
    # CONTINUE, it brings data from killed job to continue with MD, needs Infinity at Metacentrum
    cat <<< '#!/bin/bash
machine=`grep -oP '\''(?<=INF_MAIN_NODE" value=")[^"]+'\'' md-run.sh.info`
directory=`grep -oP '\''(?<=INF_WORK_DIR" value=")[^"]+'\'' md-run.sh.info`
scp ${machine}:${directory}/* .
mkdir md-run1
mv md-run.sh.* ___JOB_IS_RUNNING___ md-run.sh md-run1/
cp md-run-append.sh md-run.sh
psubmit gpu md-run.sh ngpus=1,ncpus=1,mem=2gb,scratch=6gb,walltime=24h,props=^cl_konos <<< YES
' > continue.sh

    chmod a+x md-run.sh
    chmod a+x md-run-append.sh
    chmod a+x convertTRR.sh
    chmod a+x continue.sh
    
    cd ..

done

