#!/usr/local/bin/gawk -f
#   nmr2gmx.awk
#@      make gromacs topology files for distance restraints 
#@      and dipolar coupling data                                          
#
#
# to run do one of: 
# gawk -f nmr2gmx.awk args top-or-gro.file   restraint-files ...
# gawk -f `which nmr2gmx.awk` args top-or-gro-file  restraint-files ......
# nmr2gmx.awk args files..... [ if you have /usr/local/bin/mawk ]
#
# Command-line options:
#
# Ambiguous Distance Restraining options:
#
# -v lengthening=0   : NOE distances will be extend by x angstrems (default 0)
#
# -v wtype=R0  : default
# -v wtype=R-1
# -v wtype=R-3
# -v wtype=R-6
# -v wrefdist=3.0 [default] set reference distance Dref 
#	 weight = (D/wrefdist)^-1 for wtype  =R-1 	
#
#  # dup : set fixed distamce between up0 and up1 , force is constant if d>up1
#          default not set 
# -v dup=0.1        : difference between up0 and up1 in nanometer 
# -v dup=0.1nm      : difference between up0 and up1 in nanometer
# -v dup=1angstr  : difference between up0 and up1 in angstrom
# 
# XV: cross-validation stuff ; default XV=0
# -v XV=1       : make files for cross-validation
# -v nsplit=10  : number of groups to split data in for XV 
# -v seed=1999  : seed for random number generator [gawk/mawk/awk/awka might vary]
#
#
# -v ensfunc=[1|2]  
#      ensfunc=1 :  normal  
#      ensfunc=2 :  no ensemble averaging [in case of ensemble]
#
# 
# -v weight=1 [default] :  value of weight of restraint 
#
# -v rsttype=1 [default]  restraint type , currently only 1 
#                        [ 3rd number in top-tile ]
#
# -v resoffset=i [default=0]   
#                 offset to match residue number in nmr data file and 
#                 gromacs topology 
#                 nmr-res-nr = gmx-res-nr + resoffset  
#              so:  resoffset is (pdb) res number at which gmx structure starts.- 1
# 
# the script needs the atom numbers ,it can reed it from topology and .gro
# 	coordinate files.
#
# restraint format [ free , whitespace separated ]
# constraint res-nr res-nam atom  res-nr res-nam atom upperbound lowerbound nr
# e.g. :
# constraint 5 ala ha 8 leu  hg  4.5  3.5  1
# or for ambiguous restraints 
# constraint 5 ala hb1 8 leu  hg  4.5  3.5  3
# constraint 5 ala hb2 8 leu  hg  4.5  3.5  2
# constraint 5 ala hb3 8 leu  hg  4.5  3.5  1
#
# there is experimental support for cyana type .upl files 
#
# the following can also be read.
# constraint 5 ala hb# 8 leu  hg  4.5  3.5  1
#
# dipolar coupling data can be read  in the following form:
# dipolar coupling data for phospho-iib-chb
#    _not_ used for pdb file : 1h9c.pdb
#    all couplings in hz
#    data without reference value : err = 0.8 hz
#
# #key  nr1 res1 atom1  nr2 res2 atom2  coupling    error   dcr-nr
# dip     5  res  hn         5 res  n      -12.406    0.026     5
# dip     6  res  hn         6 res  n       -0.263    0.027     6
# dip     7  res  hn         7 res  n        3.055    0.800     7
# dip     8  res  hn         8 res  n       10.008    0.801     8
#
# end of help

BEGIN {
	if (!stderr) stderr="/dev/stderr"
    progname="nmr2gmx.awk"
    if ( ARGV[1]~/^(h|help)$/ || help!="" ) { 
        print ARGV[0]
        # because of the which command, help is only available if 
        # this script is in the PATH
        system("gawk '/^#/{print}/^BEGIN/{exit}' `which "progname"`")
        help=1
        exit
    }
    IGNORECASE=1
    if (ARGC==1) help=1
    for ( i in ARGV ) if (ARGV[i]=="help") help=1
    if (help) { 
        system("grep '^#' `which nmr2gmx.awk`") 
        exit
    }
    
    top="top";gro="gro";itp="itp"
    
    if (!lengthening) lengthening=0
	if (!rsttype) rsttype=1
	if (!weight) weight=1
	if (!dipco_viol_est) dipco_viol_est=1/2 # Hz
	if (!type) type="upl"
	if (type=="lol") {
		print "Lower bounds not yet implemented in Gromacs??" > stderr
		print "expecting Angstroms in input file" > stderr
	}

    pwd=ENVIRON["PWD"]
    
	if (!wrefdist) wrefdist=0.3 # 3.0 Angstrom
	
	# distance larger than maxup nm. are considered lowerbounds (nonnoes)
	if (!maxup) maxup=2  # distance larger than maxup nm. are considered lowerbounds (nonnoes)
	# XVAl : Cross Validation stuff
	if (!XV ) noXV=1

	else { 
		if (!nsplit) nsplit=10
		if (!seed) seed=1999
		srand(seed)
		print "# generating random numbers for XV [cross-validation]" > stderr
		for (ra=0;ra<=1000;ra++) ran[ra]=rand()
		print "# -done"> stderr
	}

	if (dgxvout ) "print "nsplit"constraint-files" 
	
	# make list of fixed pairs hb1-hb2 etc
	fixedpairs["hb1"] = "hb2"
	fixedpairs["hg1"] = "hg2"
	fixedpairs["hg11"]= "hg12"
	fixedpairs["hd1"] = "hd2"
	fixedpairs["he1"] = "he2"
	for ( p in fixedpairs ) fixedpairs[fixedpairs[p]]=p
	
	mult=0;	prevamb=1
}

#
#  residue names are disregarded!!
#

FNR==1 {
	if (FILENAME~/\.top$/) { ftype=top; print "# top:",topfile=FILENAME > stderr}
	if (FILENAME~/\.gro$/) { ftype=gro; print "# gro:",grofile=FILENAME > stderr}
	if (FILENAME~/\.itp$/) { ftype=itp; print "# itp:",itpfile=FILENAME > stderr}
	print  ";  collection atom indices ... : "
	printf ";  reading topology info from :%s\n",pwd"/"FILENAME
}

/^#/ {
    if (FILENAME!=topfile && FILENAME!=grofile && FILENAME!=itpfile) 
        print ";",$0
    next
}

 # dipolar coupling restraints
/^dip / {
	if (!topfile && !grofile ) { 
		print "no atom indices, give *.gro or *.top file first"> stderr
		exit
    }

	if (!(FILENAME in cread)) { 
		cread[FILENAME]=1
		printf "\n;  topology-file created by nmr2gmx.awk :%s\n",FILENAME
		printf "\n;  dipolar-coupling values read from :%s\n",FILENAME
	}
	if (!dipconstraintfile) {
		print "# reading dipolar coupling values from :",pwd"/"FILENAME > stderr
		
		dipconstraintfile=FILENAME

        # mawk 'BEGIN{ h=599.6765403;printf "%12.8f %12.8f \n",60.7715167/h,150.8140440/h }'
        #  0.10134049   0.25149232 
        gamma["h"] = 1
        gamma["c"] = 0.25149232
        gamma["n"] = 0.10134049
        K = 24.18761733
#
#  K = mu0/(4pi )  h_ /(4pi) = 24.18761733 
#  const = K gamma1 gamm2
#  const[H-H]= 24.1876
#  const[H-N]=  6.083
#  const[H-C]=  2.45118
#  const[N-C]=  0.616454
# 

		print "\n[ orientation_restraints ]"
		print ";  ftype : function type = 1"  
		print ";  exp   : label for the experiment [one consistent ordering tensor]"  
		print ";  label : label for the specific coupling value "  
		print ";  power : power used for distance dependency usually r=0 or r=3"  
		print ";  const : constant depending on nuclei:"    
		print ";          const = mu0/(4pi ) * gamma1 * gamma2 * h_ /(4pi)"    
		print ";          const = K * gamma1 * gamma2 ; K = "K" ; gammaH==1"    
		print ";          const = 6.083 Hz nm^3   for N-H vector"    
		print ";  dipco : dipolar coupling in Hz"    
        print ";  weight: ~ (v^2)/(v^2+s^2)"
        print ";          v = violation estimate =",dipco_viol_est,"Hz"
        print ";          s = estimated error in measured coupling value"
		print ";  ai    aj ftype exp label power  const    dipco   weight ;& source"  
		#  
	}
    
	# generate proper $0
    source=$0
    comment = ( $0~"#" ? gsub(/^.*#/,"") : "" )
    # remove comment from $0
    gsub(/#.*/,"")
    $0=tolower($0)

	$1=$1; # shrink whitespace 
    
    $4 = rename_atom($4,$3)
    $7 = rename_atom($7,$6)

	# check if atom names exist
	if ( !($2" "$4 in num )) { 
		printf "# %4s %4s %4snot found in [%s]\n",$2,$3,$4,$0 > stderr
	    next
	} 
	if ( !($5" "$7 in num )) { 
		printf "# %4s %4s %4s not found in [%s]\n",$5,$6,$7,$0 > stderr
        next
	}

	if ( $2$4==$5$7 ) {
		selferr++
		print "ERROR in dcr : same atoms"$2$4"=="$5$7"\n["$0"]" > stderr 
		exit
	} 

	dico=$8+0  # Angstr
    dico_err=$9+0   # error estimate in same units as coupling
	dico_weight=dipco_viol_est^2/(dipco_viol_est^2+dico_err^2)
    #dico_nr=$10+0  # Angstr
	
	functype=1  	    # default dcr-restraints
    if (!experiment_nr) experiment_nr = 1 
    
    max_restr_id = ( max_restr_id < restr_id ? restr_id : max_restr_id )
    
    #|| (restr_id in idlist)
    next_id=max_restr_id+1
    if ( $10!="" && !(nextid in idlist) ) 
        restr_id = $10
    else if (!($2 in idlist)) 
        restr_id = $2
    else
        restr_id=next_id
    
    idlist[restr_id]++

    #  power for r-dependence 
    #    r=0 [fixed] or  
    #    r=3 not fixed both for dipolar coupling expts
    if (pow=="") pow = 3   
    
    # const = mu / (4 pi ) g1 g2 h_ /4
    #   depends on nuclei

    nucl1_type=nucl2_type=""
    if (match($4,/[a-z]/)>0) nucl1_type=substr($4,RSTART,RLENGTH)
    if (match($7,/[a-z]/)>0) nucl2_type=substr($7,RSTART,RLENGTH) 
    if (nucl1_type=="") { print "error: gamma missing for :",nucl1_type; exit} 
    if (nucl2_type=="") { print "error: gamma missing for :",nucl1_type; exit} 
    #print nucl1_type,gamma[nucl1_type],gamma[nucl2_type],nucl2_type
    const=K * gamma[nucl1_type] * gamma[nucl2_type]
    #    const = 6.083 # in  Hz nm^3 :  mu0/(4 pi ) * gamma1 gamma2 * h_/(4pi) 

    dcrwr++

 	printf "%5i %5i %5i %3i %5i %5i %6.3f %8.4f %8.5f ; # %s\n",
		num[$2" "$4],num[$5" "$7],functype,experiment_nr,restr_id,
		pow,const,dico,dico_weight,source

	if (dgxvout && XV) print source > FILENAME"_"tag
	#if ( $10==1 && XV ) print "#endif"
    	
	next
}

 # read constraintfile
/^constraint / { $0=tolower($0)
	#print "#["FNR":"$0"]">stderr
    if (!topfile && !grofile ) { 
		print "no atom indices, give *.gro or *.top file first"> stderr
		exit
	}
	if (!(FILENAME in cread)) { 
		cread[FILENAME]=1
		printf "\n;  dg90-distance-constraints from :%s\n",pwd"/"FILENAME
		printf "; Wtype = :'%s'\n",wtype
		printf "; type  = :'%s'\n",type
		printf "; dup   = :'%s'\n",dup
	
	}
	if (!constraintfile) {
		
		constraintfile=FILENAME

		print "# reading dg90-constraints from :",constraintfile > stderr

		print "\n[ distance_restraints ]"
		print ";   ai    aj  type1 index type2 low  up1 up2  fac  ; & source"  
		#  
	}
	
	# check constraints 
	if ($10==1) ccn++ 

	curramb=$10 ; curline=$0
	
	#printf "#[%s][%s][%s]\n",prevamb,prevamb-1,curramb
    if ($4=="-" || $7=="-") next

	if ( prevamb!=1 && prevamb-curramb!=1 ) {
		amberr=1
		print "ERROR in constraint ambiguity",prevamb"-1 != "curramb > stderr 
		print " previous line:["prevline"]" > stderr 
		print "  current line:["$0"]" > stderr 
		exit
	} 
	# skip fixed distances and diagonals
	if ( fixwarn ) {
		if ( $2==$5 ) {
			if ( $4==$7 ) print "Warning : self restraint ",$0 > stderr  
			if ( $4 in fixedpairs ) if ( fixedpairs[$4]==$7 )
				print "Warning : fixed distance ",$0 > stderr
			if ( $3~/^(phe|tyr)/ && $4~/^h[dez]/ && $7~/^h[dez]/ )
				print "Warning : fixed distance ",$0 > stderr
		}
	}
	
	if (mult==0) { mult=ccn ; countc[++c]=$0 }

	if (prevamb==1) {
		upperbound=$8+lengthening  # Angstr
		lowerbound=$9+0  # Angstr
	}

	#if (type=="lol") bound=$9
	#weight=1/mult
	source=$0
	$1=$1; # shrink whitespace 
	
    $4 = rename_atom($4,$3)
    $7 = rename_atom($7,$6)
    
    # check if atom names exist
    
    if ( nm1=match_atom($2-dres,$4,num,att1,nmat1) ) {
    } else {
		printf "# %4s %4s %4s not found in [%s]\n",$2,$3,$4,$0> stderr
        next
    }
    if ( nm2=match_atom($5-dres,$7,num,att2,nmat2) ) {
    } else {
		printf "# %4s %4s %4s not found in [%s]\n",$5,$6,$7,$0 > stderr
        next
    }
    #for ( atid in att1 ) { natt1[atid]=num[atid] ; print ">>1",natt1[atid],"["$2" "$4"]" }
    #for ( atid in att2 ) { natt2[atid]=num[atid] ; print ">>2",natt2[atid],"["$5" "$7"]" }
    #for ( i in nmat1 ) print "]"i,nmat1[i]
    #for ( i in nmat2 ) print "]"i,nmat2[i]
	#asort(natt1)
    #asort(natt2)

	if ( $2$4==$5$7 ) {
		selferr=1
		print "ERROR in constraint : same atoms" > stderr 
		print "["FNR":"$0"]" > stderr 
		exit
	} 

	
	functype=1  	     # default  NOE-restraints
	low=lowerbound*0.10  # Angstrom -> nanometer
	UP0=upperbound*0.10  # Angstrom -> nanometer
    # was: UP1=UP0+0.1*0.4/UP0
	if (dup) { 
	 	if (dup~/nm$/) dup+=0
	 	else if (dup~/angstr$/) dup/=10
		else dup+=0 #default nm		
		UP1=UP0+dup              # was UP1=UP1+dup 
	} else UP1=UP0+0.1*0.4/UP0   
	
    
    
	fac=weight     
	# hhmm? next is superfluous
    if (!wtype ||wtype=="R0") { 
        fac=weight 
	}
    
    #wrefdist=0.3 # 3.0 Angstrom ## moved up
	if (wtype=="R-1") { 
		fac=(UP0/wrefdist)^(-1) ; UP1=UP0+0.1;fac=sprintf("%.4f",fac) }
	if (wtype=="R-3") { 
		fac=(UP0/wrefdist)^(-3) ; UP1=UP0+0.1;fac=sprintf("%.4f",fac) }
	if (wtype=="R-6") { 
		fac=(UP0/wrefdist)^(-6) ; UP1=UP0+0.1;fac=sprintf("%.4f",fac) }

	#
    # setting functype
    # functype =2 ; keep restraint for each struc of ensemble
	#               thus no ens-averaging in this case for 
    #               H-bonds
	# in case of hbonds [ one of the atoms is oxygen ]
	if ( $4~/^o/ || $7~/^o/ ) {
		fac=1
		functype=2 
	}

	if ( ensfunc ) functype=ensfunc  
	#if ( setens ) if ( $0 !~ setens ) functype=2
	
    # check if it is a lower bound [ if upperbound > 20 angstr ]
    # not really necessary
    islow = ( upperbound > UP0 ? "1" : "" )
	if ( prevamb==1 ) { cwr++ ; countcwr[c]=$0 } 
	#fc = (islow ?  ";" : " " )

	if ( $10==mult ) { 
		# print to itp file with or without XV ( Cross Validation )
		if (XV) printf "#ifndef SKIP_%s\n",tag=randtag(cwr-1)
		if (XV ) xg="xv-group "tag
		if (mult>1) printf "; restraint %i multiplicity %i %s\n",cwr-1,mult,xg
	}
    #print nm1,nm2
    for ( i=1;i<=nm1;i++) for ( j=1;j<=nm2;j++) {
	    printf fc"%5i %5i %4i %4i %4i %6.4f %6.4f %6.4f %s ; # %s\n",
		    nmat1[i],nmat2[j],rsttype,cwr-1,functype,
		    low,UP0,UP1,fac,(i*j==1? $0 : "" )
    }
    
    if (XV) {
        if ( dgxvout ) print source > FILENAME"_"tag
	    if ( $10==1  ) print "#endif"
    }
    
	if (UP0==0) { print "up=0" > stderr ; exit}

	if (curramb==1) mult=0
	prevamb=curramb+0
    prevline = curline
	
	next
}


# read cyana .upl
FILENAME ~ /[.]upl$/ { 
    if (echo) print "#",$0
    $0=tolower($0)
	if (!topfile && !grofile ) { 
		print "no atom indices, give *.gro or *.top file first"> stderr
		exit
	}
	if (!(FILENAME in cread)) { 
		cread[FILENAME]=1
		printf "\n;  cyana-distance-constraints from :%s\n",pwd"/"FILENAME
		printf "; Wtype = :'%s'\n",wtype
		printf "; type  = :'%s'\n",type
		printf "; dup   = :'%s'\n",dup
	
	}
    
	if (!constraintfile) {
		
		constraintfile=FILENAME

		print "# reading cyana-constraints from :",constraintfile > stderr

		print "\n[ distance_restraints ]"
		print ";   ai    aj  type1 index type2 low  up1 up2  fac  ; & source"  
		#  
	}
	
	# check constraints 
	if ( $7+0>0 ) { 
        ccn++ 
		upperbound= $7 + lengthening  # Angstr
		lowerbound= $7 + 0  # Angstr
    }

	curline=$0
	
	# skip fixed distances and diagonals
	if ( fixwarn ) {
		if ( $1==$4 ) {
			if ( $3==$6 ) print "Warning : self restraint ",$0 > stderr  
			if ( $3 in fixedpairs ) if ( fixedpairs[$3]==$6 )
				print "Warning : fixed distance ",$0 > stderr
			if ( $2~/^(phe|tyr)/ && $3~/^h[dez]/ && $6~/^h[dez]/ )
				print "Warning : fixed distance ",$0 > stderr
		}
	}
	
	source=$0
	$1=$1; # shrink whitespace 
	
    
#    $3 = rename_xeasy_atom($3,$2)
#    $6 = rename_xeasy_atom($6,$5)
    
    # check if atom names exist
    
    if ( nm1=match_atom($1-dres,$3,num,att1,nmat1) ) {
    } else {
		printf "# %4s %4s %4s not found in [%s]\n",$1,$2,$3,$0> stderr
        next
    }
    if ( nm2=match_atom($4-dres,$6,num,att2,nmat2) ) {
    } else {
		printf "# %4s %4s %4s not found in [%s]\n",$4,$5,$6,$0 > stderr
        next
    }

	if ( $1$3==$4$6 ) {
		selferr=1
		print "ERROR in constraint : same atoms" > stderr 
		print "["$0"]" > stderr 
		exit
	} 

	
	functype=1  	     # default  NOE-restraints
	low=lowerbound*0.10  # Angstrom -> nanometer
	UP0=upperbound*0.10  # Angstrom -> nanometer
    # was: UP1=UP0+0.1*0.4/UP0
	if (dup) { 
	 	if (dup~/nm$/) dup+=0
	 	else if (dup~/angstr$/) dup/=10
		else dup+=0 #default nm		
		UP1=UP0+dup              # was UP1=UP1+dup 
	} else UP1=UP0+0.1*0.4/UP0   
	
	fac=weight     
	# hhmm? next is superfluous
    #if (!wtype ||wtype=="R0") { 
    #    fac=weight 
	#}
    
    #wrefdist=0.3 # 3.0 Angstrom ## moved up
	if (wtype=="R-1") { 
		fac=(UP0/wrefdist)^(-1) ; UP1=UP0+0.1;fac=sprintf("%.4f",fac) }
	if (wtype=="R-3") { 
		fac=(UP0/wrefdist)^(-3) ; UP1=UP0+0.1;fac=sprintf("%.4f",fac) }
	if (wtype=="R-6") { 
		fac=(UP0/wrefdist)^(-6) ; UP1=UP0+0.1;fac=sprintf("%.4f",fac) }

	#
    # setting functype
    # functype =2 ; keep restraint for each struc of ensemble
	#               thus no ens-averaging in this case for 
    #               H-bonds
	# in case of hbonds [ one of the atoms is oxygen ]
	if ( $3~/^o/ || $6~/^o/ ) {
		fac=1
		functype=2 
	}

	if ( ensfunc ) functype=ensfunc  
	
	if ( $7+0>0 ) { 
		# print to itp file with or without XV ( Cross Validation )
		if (XV) { 
            printf "#ifndef SKIP_%s\n",tag=randtag(cwr-1)
		    xg="xv-group "tag
        }
		if (mult>1) printf "; restraint %i multiplicity %i %s\n",cwr-1,mult,xg
	}

    cwr=ccn
    for ( i=1;i<=nm1;i++) for ( j=1;j<=nm2;j++) {
	    printf fc"%5i %5i %4i %4i %4i %6.4f %6.4f %6.4f %s ; # %s\n",
		    nmat1[i],nmat2[j],rsttype,cwr-1,functype,
		    low,UP0,UP1,fac,(i*j==1? $0 : "" )
    }

    if (XV) {
        if ( dgxvout ) print source > FILENAME"_"tag
	    if ( $10==1  ) print "#endif"
    }

	if (UP0==0) { print "up=0" > stderr ; exit}

	if (curramb==1) mult=0
	prevamb=curramb+0
    prevline = curline

	next
}


FILENAME==topfile || FILENAME==itpfile {
	if (/[[] .* \]/) {  
		itype=$2 ; 
		if (verbose!=0) print FNR,$0,itype > stderr
        next
	}
    if ($1==";") next
	if (itype=="atoms" && NF>5 ) {
		atomname=tolower($5)
		resname=tolower($4)
		resnum=tolower($3)
	    atomnum=$1+0
		# 1met h1   -> 1 met h
		if (resnum==1 && atomname=="h1" ) {
			atomname="h"
			print "# 1 h1 renamed to 1 h" >  stderr
		}

        #addjust residue offset
        resnum+=resoffset

		num[resnum" "atomname]=atomnum
		name[atomnum]=resnum" "atomname
        if (logfile) print atomnum,resnum,atomname > "top.log"
		next
	} else 
        next
}

 # read coodinates
FILENAME==grofile && /^ +[0-9\-]+[0-9A-Z\-\+]+ +[A-Z0-9]+ +[0-9]+/ {
	if (!gro++) print "reading gromos coordinates from:",FILENAME > stderr
	if (FILENAME!~/\.gro$/) print "Warning no gro file "> stderr
	
	resnum=substr($0,1,5);sub(/^ +/,"",resnum)
    
	resname=tolower(substr($0,6,5));sub(/ +$/,"",resname)
	atomname=tolower(substr($0,11,5));sub(/^ +/,"",atomname)
	atomnum=substr($0,16,5)+0
	
	# 1met h1   -> 1 met h
	if (resnum==1 && atomname=="h1" ) {
		atomname="h"
		print "# 1 h1 renamed to 1 h" > stderr	
	}
    #addjust residue offset
    resnum+=resoffset
	
	num[resnum" "atomname]=atomnum
	name[atomnum]=resnum" "atomname
    if (logfile) print atomnum,resnum,atomname > "gro.log"
	next
}

NF>0 {print "not processed: "$0>stderr ; next}

END {
    if (help)exit
    if(curramb!="" && curramb!=1) {
		amberr=1
		print "ERROR in constraint ambiguity;",prevamb-1,"more lines expected"> stderr 
		print "  after line:"FNR"["curline"]" > stderr 
		exit
	} 

	if ( selferr ) print "encountered selferr [restraint to same atom]" 
	if ( amberr ) print  "encountered ambiguity erorr"
	if (dipconstraintfile) print "# orientation restraints written :",dcrwr > stderr 
	if (constraintfile) { 
        print "# distance constraints read :",ccn+0 > stderr 
	    print "#   lowerbound constraints skipped :",lbc+0 > stderr 
	    print "#   constraints written :",cwr+0 > stderr 
	    print "#   check:total",ccn+0"-"lbc+0"="ccn-lbc"="cwr+0 > stderr 
    }
	for ( i in countc ) if (!countcwr[i]) print i,countc[i] > stderr 
	for ( i in countcwr ) if (!countc[i]) print i,countcwr[i] > stderr 
	
}

function rename_atom(atom,res) {
	
    renamed=0
    outatom=atom
    # renaming
    if (atom=="hn") { outatom="h" ; renamed++ }
    if (atom~"^hd1[123#%*]$" && res~/^ile$/) {
        sub(/1/,"",outatom) ; renamed1++ 
    }
    if (renamed && !action["rename"]++ ) {  
        printf "# renaming atoms in %s to %s ; rest silent..\n",FILENAME,$2" "$4 > stderr
    }

    return outatom
}
function rename_xeasy_atom(atom,res) {
	
    renamed=0
    outatom=atom
    # renaming
    if (atom=="hn") { outatom="h" ; renamed++ }
    if (atom~/^q/)  { 
        sub("q","h",outatom)
        outatom=outatom"#"
        renamed++ 
    }
    if (atom~"^hd1[123#%*]$" && res~/^ile$/) {
        sub(/1/,"",outatom) ; renamed1++ 
    }
    if (renamed && !action["rename"]++ ) {  
        printf "# renaming atoms in %s to %s ; rest silent..\n",FILENAME,$2" "$4 > stderr
    }

    return outatom
}

function randtag(nseed,   id) {
 # global nsplit 
 # global ran
	if (!( nseed in ran )) { 
		#print "restarting rand" > stderr
		for ( ra=ra ; ra<=nseed+1000 ;ra++) ran[ra]=rand()
	}
	#print "; rand "ran[nseed]" seed "nseed
	return int(ran[nseed]*nsplit)+1
}

function match_atom(resnum,atom,num,matches,nmatches,nm,i) {
    #for ( atid in matches ) delete matches[atid]
    #for ( atid in nmatches ) delete nmatches[atid]
    delete matches ; delete nmatches
    
    nm=0
    if  (resnum" "atom in num) { 
        #print resnum" "atom,"matched"
        matches[resnum" "atom]++
    } else if (atom~/[%*#]$/) { 
        qexp="^"resnum" "atom;sub(/[%*#]$/,"[123]$",qexp)
        for ( atid in num) if (atid~qexp) matches[atid]++ 
    }
    for ( atid in matches ) {
        nmatches[++nm]=num[atid];
        if (verbose) print "#",nm,resnum" "atom,"matched to ",atid 
    }
#    printf "in"; 
#    for (i in nmatches) printf " %i:%i",i,nmatches[i] ; printf "\n"
    asort(nmatches)
#    printf "in"; 
#    for (i in nmatches) printf " %i:%i",i,nmatches[i] ; printf "\n"
    
    if (nm && verbose) print "### total matches:",nm
    return nm
}
