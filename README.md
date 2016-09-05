# Protein Structure from assigned NMR spectra 

This project is mostly collection of simple BASH and AWK scripts, which helps to transform your data to correct input for next step of calculations... To explain, you usually need to combine several programs to reach your final goal (e.g., protein structure) and each program has different expectation for input - format, atom names, units... 

I went through this hell and prepared several scripts, which can help anyone in similar situation. Idea is to provide small simple scripts, that can be used separately and can be easily modify for anyone needs, because every one has own unique project with some specific hitch. 

## What can you find here

Starting point is an assignment of NMR (Nuclear Magnetic Resonance) spectra in the program [SPARKY](https://www.cgl.ucsf.edu/home/sparky/) or [NMRFAM-SPARKY](http://www.nmrfam.wisc.edu/nmrfam-sparky-distribution.htm) 

Scripts for...
* conversion to [CYANA](http://www.cyana.org/wiki/index.php/Main_Page) / XPLOR format. CYANA can make an automatic assignment of NOESY spectra, extract distance restraints, and calculate starting structure of protein. 
* conversion CYANA outputs to [GROMACS](http://www.gromacs.org), software for molecular dynamics. You should refine CYANA structures in proper force field.
* extracting results and running some statistics at [iCING](https://nmr.le.ac.uk/icing)
* conversions to NMRStar 3.1 format (.cif) for uploading your results to the databases [BMRB](http://www.bmrb.wisc.edu/) and [PDB](http://www.rcsb.org/pdb/home/home.do) 

## Future plans

* I am going to add examples. 

## Start with GIT:
[Practical short online course on CodeAcademy](https://www.codecademy.com/learn/learn-git)

