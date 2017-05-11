#!/usr/local/bash

#The original code was taken from Bevan LAB from Virginia TECH
#Some modifications are by Mahmoud Agha
#No permission for copying the code was granted by Bevan lab
#The mere reason that we used thes configurations simply because they worked, we tried other things, but they were no good and wasted our time and effort

### Generating Topology
# Converting the pdb file to gmx
gmx pdb2gmx -f 2gaq.pdb -o 2gaq_processed.gro -water spce -ignh
# we chose 15 for OPLS

### Defining the Box (We use Periodic boundary conditions of a dodecahecron to save memory space and speed up calcultions)
# The protein is centered and 1.0 nm is enough to guarantee that it does not see its mirror images
gmx editconf -f 2gaq_processed.gro -o 2gaq_newbox.gro -c -d 1.0 -bt dodecahedron
#Solvate the protein with scp216
gmx solvate -cp 2gaq_newbox.gro -cs spc216.gro -o 2gaq_solv.gro -p topol.top

### Adding Ions
gmx grompp -f ions.mdp -c 2gaq_solv.gro -p topol.top -o ions.tpr
#The charge on the protein was a fraction and equilibrating it was hard, we settled for the value of 1, though the charge on the protein was around 1.3 ~ 1.4
gmx genion -s ions.tpr -o 2gaq_solv_ions.gro -p topol.top -pname NA -nname CL -np 1

### Energy Minimization
gmx grompp -f minim.mdp -c 2gaq_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
mkdir graphs                #We make a seperate file for graphs
gmx energy -f em.edr -o ./graphs/potential.xvg
# We chose 10 0 for potential and terminate
#We ran another energy minimization process
gmx energy -f em.edr -o ./graphs/potential2.xvg

### Equilibration
#temperature
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx energy -f nvt.edr -o ./graphs/temperature.xvg
#We chose 15 0 for temperature and terminate
#pressure
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
gmx energy -f npt.edr -o ./graphs/pressure.xvg
# we seleced 16 0 for pressure and terminate
gmx energy -f npt.edr -o ./graphs/density.xvg
#Type 22 0 for density and terminate

### Productino of MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -v    #-v to be verbose and display how long will it take
#The simulation did not run to full due to hardware malfunction and lack of advanced hardware equipment (The laptop was a mess acutally)

### Primary analysis
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -ur compact
# 0 for system output
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
# 4 for backbone
gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o ./graphs/gyrate.xvg




