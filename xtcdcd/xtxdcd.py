#!/usr/bin/env python

import MDAnalysis as mda 

TPR = "npt.gro"
XTC = "md_0_1_noPBC.xtc"
u = mda.Universe(TPR, XTC)   # or use PDB or GRO instead of TPR 
with mda.Writer("2gaq.dcd", n_atoms=u.atoms.n_atoms) as W: 
   for ts in u.trajectory: 
        W.write(ts) 