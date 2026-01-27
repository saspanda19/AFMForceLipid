# AFMForceLipid

Python scripts used to construct a coarse-grained (CG) amorphous silica AFM tip and to analyze AFM tip–lipid membrane interaction simulations performed using GROMACS and the MARTINI 3 force field.

**Requirements:**  

Python3 or greater
MDAnalysis 2.6.1 or greater

**Repository Structure:**  

AFMForceLipid/
├── AFM Tip/
│   ├── cgsilica.py        # Build CG amorphous silica from atomistic structure
│   └── tipshaping.py      # Generate a square pyramidal AFM tip geometry
├── Analysis/
│   ├── force_analysis.py     # Analyze AFM tip force data (.xvg)
│   └── position_analysis.py  # Analyze AFM tip position data (.xvg)








