# **Mixed Coarse-Grained Anisotropic Network Model (mcgANM)**

This repository contains `mcganm.f` script for **Mixed Coarse-Grained Anisotropic Network Model (mcgANM)**, which describes the protein structure as an elastic network, formed of nodes connected by elastic springs (Kurkcuoglu et al., 2009). Nodes can be placed at the center of mass (coarse-grained, low-resolution region) or the heavy atoms (atomistic, high resolution region) of the residues/nucleotides. The neighboring nodes (residues/nucleotides) are connected if their heavy atoms are within a distance of 7-10 Å. The strength of the elastic spring is proportional to the number of heavy atom contacts.

Kurkcuoglu O, Turgut OT, Cansu S, Jernigan RL, Doruker P (2009), Focused Dynamics of Supramolecules by Use of a Mixed-Resolution Elastic Network Model, Biophysical Journal, 97(4), 1178-1187

## 📖 **Overview**  
-	**High-resolution regions** are modeled with atomic detail.
-	**Low-resolution regions** are coarse-grained to reduce computational cost.
-	**Hessian matrix** is constructed, and **normal mode analysis** is applied.
-	**Principal Axes Calculation**: Computes the principal axes of the molecule for structural analysis.
-	**Mode Analysis**: Identifies slowest and fastest modes and calculates residue fluctuations.
-	**B-Factor Calculation**: Computes theoretical B-factors and compares them with experimental values.

## 🛠 Dependencies

- **Fortran Compiler** (e.g., `gfortran`) → Required for `mcganm.f`.  



## 🚀 **Usage**  

**Input:**  
- PDB file containing atomic coordinates (e.g., `7aku-drugA.pdb`).

**Key Outputs:**  
- `coordinates_com.out`: Coordinates relative to the center of mass.
- `principal_axes2.out`: Principal axes of the protein structure.
- `min_coordinates.out`: Transformed coordinates according to the principal axes of the protein.
- `eigenvalues.out`: Eigenvalues sorted from highest to lowest. (6 eigenvalues should be approximately equal to zero for convergence)
- `bfactors.out`: Experimental (2nd column) and theoretical B-factors (4th column).
- `mod1a.pdb`, `mod1b.pdb`, `mod2a.pdb`, `mod2b.pdb`, `mod3a.pdb`, `mod3b.pdb`, `mod4a.pdb` and `mod4b.pdb`: Distorted coordinates for visualization of normal mode displacements at the first four slowest modes.
- `slowest.out`: Residue fluctuations at the 20 slowest modes.
- `slowest_cumulative.out`: Residue fluctuations averaged over the 20 slowest modes.
- `fastest_cumulative.out`: Residue fluctuations averaged over the 20 fastest modes.
- `first100modes.out`: Displacement vectors for the first 100 modes.

**Key Parameters that should be specified in the script:**  
- `rcutt`: Cutoff distance for residue interactions (default: 10 Å, it can be taken 7-10 Å). 
- `rno`: Number of residues in the system.  
- `nresidue`: Number of atoms in the system.  
- `cg and cgx`: Residue indices defining the coarse-grained and high-resolution regions. If cg is taken to be equal to cgx, the protein structure is coarse-grained, with one node per residue.

**Command-line Usage:**  
```bash
# Compile
gfortran -o mcganm mcganm.f  

# Run
./mcganm
```

(Ensure the input PDB file is in the same directory.)


## 📂 **Example Input & Output**

The example_data/ directory provides sample files:

   * 7aku-drugA.pdb (example PDB structure).
   * mcganm.f (Fortran script).
   * Output files.

  
## 📜 **License**

This project is licensed under the **MIT License**.
See the `LICENSE` file for details.

## 📌 **Citation**

If you use this work in your research, please cite this study and the associated publication:

Kurkcuoglu O, Turgut OT, Cansu S, Jernigan RL, Doruker P (2009), Focused Dynamics of Supramolecules by Use of a Mixed-Resolution Elastic Network Model, Biophysical Journal, 97(4), 1178-1187
DOI: 10.1016/j.bpj.2009.06.009
  
## ✉️ **Contact**

For questions or issues, please contact:

Dr. Ozge Kurkcuoglu — 📧 olevitas@itu.edu.tr

