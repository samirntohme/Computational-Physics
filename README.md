# Molecular Calculation from the First Principles (Ab Initio) with MOLPRO

## Overview
This repository contains a MOLPRO input script for quantum chemical calculations of the titanium nitride (TiN) molecule. The script is designed to compute electronic structure properties using multi-configurational self-consistent field (MCSCF) and configuration interaction (CI) methods, among others.

## Features
- **Geometry Specification**: The molecular geometry is defined in Cartesian coordinates (in Angstroms).
- **Basis Set Definition**: The script uses high-quality basis sets for titanium (Ti) and nitrogen (N) to ensure accurate calculations.
- **Effective Core Potential (ECP) for Ti**: To reduce computational cost, an ECP is applied to the titanium atom.
- **Wavefunction Optimization**: The script prints basis set details and orbital information to facilitate analysis.

## File Contents
- **MOLPRO Input File**: Contains commands for geometry setup, basis set selection, and electronic structure calculations.

## Dependencies
To execute this script, you need:
- **MOLPRO Quantum Chemistry Software**
- Adequate computational resources (recommended: at least 16 GB RAM and multiple CPU cores)

## Running the Calculation
1. **Ensure MOLPRO is installed** and accessible from your command line.
2. **Save the script** as `TiN.inp`.
3. **Run the calculation** using:
   ```sh
   molpro TiN.inp > TiN.out
   ```
   This command runs the MOLPRO job and redirects the output to `TiN.out`.

## Interpretation of Results
- The **SCF and MCSCF orbital information** will provide insight into the electronic structure.
- **Total energy values** indicate the stability of the molecule.
- **Dipole moments and transition dipole moments** can be extracted for spectroscopic analysis.

## Contact & Contribution
For any issues or improvements, please submit a pull request or contact the repository maintainer.

---
This repository serves as a reference for computational chemistry studies on TiN and can be adapted for similar transition metal-containing systems.

