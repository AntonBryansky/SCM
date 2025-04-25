# Spectral Collocation Method for dispersion curves plotting
Calculate guided wave dispersion curves for orthotropic (anisotropic) laminates.

## Description
The MATLAB-based script for calculation of dispersion curves for composite plates as the case of multilayer orthotropic media using (Pseudo-)Spectral Collocation Method.

## Papers
A. A. Bryansky, S. V. Panin, Dispersion curves calculation for all-CFRP sandwich composite using Spectral Collocation Method, Acta Astronautica (2025). [DOI](https://doi.org/10.1016/j.actaastro.2025.04.042)

## Files

`SCM_SL.m` - Spectral Collocation Method script for a single-layer case

`SCM_ML_HCS` - Spectral Collocation Method script for a multilayer case (example for honeycomb sandwich composite)

`chebdiff.m` - function for creation of the Chebyshev differentiation matrices

`balance2.m` and `gcg.m` - functions required for balancing algorithm (look [Important links](#important-links)).

`SM_HC.m` - function for calculation of the stifness matrix for honeycomb core

`SM_orthotropic.m` - function for calculation of the stiffness matrix for orthotropic materials (skin sheets)

`SM_isotropic.m` - function for calculation of the stiffness matrix for isotropic materials (adhesive layers)

## Literature
1. F. Hernando Quintanilla, M. J. S. Lowe, R. V. Craster, Modeling guided elastic waves in generally anisotropic media using a spectral collocation method, J. Acoust. Soc. Am. 137.3 (2015) 1180-1194. [DOI](https://doi.org/10.1121/1.4913777)
2. F. Hernando Quintanilla, M. J. S. Lowe, R. V. Craster, The symmetry and coupling properties of solutions in general anisotropic multilayer waveguides, J. Acoust. Soc. Am. 141(1) (2017) 406-418. [DOI](https://doi.org/10.1121/1.4973543)
3. M. Mekkaoui, S. Nissabouri, H. Rhimini, Towards an Optimization of the Spectral Collocation Method with a New Balancing Algorithm for Plotting Dispersion Curves of Composites with Large Numbers of Layers, J. Appl. Comput. Mech. 10(4) (2024) 801-816. [DOI](https://doi.org/10.22055/jacm.2024.45578.4390)
4. I. Zitouni, H. Rhimini, A. Chouaf, Comparative study of the spectral method, DISPERSE and other‎ classical methods for plotting the dispersion curves in‎ anisotropic plates, J. Appl. Comput. Mech. 9(4) (2023) 955-973. [DOI](https://doi.org/10.22055/jacm.2023.42530.3941)
5. A. Huber, M. G. Sause, Classification of solutions for guided waves in anisotropic composites with large numbers of layers, J. Acoust. Soc. Am. 144(6) (2018) 3236-3251. [DOI](https://doi.org/10.1121/1.5082299)

## Important links:
[Dispersion Calculator](https://github.com/ArminHuber/Dispersion-Calculator) by Armin Huber

About balancing algorithm for eigenvalue problem:
[BALANCE2 Balancing generalized eigenvalue problem](https://www.mathworks.com/matlabcentral/fileexchange/49719-balance2-balancing-generalized-eigenvalue-problem) (requires
[GCG Generalized conjugate gradient method](https://www.mathworks.com/matlabcentral/fileexchange/49720-gcg-generalized-conjugate-gradient-method))

## Thanks to
* Dr. Armin Huber, German Aerospace Center (DLR), Augsburg, Germany
* Prof. Pavel Utkin, Harbin Institute of Technology, Harbin, China
* Abir Rous (woah!)
