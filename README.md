# Spectral Collocation Method for dispersion curves plotting
## Description
The main description will be uploaded later. 


## Files

`SCM_SL.m` - Spectral Collocation Method for a single-layer case

`SCM_ML_HCS` - Spectral Collocation MEthod for a multilayer case (honeycomb sandwich composite)

`chebdiff.m` - function for creation of the Chebyshev differentiation matrices

`balance2.m` and `gcg.m` - functions required for balancing algorithm (look [Important links](##Important-links)).

`SM_HC.m` - function for calculation of the stifness matrix for honeycomb core

`SM_orthotropic.m` - function for calculation of the stiffness matrix for orthotropic materials (skin sheets)

`SM_isotropic.m` - function for calculation of the stiffness matrix for isotropic materials (adhesive layers)

## Literature
1. Link 1
2. Link 2
3. Link 3
4. Will be added later...

## Important links:
[Dispersion Calculator](https://github.com/ArminHuber/Dispersion-Calculator) by Armin Huber

About balancing algorithm for eigenvalue problem:
[BALANCE2 Balancing generalized eigenvalue problem](https://www.mathworks.com/matlabcentral/fileexchange/49719-balance2-balancing-generalized-eigenvalue-problem) (requires
[GCG Generalized conjugate gradient method](https://www.mathworks.com/matlabcentral/fileexchange/49720-gcg-generalized-conjugate-gradient-method))
