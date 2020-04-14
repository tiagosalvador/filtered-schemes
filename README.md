# Filtered Schemes
This repo contains code for the paper [Filtered schemes for Hamilton-Jacobi equations: a simple construction of convergent accurate difference schemes](http://dx.doi.org/10.1016/j.jcp.2014.12.039).
The idea of filtered schemes is to blend a stable monotone convergent scheme with an accurate scheme and retain the advantages of both: stability and convergence of the former, and higher accuracy of the latter.
The only requirement on the accurate scheme is consistency, thus allowing for a wide range of possibilities.

## Code
In order to use this code, start by running the script `compile_mex_files.m`.
This will compile the MEX-files needed for the solver.
The script `example_script.m` contains examples on how to call the different solvers and schemes.

## Citation
If you use this paper or code in your scientific work, please cite as
```
@article{filteredObermanSalvador,
  author    = {Adam M. Oberman and Tiago Salvador},
  title     = {Filtered schemes for Hamilton–Jacobi equations: A simple construction of convergent accurate difference schemes},
  journal   = {Journal of Computational Physics},
  volume    = {284},
  year      = {2015},
  url       = {https://doi.org/10.1016/j.jcp.2014.12.039},
  pages     = {367--388},
  publisher = {Elsevier}
}
```
