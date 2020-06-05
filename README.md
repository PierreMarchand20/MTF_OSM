# Boundary integral multi-trace formulations and Optimised Schwarz Methods

This repository contains the code associated with:

Xavier Claeys and Pierre Marchand (2020). [_Boundary integral multi-trace formulations and Optimised Schwarz Methods_](https://hal.inria.fr/hal-01921113/document). Computers & Mathematics with Applications, June 2020, Vol. 79, Issue 11, Pages 3241-3256.

## Dependencies

C++:

1. CMake
2. Boost

Python:

1. scipy
2. pandas
3. matplotlib

For some reason, computing the spectrum with scipy did not always work and we used matlab in these cases.

## Figures

To reproduce the results of the article:

1. Compile the code
   1. `mkdir build & cd build`
   2. `cmake ../`
   3. `make`
2. Create the input data
   1. `cd ../scripts`
   2. `./gen_meshes.zsh` generates the meshes in `meshes`
   3. `./gen_mat.zsh` generates the matrices in `output/matrices`
3. Create the results
   1. `./gen_iterations.zsh` generate the data (number of iterations with number of interfaces )in `output/csv`
   2. `./gen_residus.zsh` generate the data (residual with number of iterations in GMRes )in `output/csv`
4. Create the figures
   1. `./gen_plot_iterations.zsh` generate the figures (number of iterations with number of interfaces )in `output/figures`
   2. `./gen_plot_residus.zsh` generate the figures (residual with number of iterations in GMRes )in `output/figures`

Concerning the spectrum:

- For Configuration II, you can use
  - `python3 compute_spectrum.py --ni 3 --geo non_emboite --type 1`
  - `python3 plot_spectrum.py --ni 3 --geo non_emboite --type 1`

where `ni` is the number of interfaces and `type` the setting for $\kappa$

- For Configuration I, you need to 
  - convert the matrices to matlab format with `python3 to_matlab.py --ni 3 --geo emboite --type 1`
  - run `compute_spectrum.m` with matlab
  - `python3 plot_spectrum.py --ni 3 --geo emboite --type 1`