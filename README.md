# Fortran-parallel-sort
A fortran library to perform parallel sorts.

## Install:
Download the library in any working directory.

## Compile:
You need to compile both `mrgrnk.f90` and `sort.f90`. To use, simply use `mod_sort` in your fortran files.

## Doc:
### parallel_sort:
- A (in): an integer array, containing the element to sort 
- order(out): an integer array, containing the positions so that `A(order(:))` is ordered
