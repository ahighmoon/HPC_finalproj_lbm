# Kármán vortex street simulation using Lattice Boltzmann Method

Final project for NYU Graduate course MATH-GA 2012 Advanced Topics in Numerical Analysis: High Performance Computing.

Collaborative work of Yifei Zhu, Tiffany Li, Kitty Li.

## Usage

To compile type `make`. Input parameter and obstacle files are all specified on the command line of the `d2q9-bgk` executable.

Usage:

    $ ./d2q9-bgk <paramfile> <obstaclefile>
eg:

    $ ./d2q9-bgk input_256x128.params obstacles_256x128.dat

## Script for running on Greene cluster

    $ sbatch d2q9-bgk.sbatch
    $ sbatch cuda-d2q9-bgk.sbatch

If you wish to run a different set of input parameters, you should modify `d2q9-bgk.sbatch` or `cuda-d2q9-bgk.sbatch.

## Result Animation

The final animations are saved in subdirectory `png`. Note that they will be overwritten everytime you run the program.
