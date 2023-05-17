# Kármán vortex street simulation using Lattice Boltzmann Method

Final project for NYU Graduate course MATH-GA 2012 Advanced Topics in Numerical Analysis: High Performance Computing.

Collaborative work of Yifei Zhu, Tiffany Li, Kitty Li.

* Source code is in the `d2q9-bgk.c` file
* Results checking scripts are in the `check/` directory

## Usage

To compile type `make`. Input parameter and obstacle files are all specified on the command line of the `d2q9-bgk` executable.

Usage:

    $ ./d2q9-bgk <paramfile> <obstaclefile>
eg:

    $ ./d2q9-bgk input_256x128.params obstacles_256x128.dat

## Running on BlueCrystal Phase 4

When you wish to submit a job to the queuing system on BlueCrystal, you should use the job submission script provided.

    $ sbatch job_submit_d2q9-bgk

This will dispatch a job to the queue, which you can monitor using the
`squeue` command:

    $ squeue -u $USER

When finished, the output from your job will be in a file called
`d2q9-bgk.out`:

    $ less d2q9-bgk.out

If you wish to run a different set of input parameters, you should
modify `job_submit_d2q9-bgk` to update the value assigned to `options`.

## Checking submission content

Before handing in the coursework, you can use the `check_submission.sh` script to make sure that your code builds in a clean environment. This will reduce the chances of the automarker failing to build or run your code.

To use the script, simply run it from the directory containing the files you intend to submit:

    $ /path/to/check_submission.sh

The script will:

1. Unload all the modules currently loaded.
2. Load your modules and environment variables specified in `env.sh`.
3. Use `make` to build your code and verify that an executable with the expected name is produced.

If the submission checking script prints any errors, you should try to address those before you hand in. 

Note that `check_submission.sh` does _not_ run your code, and so you _cannot_ verify that the results produced by your application validate just by running this script. You should check the correctness of your results separately, e.g. using `make check`.


# Serial output for sample inputs
Run times were taken on a Phase 4 node using the base (gcc) compiler and base compiler flags as found in the Makefile:

- 128x128
```
./d2q9-bgk  input_128x128.params obstacles_128x128.dat
==done==
Reynolds number:		9.751927375793E+00
Elapsed time:			38.387577 (s)
Elapsed user CPU time:		38.388736 (s)
Elapsed system CPU time:	0.003000 (s)
```

- 128x256
```
./d2q9-bgk  input_128x256.params obstacles_128x256.dat
==done==
Reynolds number:		3.715003967285E+01
Elapsed time:			77.446019 (s)
Elapsed user CPU time:		77.450619 (s)
Elapsed system CPU time:	0.003000 (s)
```

- 256x256
```
./d2q9-bgk  input_256x256.params obstacles_256x256.dat
==done==
Reynolds number:		1.005141162872E+01
Elapsed time:			309.040200 (s)
Elapsed user CPU time:		309.061111 (s)
Elapsed system CPU time:	0.004000 (s)
```

- 1024x1024
```
./d2q9-bgk  input_1024x1024.params obstacles_1024x1024.dat
==done==
Reynolds number:		3.375851392746E+00
Elapsed time:			1287.501875 (s)
Elapsed user CPU time:		1287.568113 (s)
Elapsed system CPU time:	0.029001 (s)
```

# Visualisation

You can view the final state of the simulation by creating a .png image file using a provided Gnuplot script:

    $ gnuplot final_state.plt
