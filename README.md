# MPI-numerical-solution-to-the-poisson-equation

The code along with a makefile can be found within the Tarball file. The codes submitted in the
assignment are rma_blocking_2d_decomp.c, which contains the code for parts 1, rma_pswc_2d_decomp.c,
which the code for part 2, that is the pscw procedure.
These codes call jacobi.c, which contains the rma exchange functions 2d while decomp1d.c
which contains the 1d decompositon function.

On chuck, the MPI modules must be loaded with the command 'module load cports openmpi'

Run the makefile with 'make'

This will generate the execuatables rma_blocking_2d_decomp, rma_pscw_2d_decomp,
which are the execuatbles for the corresponding c files.

Each executable can be run with 'mpirun -np <number of processors> <executable name> <size of grid>

This will print the analytic solution for that grid size and the numerical grid which has been
gather onto the rank zero process and print from there

Textfiles called 'rma_blocking_2D.txt' and 'rma_pscw_2D.txt' will also be generated. These can be compared with the textfiles from assignment 2 by 
running the command 'diff sendrecv_2D_9x9.txt rma_pscw_2D.txt'. When the textfile 'rma_pscw_2D.txt' has been generated with the same number of grid point as the 
file 'sendrecv_2F_9x9.txt', in this case 9, there should be no difference between the two files and the solution should be the same.

