#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"


void sweep2d(double a[][maxn], double f[][maxn], int nx,
	     int s[2], int e[2], double b[][maxn])
{
  double h;
  int i,j;

  h = 1.0/((double)(nx+1));

  for (i = s[0]; i <= e[0]; i++) {
    for (j = s[1] ; j <= e[1] ; j++) {
      b[i][j] = 0.25 * ( a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j] );
    }
  }
}

void exchangrma_2d(double x[][maxn], int nx, int s[2], int e[2], MPI_Win win,
                 int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Comm comm) {
  /* offset: avoid left ghost col and boundary condition value */
  int num_vertical_elements = e[1] - s[1] + 1, num_horz_elements = e[0] - s[0] + 1; //This is the number of points shared across the boundary 

  MPI_Datatype row_type;
  MPI_Type_vector(num_horz_elements, 1, maxn, MPI_DOUBLE, &row_type); // We have to skip maxn as this is actually the size of the grid, not nx + 2.
  MPI_Type_commit(&row_type);

  //between two processors.
  MPI_Aint offset;
  MPI_Win_fence(MPI_MODE_NOPRECEDE, win);//We can do all of these within one win fence as all act on different parts of memory.
  // GET  from the right and above the ghost columns as we don't know the correct offsets to use puts.

  offset = maxn + s[1];//same offset as 1d case
  MPI_Get(&x[e[0]+1][s[1]], num_vertical_elements, MPI_DOUBLE, nbrright, offset, num_vertical_elements, MPI_DOUBLE, win);

  offset = maxn + e[1] + 1;//We want to move 1 column in hence an offset of max n. We then want to get from e[1] + 1 up the column.
  MPI_Get(&x[s[0]][e[1]+1], 1, row_type, nbrup, offset, 1, row_type, win);
  
  //Use MPI_Put to get the ghost columns for the left and the down processors as we know the correct address to put the data into the ghost. 
  offset = s[1];
  MPI_Put(&x[e[0]][s[1]], num_vertical_elements, MPI_DOUBLE, nbrright, offset, num_vertical_elements, MPI_DOUBLE, win);

  offset = maxn + e[1];//We want to move 1 column in hence an offset of max n. We then want to place this e[1] up the column.
  MPI_Put(&x[s[0]][e[1]], 1, row_type, nbrup, offset, 1, row_type, win);//we make use of row_type and only pass 1 datatype.

  MPI_Win_fence(MPI_MODE_NOSUCCEED, win);
}



void exchangrma_2d_pscw(double x[][maxn], int nx, int s[2], int e[2], MPI_Win win,
                 int nbrleft, int nbrright, int nbrup, int nbrdown, int myid, MPI_Group edits[2], MPI_Group edited_by[2]) {
  int num_vertical_elements = e[1] - s[1] + 1, num_horz_elements = e[0] - s[0] + 1; //This is the number of points shared across the boundary 

  MPI_Datatype row_type;
  MPI_Type_vector(num_horz_elements, 1, maxn, MPI_DOUBLE, &row_type); // We have to skip maxn as this is actually the size of the grid, not nx + 2.
  MPI_Type_commit(&row_type);
  MPI_Aint offset;

  if (nbrdown != MPI_PROC_NULL) {//This means there is a process below  myid which will access the local memory on myid
  	 MPI_Win_post(edited_by[1], 0, win);//This group is the process nbrdown which edits the window on myid
     MPI_Win_wait(win);
  }
 
  if (nbrup != MPI_PROC_NULL) {//This means it only excetutes if there is no process above myid which it needs to access the memory of
     MPI_Win_start(edits[1], 0, win);//these are the processes that myid will edit
     offset = maxn + e[1] + 1;//We want to move 1 column in hence an offset of max n. We then want to get from e[1] + 1 up the column.
     MPI_Get(&x[s[0]][e[1]+1], 1, row_type, nbrup, offset, 1, row_type, win);

     offset = maxn + e[1];//We want to move 1 column in hence an offset of max n. We then want to place this e[1] up the column.
     MPI_Put(&x[s[0]][e[1]], 1, row_type, nbrup, offset, 1, row_type, win);//we make use of row_type and only pass 1 datatype.
     MPI_Win_complete(win);
  }

  if (nbrleft != MPI_PROC_NULL) {//This means it won't execute unless there is a process to myid left which will access the local memory on myid
  	 MPI_Win_post(edited_by[0], 0, win);//This group is the processes, down and left each edit the window on myid
     MPI_Win_wait(win);
  }

  if (nbrright != MPI_PROC_NULL) {//This means it won't execute unless there is no process to myid right which it needs to access the memory of
     MPI_Win_start(edits[0], 0, win);//these are the processes that myid will edit
     offset = maxn + s[1];//same offset as 1d case
     MPI_Get(&x[e[0]+1][s[1]], num_vertical_elements, MPI_DOUBLE, nbrright, offset, num_vertical_elements, MPI_DOUBLE, win);

     offset = s[1];
     MPI_Put(&x[e[0]][s[1]], num_vertical_elements, MPI_DOUBLE, nbrright, offset, num_vertical_elements, MPI_DOUBLE, win);
     MPI_Win_complete(win);
  }
}


double griddiff_2d(double a[][maxn], double b[][maxn], int nx, int s[2], int e[2]) {
  double sum;
  double tmp;
  int i, j;

  sum = 0.0;

  for (i = s[0]; i <= e[0]; i++) {
    for ( j = s[1]; j <= e[1]; j++) {
      tmp = (a[i][j] - b[i][j]);
      sum = sum + sqrt(tmp*tmp);
    }
  }

  return sum;

}
