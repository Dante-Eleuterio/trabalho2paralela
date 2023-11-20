//////////////////////////////////////////////////////////////////////////
// Traveling Salesman Problem with MPI
// Note: this is a C++ program.
//
// Process 0 manages the queue of incomplte paths, and keeps
// track of the best path.
//
// All other processes are workers that get and put jobs from and into 
// the queue. Each time they get a path, they are also informed
// about the best length so far.
//
// Note that, unlike previous examples, this one does not work with
// only one process.
//
// Starting city is assumed to be city 0.
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <mpi.h>
#include <memory.h>
int myrank, NumProcs, NumCities;
void Worker(){
    printf("COMUNISMO PORRA\n");
}

void Coordinator(){
    printf("EU SOU UM PRODUTOR\n");
}
void Fill_Dist(){
    printf("TESTE\n");
}

int main(int argc, char *argv[])
{
  Fill_Dist();  // process 0 read the data and broadcast it to the others
  MPI_Init (&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  printf("%d\n",myrank);

  if (NumProcs<2) {
    printf("At least 2 processes are required\n");
    exit(-1);
  }  


  // Initialize distance matrix. Ususally done by one process 
  // and bcast, or initialized from a file in a shared file system.

  if (myrank==0) 
    Coordinator();
  else
    Worker();
  
  MPI_Finalize();
  return 0;
}
