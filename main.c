#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <mpi.h>

int main(int argc, char **argv) {
    
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    MPI_Finalize();

    return 0;
}
