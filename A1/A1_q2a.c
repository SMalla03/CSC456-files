#include "mpi.h"
#include <stdio.h>
int main(argc, argv)
int argc;
char **argv;
{
    int rank, size;
    int num_proc = 32;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int a = rank + 1;
    int t1 = 0;
    for(int i = 0; i < num_proc; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        int result;
        
        double time1 = MPI_Wtime();
        MPI_Allreduce(&a, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        double time2 = MPI_Wtime();
        double dt = time2 - time1;
        double maxtime;
        MPI_Reduce(&dt, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        t1 += maxtime;
    }
    t1 = t1 / num_proc;
    int t2 = 0;
    for (int it = 0; it < nit; it++) { 
        MPI_Barrier(MPI_COMM_WORLD); 
        double t_start = MPI_Wtime(); 
        int a_2 = rank + 1;
        int b; 
        for (int i = 0; i < num_proc; i++) { 
            int partner = rank ^ (1 << i); 
            MPI_Sendrecv(&a_2, 1, MPI_INT, partner, 0, &b, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            a_2 += b; 
        } 
        double t_end = MPI_Wtime(); 
        double dt = t_end - t_start; 
        double dt_max; 
        MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
        t2 += dt_max;
    }
    t2 = t2 / num_proc;

    MPI_Finalize();
    if (rank == 0) {
        printf("Average time for MPI_Allreduce: %d\n", t1);
        printf("Average time for global sum up: %d\n", t2);
    }
    return 0;
}

