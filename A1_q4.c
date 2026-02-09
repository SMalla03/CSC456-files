#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main(argc, argv)
int argc;
char **argv;
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int i = 1; i < 6; i++) {
        int n = 1000 * 1 << (i-1);
        int *v = malloc(n*n * sizeof(int));
        MPI_Barrier(MPI_COMM_WORLD);
        double t_start = MPI_Wtime();
        for (int j = 1; j <= (n-1)*n+1; j+=n) {
            MPI_Sendrecv(&v[j], 1, MPI_INT, (rank+1)%size, 0, &v[j], 1, MPI_INT, (rank-1+size)%size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        double t_end = MPI_Wtime();
        double dt = t_end - t_start;
        double dt_max;
        MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("p = %d,Average time in method 1 for n = %d: %.9f\n", size, n, dt_max);
        }
        free(v);

        int *v2 = malloc(n*n * sizeof(int));
        int *v2_short = malloc(n * sizeof(int));
        MPI_Barrier(MPI_COMM_WORLD);
        double t_start2 = MPI_Wtime();
        for (int j = 0; j < n; j++) {
            v2_short[j] = v2[j*n + 1];
        }
        MPI_Sendrecv(v2_short, n, MPI_INT, (rank+1)%size, 0, v2_short, n, MPI_INT, (rank-1+size)%size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int j = 0; j < n; j++) {
            v2[j*n + 1] = v2_short[j];
        }
        double t_end2 = MPI_Wtime();
        double dt2 = t_end2 - t_start2;
        double dt2_max;
        MPI_Allreduce(&dt2, &dt2_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("p = %d,Average time in method 2 for n = %d: %.9f\n", size, n, dt2_max);
        }
        free(v2);
        free(v2_short);
    }
    MPI_Finalize();
    return 0;
}


