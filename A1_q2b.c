#include "mpi.h"
#include <stdio.h>

int main(argc, argv)
int argc;
char **argv;
{
    int rank, size;
    int nit = 50;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int d = 0; while ((1 << d) < size) d++;
    int a = rank + 1;
    double t1 = 0;
    for (int i = 0; i < nit; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        double t_start = MPI_Wtime(); 
        int result_sum = a;
        int result_prod = a;
        int max_val = a;
        int min_val = a;
        
        double time1 = MPI_Wtime();
        MPI_Allreduce(&a, &result_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&a, &result_prod, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
        MPI_Allreduce(&a, &max_val, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&a, &min_val, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        double time2 = MPI_Wtime();

        double dt = time2 - time1;
        double maxtime;
        MPI_Allreduce(&dt, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        t1 += maxtime;
    }
    t1 = t1 / nit;
    double t2 = 0;
    for (int i = 0; i < nit; i++) { 
        MPI_Barrier(MPI_COMM_WORLD); 
        double t_start = MPI_Wtime(); 
        int result_sum = a;
        int result_prod = a;
        int max_val = a;
        int min_val = a;
        int b; 
        for (int j = 0; j < d; j++) { 
            int partner = rank ^ (1 << j); 
           // sum
            MPI_Sendrecv(&result_sum, 1, MPI_INT, partner, 0, &b, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            result_sum += b; 
            // product 
            MPI_Sendrecv(&result_prod, 1, MPI_INT, partner, 0, &b, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            result_prod *= b; 
            // max/min 
            MPI_Sendrecv(&max_val, 1, MPI_INT, partner, 0, &b, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (b > max_val) max_val = b; 
            if (b < min_val) min_val = b;
        }
        double t_end = MPI_Wtime(); 
        double dt = t_end - t_start; 
        double dt_max; 
        MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
        t2 += dt_max;
    }
    t2 = t2 / nit;

    if (rank == 0) {
        printf("Average time for MPI allreduce, prod, max, min: %.9f\n", t1);
        printf("Average time for global sum up, prod, max, min: %.9f\n", t2);
    }
    MPI_Finalize();
    return 0;
}


