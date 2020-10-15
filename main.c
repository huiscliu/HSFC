
#include "partition-sfc.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char **argv)
{
    int i, n, n_total = 0, n_max = 0, rank, nprocs;
    SFC_ELEM *x;
    MPI_Comm comm = MPI_COMM_WORLD;
    double real_lif;
    size_t seed;
    SFC_PTNS *save_ptn = NULL;
    int *part = NULL;
    double lif = 1.02;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    seed = time(0) * rank;
    srand(seed);
    srand(rand() + rank);

    n = 10 + rand() % 300000;
    MPI_Allreduce(&n, &n_total, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&n, &n_max, 1, MPI_INT, MPI_MAX, comm);

    printf("elements in rank %2d: %d\n", rank, n);

    x = malloc(n * sizeof(*x));
    for (i = 0; i < n; i++) {
        /* coordinate */
        x[i].co[0] = rand() * 1. * 10 / RAND_MAX + 1;
        x[i].co[1] = rand() * 1. * 10 / RAND_MAX + 3;
        x[i].co[2] = rand() * 1. * 10 / RAND_MAX + 5;

        /* weight */
        x[i].weight = 1.;
    }

    real_lif = phgPartitionSFC(x, n, comm, comm, lif, SFC_TRUE, &save_ptn, &part);
    if (rank == 0) {
        printf("\ntotal elements: %d, max elements: %d\n", n_total, n_max);
        printf("input load balance factor tolerance: %f\n\n", n_max * nprocs * 1. / n_total);
        printf("load balance factor tolerance: %f\n", lif);

        printf("real load balance factor after partition: %f\n\n", real_lif);

        if (save_ptn != NULL) {
            for (i = 0; i < nprocs; i++) {
                printf("interval %d, [%g, %g], rank: %d\n", i, save_ptn[i].low, save_ptn[i].high, part[i]);
            }
        }
    }

    MPI_Finalize();

    free(x);
    free(save_ptn);
    free(part);

    return 0;
}
