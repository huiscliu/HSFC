
#include "partition-sfc.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char **argv)
{
    int i, n, gs = 0, rank, nprocs;
    SFC_ELEM *x;
    MPI_Comm comm = MPI_COMM_WORLD;
    double lif;
    size_t seed;
    SFC_PTNS *save_ptn = NULL;
    int *part = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    seed = time(0) * rank;
    srand(seed);
    srand(rand() + rank);

    n = 10 + rand() % 300000;
    MPI_Scan(&n, &gs, 1, MPI_INT, MPI_SUM, comm);
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

    lif = phgPartitionSFC(x, n, comm, comm, SFC_TRUE, &save_ptn, &part);
    if (rank == 0) {
        printf("\nload balance factor: %f\n\n", lif);

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
