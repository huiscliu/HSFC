
#ifndef PHG_PARTITION_SFC_H

#include <mpi.h>

#if 1
typedef long int SFC_INT;
#define MPI_SFC_INT MPI_LONG
#else
typedef int SFC_INT;
#define MPI_SFC_INT MPI_INT
#endif

/* float point number used by SFC */
typedef double SFC_FLOAT;

#define MPI_SFC_FLOAT MPI_DOUBLE
#define HSFC_EPSILON  1.6e-7

#define SFC_TRUE  (1)
#define SFC_FALSE (0)

typedef struct {
    SFC_FLOAT low;
    SFC_FLOAT high;
    int pn;

} SFC_PTNS;

typedef struct {
    SFC_FLOAT key;   /* key \in [0, 1]      */
    SFC_FLOAT w;     /* weight              */
    int pn;          /* partition number    */

} SFC_DOTS;

/* Struct for a 3-D inverse SFC */
typedef struct {
    SFC_FLOAT co[3];	/* coordinates, input */
    SFC_FLOAT sfc;  	/* key (output, in the range (0,1)) */

    SFC_FLOAT weight;   /* weight, input */
    int part;           /* partition, output, returned by algorithm */

} SFC_ELEM;

#ifdef __cplusplus
extern "C" {
#endif

/* Hilbert space-filling curve partitioner
 *
 * nleaf should be greater than 0
 *
 * save_ptn is optional. If save_ptn is non-NULL, part should be non-NULL.
 * The length of save_ptn and part is the size of newcomm.
 * save_ptn stores the partion intervals, which is partitioned as [low, high)
 *
 * part[i] is the rank of i-th interval.
 *
 * */

double phgPartitionSFC(SFC_ELEM *x, SFC_INT nleaf, MPI_Comm oldcomm, MPI_Comm newcomm, double lif, int remap,
        SFC_PTNS **save_ptn, int **part);

/*------------------------- Internal Functions ------------------------*/
/* n <= length of array hsfc             */
/* in general, n == length of array hsfc */
/* inverse Hilbert SFC                   */
void phgSFCInvHilbert3D(SFC_ELEM *x, SFC_INT n);

double phgPartition1DV1(SFC_DOTS *dts, SFC_INT lx, int p, MPI_Comm comm, double lif, int *itr, SFC_PTNS **save_ptn);

int phgPartitionRemap(MPI_Comm comm, const double datasize[], int perm[], MPI_Comm newcomm);

/* search which interval includes key.
 * p is the number of intervals. ptn is a partition of (0, 1) */
int phgPartitionSearch(SFC_PTNS *ptn, int p, SFC_FLOAT key);

#ifdef __cplusplus
}
#endif

#define PHG_PARTITION_SFC_H
#endif
