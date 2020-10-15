
#include "partition-sfc.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

static int phgVerbosity = 0;

/* max loops */
static int MAXLOOPS = 20;

static const SFC_FLOAT SFC_LBE = 1.004;

static SFC_INT CUTS = 16;

typedef SFC_FLOAT HSFC_FARY[2];

/* **************************************************************
 * Inverse Hilbert Space-Filling curve (3d)
 * converte coordinates in [0,1]x[0,1]x[0,1] to [0, 1]
 * n <= the length of array x 
 * 
 * Here we give two versions, the first one is general method
 * that based on table, which is based on Zoltan's.
 * We extented original implement, and we fixed some bugs.
 * The second one is based on an improved algorithm, which we
 * modified the first one. The time complexity is O(r), 
 * r = log2(max(x, y, z)) + 1. In general, r < n;
 * Both are fast. User can choose any one.
 * If grid has big coordinates much more than small coordinates, the 
 *   first maybe better.
 * If grid has small coordinats much more than bigs ones, the second
 *   maybe better.
 *
 * assume sizeof(SFC_INT) >= 4
 * **************************************************************/

/* data: idata2d, idata3d, istat2d, istat3d
 * are borrowed from zoltan.
 * They are used for generating inverse Hilbert SFC */

/* 3 dimension to nkey conversion */
static unsigned const int idata3d [] = {
    0,  7,  3,  4,  1,  6,  2,  5,
    0,  1,  3,  2,  7,  6,  4,  5,
    0,  3,  7,  4,  1,  2,  6,  5,
    2,  3,  5,  4,  1,  0,  6,  7,
    4,  5,  3,  2,  7,  6,  0,  1,
    4,  7,  3,  0,  5,  6,  2,  1,
    6,  7,  5,  4,  1,  0,  2,  3,
    0,  1,  7,  6,  3,  2,  4,  5,
    2,  1,  5,  6,  3,  0,  4,  7,
    6,  1,  5,  2,  7,  0,  4,  3,
    0,  7,  1,  6,  3,  4,  2,  5,
    2,  1,  3,  0,  5,  6,  4,  7,
    4,  7,  5,  6,  3,  0,  2,  1,
    4,  5,  7,  6,  3,  2,  0,  1,
    6,  1,  7,  0,  5,  2,  4,  3,
    0,  3,  1,  2,  7,  4,  6,  5,
    2,  3,  1,  0,  5,  4,  6,  7,
    6,  7,  1,  0,  5,  4,  2,  3,
    2,  5,  1,  6,  3,  4,  0,  7,
    4,  3,  7,  0,  5,  2,  6,  1,
    4,  3,  5,  2,  7,  0,  6,  1,
    6,  5,  1,  2,  7,  4,  0,  3,
    2,  5,  3,  4,  1,  6,  0,  7,
    6,  5,  7,  4,  1,  2,  0,  3
};


/* 3 dimension to nkey state transitions */
static unsigned const int istate3d [] ={
    1,  6,  3,  4,  2,  5,  0,  0,
    0,  7,  8,  1,  9,  4,  5,  1,
    15, 22, 23, 20,  0,  2, 19,  2,
    3, 23,  3, 15,  6, 20, 16, 22,
    11,  4, 12,  4, 20,  1, 22, 13,
    22, 12, 20, 11,  5,  0,  5, 19,
    17,  0,  6, 21,  3,  9,  6,  2,
    10,  1, 14, 13, 11,  7, 12,  7,
    8,  9,  8, 18, 14, 12, 10, 11,
    21,  8,  9,  9,  1,  6, 17,  7,
    7, 17, 15, 12, 16, 13, 10, 10,
    11, 14,  9,  5, 11, 22,  0,  8,
    18,  5, 12, 10, 19,  8, 12, 20,
    8, 13, 19,  7,  5, 13, 18,  4,
    23, 11,  7, 17, 14, 14,  6,  1,
    2, 18, 10, 15, 21, 19, 20, 15,
    16, 21, 17, 19, 16,  2,  3, 18,
    6, 10, 16, 14, 17, 23, 17, 15,
    18, 18, 21,  8, 17,  7, 13, 16,
    3,  4, 13, 16, 19, 19,  2,  5,
    16, 13, 20, 20,  4,  3, 15, 12,
    9, 21, 18, 21, 15, 14, 23, 10,
    22, 22,  6,  1, 23, 11,  4,  3,
    14, 23,  2,  9, 22, 23, 21,  0
};

/*  maximum level */
/*  in this file, we force that MAXLEVEL >= 20 && MAXLEVEL <= 30 */
#define  MAXLEVEL  25

void phgSFCInvHilbert3D(SFC_ELEM *x, SFC_INT n)
{
    static unsigned const int *d[] = {
        idata3d,         idata3d + 8,   idata3d + 16,  idata3d + 24,
        idata3d + 32,    idata3d + 40,  idata3d + 48,  idata3d + 56,
        idata3d + 64,    idata3d + 72,  idata3d + 80,  idata3d + 88,
        idata3d + 96,    idata3d + 104, idata3d + 112, idata3d + 120,
        idata3d + 128,   idata3d + 136, idata3d + 144, idata3d + 152,
        idata3d + 160,   idata3d + 168, idata3d + 176, idata3d + 184
    };

    static unsigned const int *s[] = {
        istate3d,       istate3d + 8,   istate3d + 16,  istate3d + 24,
        istate3d + 32,  istate3d + 40,  istate3d + 48,  istate3d + 56,
        istate3d + 64,  istate3d + 72,  istate3d + 80,  istate3d + 88,
        istate3d + 96,  istate3d + 104, istate3d + 112, istate3d + 120,
        istate3d + 128, istate3d + 136, istate3d + 144, istate3d + 152,
        istate3d + 160, istate3d + 168, istate3d + 176, istate3d + 184
    };

    int level, EL;
    unsigned int key[3], c[3], temp, state;
    SFC_INT i;

    unsigned IMAX;
    unsigned EfBit;
    int k0 = 0, k1 = 0, k2 = 0;

    assert(sizeof(int) >= 4);
    assert(sizeof(unsigned int) >= 4);

    /* Init IMAX such that all effective bits are 1 */
    /* IMAX = 2^32 - 1 */
    IMAX = 4294967295U;

    /* 30 effective bits, all 1 */
    /* EfBit = 2^30 - 1 */
    EfBit = IMAX >> 2;


    k0 = 60 - MAXLEVEL * 3;
    k1 = 30 - MAXLEVEL * 3;
    k2 = - MAXLEVEL * 3;

    for (i = 0; i < n; i++) {
        /* convert x,y,z coordinates to integers in range [0,IMAX] */
        c[0] = (unsigned int)(x[i].co[0] * (double)IMAX);     /* x */
        c[1] = (unsigned int)(x[i].co[1] * (double)IMAX);     /* y */
        c[2] = (unsigned int)(x[i].co[2] * (double)IMAX);     /* z */
        c[1] >>= 1;
        c[2] >>= 2;

        /* use state tables to convert nested quadrant's 
         * coordinates level by level */
        key[0] = key[1] = key[2] = 0;
        state = 0;
        EL = 30;
        for (level = 0; level < MAXLEVEL; level++) {
            /* extract 3 bits at current level */
            EL--;
            temp = ((c[0] >> EL) & 4)
                |((c[1] >> EL) & 2)
                |((c[2] >> EL) & 1);

            /* treat key[] as long shift register */
            /* shift in converted coordinate */
            /* every key[] has thirty effective bit */
            /* the last 30 are effective */
            key[0] = (key[0] << 3) | ((key[1] >> 27) & 7);
            key[1] = (key[1] << 3) | ((key[2] >> 27) & 7);
            key[2] = (key[2] << 3) | *(d[state] + temp);

            state = *(s[state] + temp);
        }

        key[0] = key[0] & EfBit;
        key[1] = key[1] & EfBit;
        key[2] = key[2] & EfBit;

        /* convert 3 part Hilbert key to double */
        x[i].sfc  = ldexp((double)key[2], k2);
        x[i].sfc += ldexp((double)key[1], k1);
        x[i].sfc += ldexp((double)key[0], k0);
    }
}

double phgPartition1DV1(SFC_DOTS *dts, SFC_INT nelem, int np, MPI_Comm comm, double lif, int *itr, SFC_PTNS **save_ptn)
{
    int npctss = (np - 1) * CUTS + 1;
    SFC_PTNS *ptn;

    SFC_FLOAT lbe;
    SFC_FLOAT dsum;
    SFC_FLOAT *ptsum;
    SFC_FLOAT *sum;
    SFC_FLOAT *tsum;
    HSFC_FARY *bbox;
    SFC_INT i, j;
    SFC_INT ctn;
    int rank;
    int flg;
    int *map;
    SFC_FLOAT *tmp;
    SFC_INT *ne;
    SFC_INT *npe;
    int *cts;
    SFC_FLOAT hsfc_lb_tol = SFC_LBE;

    assert((CUTS >= 8) && (CUTS <= 30));
    assert(nelem >= 0 && np >= 1);

    if (nelem > 0) {
        assert(dts != NULL);
    }

    if (np == 1) {
        for (i = 0; i < nelem; i++) {
            dts[i].pn = 0;
        }

        return 1.;
    }

    j = nelem;
    MPI_Allreduce(&j, &ctn, 1, MPI_SFC_INT, MPI_SUM, comm);
    if (ctn == 0) {
        return 1.;
    }

    lbe = 0.;
    for (i = 0; i < nelem; i++) {
        lbe += dts[i].w;
    }

    MPI_Allreduce(&lbe, &dsum, 1, MPI_SFC_FLOAT, MPI_SUM, comm);

    if (dsum == 0.) {
        for (i = 0; i < nelem; i++) dts[i].pn = 0;

        return 1.;
    }
    else if (dsum < 0.) {
        printf("negative weights.\n");
        MPI_Abort(comm, -1);
    }

    if (lif < 1 || lif > 2) {
        hsfc_lb_tol = SFC_LBE;
    }
    else {
        hsfc_lb_tol = lif;
    }

    ptsum = malloc(np * sizeof(*ptsum));
    bbox = malloc(np * sizeof(*bbox));
    lbe = 1.0 / np;
    ctn = np - 1;
    for (i = 0; i < ctn; i++) {
        ptsum[i] = dsum * (i + 1.) * lbe;
        bbox[i][0] = 0.;
        bbox[i][1] = 1.;
    }

    ptsum[ctn] = dsum;
    bbox[ctn][0] = 0.;
    bbox[ctn][1] = 1.;

    tsum = malloc(npctss * sizeof(*tsum));
    sum = malloc(npctss * sizeof(*sum));
    map = malloc(npctss * sizeof(*map));
    tmp = malloc(np * sizeof(*tmp));
    ne = malloc(npctss * sizeof(*ne));
    npe = malloc(npctss * sizeof(*npe));
    cts = malloc(npctss * sizeof(*cts));

    ptn = malloc(npctss * sizeof(*ptn));
    lbe = 1. / npctss;
    for (i = 0; i < npctss; i++) {
        ptn[i].low = i * lbe;
        ptn[i].high = (i + 1) * lbe;
        ptn[i].pn = i;
    }

    /* if save ptn */
    if (save_ptn != NULL) {
        *save_ptn = malloc(sizeof(**save_ptn) * np);
    }

    MPI_Comm_rank(comm, &rank);

    flg = 1;
    ctn = 1;
    while (flg <= MAXLOOPS && ctn) {
        SFC_FLOAT mxm;
        int pmo = np - 1;
        int kk = 0;

        for (i = 0; i < npctss; i++) {
            tsum[i] = 0.;
            npe[i] = 0;
        }

        for (i = 0; i < nelem; i++) {
            SFC_FLOAT key = dts[i].key;
            int lo, hi;
            int md;
            int ret;
            int nprocs = npctss;

            lo = 0;
            hi = nprocs - 1;
            md = (lo + hi) >> 1;
            while (SFC_TRUE) {
                if (key >= ptn[md].high)
                    ret = 1;
                else if (key - ptn[md].low < 0.)
                    ret = -1;
                else
                    ret = 0;

                if (ret == 1) {
                    lo = md;
                    md = (lo + hi) >> 1;
                    if (md == lo)
                        md = hi;
                }
                else if (ret == -1) {
                    hi = md;
                    md = (lo + hi) >> 1;
                }
                else if (hi == lo) {
                    break;
                }
                else {
                    break;
                }
            }

            dts[i].pn = md;
            tsum[md] += dts[i].w;
            npe[md] += 1;
        }

        MPI_Reduce(tsum, sum, npctss, MPI_SFC_FLOAT, MPI_SUM, 0, comm);
        MPI_Reduce(npe, ne, npctss, MPI_SFC_INT, MPI_SUM, 0, comm);

        if (rank == 0) {
            tsum[0] = sum[0];
            cts[0] = 0;
            for (i = 1; i < npctss; i++) {
                tsum[i] = sum[i];
                sum[i] += sum[i - 1];
                cts[i] = 0;
            }

            j = 0;
            for (i = 0; i < pmo; i++) {
                while (sum[j] < ptsum[i]) {
                    j++;
                }

                if (ptn[j].high < bbox[i][1])
                    bbox[i][1] = ptn[j].high;
                if (ptn[j].low > bbox[i][0])
                    bbox[i][0] = ptn[j].low;

                cts[j] += 1;

                tmp[i] = 0.;
            }
            tmp[np - 1] = 0.;

            j = 0;
            for (i = 0; i < npctss && j < np; i++) {
                if (sum[i] < ptsum[j]) {
                    map[i] = j;
                }
                else {
                    map[i] = j;
                    j++;

                    if ((ne[i] > 1) && (j < pmo)) {
                        kk = 1;
                    }
                }
            }

            for (j = i; j < npctss; j++) {
                map[j] = pmo;
            }

            for (i = 0; i < npctss; i++) {
                j = map[i];
                tmp[j] += tsum[i];
            }
            mxm = tmp[0];
            for (i = 1; i < np; i++) {
                if (mxm < tmp[i])
                    mxm = tmp[i];
            }
            lbe = mxm * np / dsum;

            if ((lbe <= hsfc_lb_tol) || (kk == 0 && np >= 3)) {
                ctn = 0;
            }
            else {
                ctn = 1;
            }
        }

        MPI_Bcast(map, npctss, MPI_INT, 0, comm);
        MPI_Bcast(&lbe, 1, MPI_SFC_FLOAT, 0, comm);
        MPI_Bcast(&ctn, 1, MPI_SFC_INT, 0, comm);

        /* save ptn */
        if (save_ptn != NULL) {
            SFC_PTNS *ps = *save_ptn;

            j = 0;
            ps[0].low = ptn[0].low;
            ps[0].high = ptn[0].high;
            for (i = 0; i < npctss; i++) {
                /* increase upper bound */
                if (j == map[i]) {
                    ps[j].high = ptn[i].high;
                }
                else {
                    j++;

                    /* record lower and upper bound */
                    ps[j].low = ptn[i].low;
                    ps[j].high = ptn[i].high;
                }
            }

            /* special case */
            if (j < np - 1) {
                j += 1;

                for (; j < np; j++) {
                    ps[j].low = ps[j].high = 1.;
                }
            }
        }

        if (ctn) {
            ctn = npctss - 1;

            if (rank == 0) {
                SFC_FLOAT len;
                SFC_INT k, jk1, j2;
                SFC_INT ctk;

                i = 0;
                j = 0;
                k = 0;
                while (i < pmo) {
                    while (cts[j] < 1) {
                        j++;
                    }

                    j2 = cts[j] * CUTS;
                    len = 1. / j2;
                    ctk = i + cts[j] - 1;

                    for (jk1 = 1; jk1 <= j2; jk1++) {
                        sum[k] = (bbox[i][0] * (j2 - jk1)
                                + bbox[ctk][1] * jk1) * len;
                        k++;
                    }

                    i += cts[j];
                    j++;
                }
            }

            MPI_Bcast(sum, ctn, MPI_SFC_FLOAT, 0, comm);

            ptn[0].low = 0.;
            ptn[0].high = sum[0];
            for (i = 1; i < ctn; i++) {
                ptn[i].high = sum[i];
                ptn[i].low = sum[i - 1];
            }
            ptn[ctn].high = 1.;
            ptn[ctn].low = sum[ctn - 1];

            ctn = 1;
        }
        else {
            break;
        }

        flg++;
    }

    for (i = 0; i < nelem; i++) {
        j = dts[i].pn;
        dts[i].pn = map[j];
    }

    if (phgVerbosity > 0 && rank == 0) {
        printf("1DV1, ITR: %d\n", flg);
        printf("1DV1, LBE: %f\n", lbe);
    }

    free(ptn);
    free(sum);
    free(ptsum);
    free(tsum);
    free(map);
    free(tmp);
    free(ne);
    free(npe);
    free(cts);
    free(bbox);

    if (itr != NULL) *itr = flg;

    return lbe;
}

/* Hilbert space-filling curve partitioner */
/* for 3d space, can be extended to 2d and 3d with slight change */
double phgPartitionSFC(SFC_ELEM *x, SFC_INT nleaf, MPI_Comm oldcomm, MPI_Comm newcomm, double lif, int remap,
        SFC_PTNS **save_ptn, int **part)
{
    int nprocs, rank;
    int onprocs;
    SFC_DOTS *dots;
    SFC_FLOAT out[6];    /* min and max coordinats */
    SFC_FLOAT ext[3];    /* scaling factor */
    SFC_FLOAT temp;
    SFC_FLOAT real_lif = 1.0;
    SFC_INT i;
    SFC_FLOAT t[3];
    int fail = SFC_FALSE;
    SFC_INT nleaf_global;
    SFC_INT idx_off = 0;

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &rank);
    MPI_Comm_size(oldcomm, &onprocs);

    /* check */
    assert(nleaf > 0);
    assert(x != NULL);

    /* total elements */
    MPI_Allreduce(&nleaf, &nleaf_global, 1, MPI_SFC_INT, MPI_SUM, oldcomm);

    if (nprocs == 1 && onprocs == 1) {
        for (i = 0; i < nleaf; i++) x[i].part = 0;
        return 1;
    }

    /* find bounding box */
    ext[0] = ext[1] = ext[2] = -1e300;

    for (i = 0; i < nleaf; i++) {
        if (x[i].co[0] > ext[0]) ext[0] = x[i].co[0];
        if (x[i].co[1] > ext[1]) ext[1] = x[i].co[1];
        if (x[i].co[2] > ext[2]) ext[2] = x[i].co[2];
    }

    MPI_Allreduce(ext, t, 3, MPI_SFC_FLOAT, MPI_MAX, oldcomm);

    out[3] = t[0];
    out[4] = t[1];
    out[5] = t[2];

    ext[0] = ext[1] = ext[2] = 1e300;

    for (i = 0; i < nleaf; i++) {
        if (x[i].co[0] < ext[0]) ext[0] = x[i].co[0];
        if (x[i].co[1] < ext[1]) ext[1] = x[i].co[1];
        if (x[i].co[2] < ext[2]) ext[2] = x[i].co[2];
    }

    MPI_Allreduce(ext, t, 3, MPI_SFC_FLOAT, MPI_MIN, oldcomm);
    out[0] = t[0];
    out[1] = t[1];
    out[2] = t[2];

    /* enlarge bounding box */
    temp = (out[3] - out[0]) * HSFC_EPSILON;
    out[3] += temp;
    out[0] -= temp;
    temp = (out[4] - out[1]) * HSFC_EPSILON;
    out[4] += temp;
    out[1] -= temp;
    temp = (out[5] - out[2]) * HSFC_EPSILON;
    out[5] += temp;
    out[2] -= temp;

    /* scaling coordinates, such that all coordinates belong to (0, 1)^3 */
    ext[0] = out[3] - out[0]; /* length of x axis of bounding box */
    ext[1] = out[4] - out[1]; /* length of y axis of bounding box */
    ext[2] = out[5] - out[2]; /* length of z axis of bounding box */

    /* preserve domain ratio */
    if (ext[0] < ext[1]) ext[0] = ext[1];
    if (ext[0] < ext[2]) ext[0] = ext[2];

    if (ext[0] == 0.) ext[0] = 1.;
    ext[0] = 1. / ext[0];
    ext[1] = ext[0];
    ext[2] = ext[0];

    for (i = 0; i < nleaf; i++) {
        x[i].co[0] = (x[i].co[0] - out[0]) * ext[0];
        x[i].co[1] = (x[i].co[1] - out[1]) * ext[1];
        x[i].co[2] = (x[i].co[2] - out[2]) * ext[2];
    }

    /* generate inverse SFC */
    phgSFCInvHilbert3D(x, nleaf);

    /* initialize dots */
    dots = (SFC_DOTS *)malloc(nleaf * sizeof(SFC_DOTS));
    if (nleaf > 0 && dots == NULL) {
        printf("memory allocation error (%s:%d)\n", __FILE__, __LINE__);
        MPI_Abort(oldcomm, -1);
        return 1;
    }

    for (i = 0; i < nleaf; i++) {
        dots[i].key = x[i].sfc;
    }

    /* weighted condition */
    for (i = 0; i < nleaf; i++) {
        assert(x[i].weight >= 0);

        dots[i].w = x[i].weight;
    }

    /* offset */
    MPI_Scan(&nleaf, &idx_off, 1, MPI_INT, MPI_SUM, oldcomm);
    idx_off -= nleaf;

    /* 1d  partition */
    if (nleaf_global <= nprocs) {
        SFC_INT gid;

        for (i = 0; i < nleaf; i++) {
            gid = idx_off + i;
            x[i].part = gid % nprocs;
        }

        real_lif = nprocs * 1. / nleaf_global;
    }
    else {
        int itr;

        real_lif = phgPartition1DV1(dots, nleaf, nprocs, oldcomm, lif, &itr, save_ptn);

        /* assign */
        for (i = 0; i < nleaf; i++) {
            x[i].part = dots[i].pn;
        }

        if (itr >= MAXLOOPS && real_lif >= 1.5) fail = SFC_TRUE;
    }

    free(dots);

    /* if 1d partitioners fail,
     * partition the grid according to its index */
    if (fail) {
        SFC_INT kk = nleaf_global / nprocs;
        SFC_INT rr = nleaf_global % nprocs;
        SFC_INT j;

        if (rank == 0) {
            printf("SFC - 1D partitioner failed.\n");
            printf("SFC - grid partitioned according to element index.\n");
        }

        j = (kk + 1) * rr;
        for (i = 0; i < nleaf; i++) {
            SFC_INT gid = idx_off + i;
            if (gid <= j) {
                x[i].part = gid / (kk + 1);
            }
            else {
                x[i].part = rr + (gid - j) / kk;
            }
        }

        if (rr == 0) {
            real_lif = 1.;
        }
        else {
            real_lif = (kk + 1) * nprocs * 1. / nleaf_global;
        }
    }

    /* remap */
    if (remap) {
        double *count;
        int *perm;
        int ret;

        count = calloc(nprocs, sizeof(*count));
        perm = calloc(nprocs, sizeof(*perm));

        for (i = 0; i < nleaf; i++) {
            count[x[i].part] += 1;
        }

        /* remap */
        ret = phgPartitionRemap(oldcomm, count, perm, newcomm);

        /* if changed */
        if (ret > 0) {
            for (i = 0; i < nleaf; i++) {
                x[i].part = perm[x[i].part];
            }
        }

        free(count);

        if (save_ptn != NULL) {
            assert(part != NULL);
            *part = perm;
        }
        else {
            free(perm);
        }
    }
    else {
        /* if save partion */
        if (save_ptn != NULL) {
            assert(part != NULL);

            *part = calloc(nprocs, sizeof(**part));
            for (i = 0; i < nprocs; i++) (*part)[i] = i;
        }

    }

    return real_lif;
}

/* phgPartitionRemap remaps the partitions such that the
 * the mount of migrated data will be minimized.
 * newcomm: the new communicator (new number of partitions)
 * return value > 0, if partitions are changed;
 * return value = 0, if partitions arn't changed;
 * return value < 0, if subroutines fail and then partitions keep the same.
 *
 * The algorithm comes from PLUM: 
 * parallel load balancing for adaptive unstructured meshed 
 * Leonid Oliker, Rupak Biswas 
 *
 * heuristic algorithm */

/* similiar matrix */
typedef struct {
    double s;	/* number of elements */
    int pid;	/* current pid */
    int ptn;	/* partition id */
} SM;

/* descending order used by qsort */
    static int
comp_sml(const void *m1, const void *m2)
{
    SM *n1, *n2;
    int flag = 0;

    n1 = (SM *)m1;
    n2 = (SM *)m2;

    if (n1->s > n2->s)
        flag = -1;
    else if (n1->s == n2->s)
        flag = 0;
    else if (n1->s < n2->s)
        flag = 1;
    return flag;
}

/* Note: either ('comm' \subset 'newcomm') or ('newcomm' \subset 'comm'),
 * and this function can be called by all processes of either 'comm' or the
 * larger of 'comm' and 'newcomm' */
int phgPartitionRemap(MPI_Comm comm, const double datasize[], int perm[], MPI_Comm newcomm)
{
    SM *local, *sml = NULL;
    MPI_Datatype type;
    MPI_Datatype types[3];
    MPI_Aint indices[3];
    SFC_INT i;
    SFC_INT nt = 0;
    int rank, size, nprocs;
    int nonz;
    int blocklens[3];
    int *len;
    int *ofst;
    int ret = 0;

    MPI_Comm_size(comm, &size);

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &rank);

    if (rank >= nprocs) return ret;

    assert((nprocs >= 1) && (datasize != NULL) && (perm != NULL));
    if (nprocs == 1) {
        perm[0] = 0;
        return ret;
    }

#if DEBUG
    /* communication datasize should not less than 0 */
    for (i = 0; i < nprocs; i++) {
        assert(datasize[i] >= 0.);
    }
#endif	/* DEBUG */

    local = (SM *) malloc(nprocs * sizeof(SM));

    if (local == NULL) {
        printf("memory allocation error (%s:%d)\n", __FILE__, __LINE__);
        MPI_Abort(comm, -1);
    }

    /* create a new datatype for exchange information */
    /* used by structure tree */
    blocklens[0] = blocklens[1] = blocklens[2] = 1;
    types[0] = MPI_DOUBLE;
    types[1] = types[2] = MPI_INT;
#if 1	/* MPI-2 */
    MPI_Get_address(&local[0].s, &indices[0]);
    MPI_Get_address(&local[0].pid, &indices[1]);
    MPI_Get_address(&local[0].ptn, &indices[2]);
#else	/* MPI-1 */
    MPI_Address(&local[0].s, &indices[0]);
    MPI_Address(&local[0].pid, &indices[1]);
    MPI_Address(&local[0].ptn, &indices[2]);
#endif
    indices[2] -= indices[0];
    indices[1] -= indices[0];
    indices[0] = 0;
#if 1	/* MPI-2 */
    MPI_Type_create_struct(3, blocklens, indices, types, &type);
#else	/* MPI-1 */
    MPI_Type_struct(3, blocklens, indices, types, &type);
#endif
    MPI_Type_commit(&type);

    /* assemble local similiar matrix */
    for (i = 0; i < nprocs; i++) {
        local[i].s = datasize[i];
        local[i].pid = rank;
        local[i].ptn = i;
    }

    /* calculate non-zero entries, whose rank is less than nprocs.
     * and move non-zero entries to the front part of the array */
    nonz = 0;
    if (rank < nprocs) {
        for (i = 0; i < nprocs; i++) {
            if (local[i].s > 0) {
                local[nonz].s = local[i].s;
                local[nonz].pid = local[i].pid;
                local[nonz].ptn = local[i].ptn;
                nonz++;
            }
        }
    }

    /* process 0 gathers the number of non-zero entries from 
     * all processes and calculate offset  */
    len = (int *)malloc(size * sizeof(int));
    ofst = (int *)malloc(size * sizeof(int));

    if ((len == NULL) || (ofst == NULL)) {
        printf("memory allocation error (%s: %d)\n", __FILE__, __LINE__);
        MPI_Abort(comm, -1);
    }

    MPI_Gather(&nonz, 1, MPI_INT, len, 1, MPI_INT, 0, comm);

    if (rank == 0) {
        ofst[0] = 0;
        for (i = 1; i < size; i++) {
            ofst[i] = ofst[i - 1] + len[i - 1];
        }

        /* allocate space for process 0 */
        nt = ofst[size - 1] + len[size - 1];
        sml = (SM *) malloc(nt * sizeof(SM));
        if (nt > 0 && sml == NULL) {
            printf("memory allocation error (%s:%d)\n", __FILE__,__LINE__);
            MPI_Abort(comm, -1);
        }
    }

    /* process 0 build similiar matrix */
    MPI_Gatherv(local, nonz, type, sml, len, ofst, type, 0, comm);

    MPI_Type_free(&type);
    free(local);
    free(len);
    free(ofst);

    /* map algorithm by PLUM */
    if (rank == 0){
        int *part_map;
        int *proc_unmap;
        int *t1;
        int count;
        int kk;
        SM *p;

        part_map = (int *) malloc(nprocs * sizeof(int));
        proc_unmap = (int *) malloc(nprocs * sizeof(int));
        t1 = (int *) malloc(nprocs * sizeof(int));

        if ((part_map == NULL) || (proc_unmap == NULL) || (t1 == NULL)) {
            printf("memory allocation error (%s:%d)\n", __FILE__,__LINE__);
            MPI_Abort(comm, -1);
        }

        /* initialize */
        for (i = 0; i < nprocs; i++) {
            proc_unmap[i] = 1;
            part_map[i] = 0;
            perm[i] = -1;
        }
        /* order the similiar matrix, descending order */
        qsort(sml, nt, sizeof(SM), comp_sml);

        /* map algorithm */
        count = 0;
        kk = 0;
        p = sml;
        nt -= 1;
        while ((count < nprocs) && (kk < nt)) {
            while ( (proc_unmap[p->pid] == 0 || part_map[p->ptn] == 1)
                    && (kk < nt)) {
                p++;
                kk++;
            }
            if ((proc_unmap[p->pid] != 0) && (part_map[p->ptn] != 1)) {
                proc_unmap[p->pid] = 0;
                part_map[p->ptn] = 1;
                perm[p->ptn] = p->pid;
                count++;
            }
        } /* map algorithm */

        count = 0;
        for (i = 0; i < nprocs; i++) {
            if (perm[i] < 0)
                count++;
            t1[i] = 0;
        }

        /* deal with special case that some arn't mapped */
        if (count > 0) {
            for (i = 0; i < nprocs; i++) {
                if (perm[i] >= 0) {
                    t1[perm[i]] = -1;
                }
            }
            kk = 0;
            for (i = 0; i < nprocs; i++) {
                if (perm[i] < 0) {
                    while (t1[kk] < 0) {
                        kk++;
                    }
                    perm[i] = kk;
                    kk++;
                }
            }
        }

        /* verification */
        /* count = 0: valid result   */
        /* count > 0: invalid result. unmap. */
        count = 0;
        for (i = 0; i < nprocs; i++) {
            t1[i] = 0;

            if ((perm[i] < 0) || (perm[i] > nprocs)) {
                count++;
            }
        }

        /* check if each index is mapped just once 
         * count = 0 if yes, > 0 if not */
        if (count == 0) {
            for (i = 0; i < nprocs; i++) {
                t1[perm[i]] += 1;
            }
            for (i = 0; i < nprocs; i++) {
                if (t1[i] != 1)
                    count++;
            }
        }

        if (count == 0) { /* succeed */
            /* default: partitions unchanged. */
            ret = 0;

            /* check if partitions are changed */
            for (i = 0; i < nprocs; i++) {
                if (perm[i] != i)
                    ret = 1;
            }
        }
        else { /* fail */
            ret = -1;
        }

        free(part_map);
        free(proc_unmap);
        free(t1);
    }

    /* broadcast results */
    MPI_Bcast(&ret, 1, MPI_INT, 0, comm);
    if (ret >= 0) {
        MPI_Bcast(perm, nprocs, MPI_INT, 0, comm);
    }
    else if (ret < 0){
        for (i = 0; i < nprocs; i++) {
            perm[i] = i;;
        }
    }

    free(sml);

    return ret;
}

int phgPartitionSearch(SFC_PTNS *ptn, int p, SFC_FLOAT key)
{
    int lo, hi;
    int md;
    int ret;

    assert(p >= 1);

    lo = 0;
    hi = p - 1;
    md = (lo + hi) >> 1;
    while (1) {
        if (key >= ptn[md].high)
            ret = 1;
        else if (key - ptn[md].low < 0.)
            ret = -1;
        else
            ret = 0;

        if (ret == 1) {
            lo = md;
            md = (lo + hi) >> 1;
            if (md == lo)
                md = hi;
        }
        else if (ret == -1) {
            hi = md;
            md = (lo + hi) >> 1;
        }
        else {
            break;
        }
    }

    return md;
}
