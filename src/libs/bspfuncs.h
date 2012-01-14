void bspmv(int p, int s, int n, int nz, int nrows, int ncols,
           double *a, int *inc,
           int *srcprocv, int *srcindv, int *destprocu, int *destindu,
           int nv, int nu, double *v, double *u);

int nloc(int p, int s, int n);

void bspmv_init(int p, int s, int n, int nrows, int ncols,
                int nv, int nu, int *rowindex, int *colindex,
                int *vindex, int *uindex, int *srcprocv, int *srcindv,
                int *destprocu, int *destindu);

double bspip(int p,int s,int nv1, int nv2, double* v1, double *v2,
        int *procv2, int *indv2);

void bspinprod();
