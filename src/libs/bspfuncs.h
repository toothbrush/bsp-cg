void bspmv(int p, int s, int n, int nz, int nrows, int ncols,
           double *a, int *inc,
           int *srcprocv, int *srcindv, int *destprocu, int *destindu,
           int nv, int nu, double *v, double *u);

int nloc(int p, int s, int n);

void bspmv_init(int p, int s, int n, int nrows, int ncols,
                int nv, int nu, int *rowindex, int *colindex,
                int *vindex, int *uindex, int *srcprocv, int *srcindv,
                int *destprocu, int *destindu);

double bspip(int p,int s,int nv1, int nv2, double* v1, int*v1index,
             double *v2, int *procv2, int *indv2);

void addvec(int nv, double *v,int*vindex, int nr, double *remote,
        int *procr, int *indr);
void copyvec(int s,
        int nv, int nu, double* v, double* u, int* uindex, int* procu, int* indu);
