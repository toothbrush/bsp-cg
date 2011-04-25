void bspinputvec(int p, int s, const char *filename,
                 int *pn, int *pnv, int **pvindex,
                 double **pvalues);
void triple2icrs(int n, int nz, int *ia,  int *ja, double *a,
                 int *pnrows, int *pncols,
                 int **prowindex, int **pcolindex);
void bspinput2triple(char*filename, int p, int s, int *pnA, int *pnz, 
                     int **pia, int **pja, double **pa);
typedef struct {int i,j;} indexpair;
#define STRLEN 100
#define DIV 0
#define MOD 1

