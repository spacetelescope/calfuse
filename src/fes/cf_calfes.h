#define FES_GOOD_PIX    0.0
#define FES_BAD_PIX   -20.0
#define FES_FLAT          1
#define FES_MASK          2
#define FES_BIAS          3

int cf_fes_apply_flat(fitsfile *, float **, float *, int, int);
int cf_fes_apply_mask(fitsfile *, float **, float *, int, int);
int cf_fes_apply_bias(fitsfile *, float **, float *, int, int);

