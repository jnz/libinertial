/** @file linalg.c
 * @brief Math functions
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <ctype.h> /* tolower() */

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "linalg.h"

/******************************************************************************
 * DEFINES
 ******************************************************************************/

/******************************************************************************
 * TYPEDEFS
 ******************************************************************************/

/******************************************************************************
 * LOCAL DATA DEFINITIONS
 ******************************************************************************/

/******************************************************************************
 * LOCAL FUNCTION PROTOTYPES
 ******************************************************************************/

// inline int min(int a, int b) { return a < b ? a : b; }
inline int max(int a, int b)
{
    return a > b ? a : b;
}

/* BLAS/LAPACK routines */
static int lsame(const char* a, const char* b);

/******************************************************************************
 * FUNCTION BODIES
 ******************************************************************************/

void matmul(const char* ta, const char* tb, int n, int k, int m, float alpha, const float* A,
            const float* B, float beta, float* C)
{
    const int ca     = lsame(ta, "T");
    const int cb     = lsame(tb, "T");
    const int branch = (ca << 1) | (cb);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < k; j++)
        {
            float s = 0.0f;
            switch (branch)
            {
            /* clang-format off */
                case 0: for (int x=0;x<m;x++) { s+=A[i+x*n]*B[x+j*m]; } break; /* N N */
                case 1: for (int x=0;x<m;x++) { s+=A[i+x*n]*B[j+x*k]; } break; /* N T */
                case 2: for (int x=0;x<m;x++) { s+=A[x+i*m]*B[x+j*m]; } break; /* T N */
                case 3: for (int x=0;x<m;x++) { s+=A[x+i*m]*B[j+x*k]; } break; /* T T */
                /* clang-format on */
            }
            if (beta == 0.0f)
            {
                MAT_ELEM(C, i, j, n, n) = alpha * s;
            }
            else
            {
                MAT_ELEM(C, i, j, n, n) = alpha * s + beta * C[i + j * n];
            }
        }
    }
}

int cholesky(float* A, const int n, int onlyWriteLowerPart)
{
    /* in-place calculation of lower triangular matrix L*L' = A */

    /* set the upper triangular part to zero? */
    if (!onlyWriteLowerPart)
    {
        for (int i = 0; i < n - 1; i++) /* row */
        {
            for (int j = i + 1; j < n; j++) /* col */
            {
                MAT_ELEM(A, i, j, n, n) = 0.0f;
            }
        }
    }

    for (int j = 0; j < n; j++) /* main loop */
    {
        const float Ajj = MAT_ELEM(A, j, j, n, n);
        if (Ajj <= 0.0f || isnan(Ajj))
        {
            return -1;
        }
        MAT_ELEM(A, j, j, n, n) = SQRTF(Ajj);

        const float invLjj = 1.0f / MAT_ELEM(A, j, j, n, n);
        for (int i = j + 1; i < n; i++)
        {
            MAT_ELEM(A, i, j, n, n) *= invLjj;
        }

        for (int k = j + 1; k < n; k++)
        {
            for (int i = k; i < n; i++)
            {
                MAT_ELEM(A, i, k, n, n) -= MAT_ELEM(A, i, j, n, n) * MAT_ELEM(A, k, j, n, n);
            }
        }
    }
    return 0;
}

/**
 * @brief Triangular solve (right hand side).
 *
 * Solve matrix equation: X*L = A
 * @param[in]     L Given lower triangular matrix (dimension m x m)
 * @param[in,out] A Matrix being overwritten by X (dimension n x m)
 * @param[in]     n Matrix dimension
 * @param[in]     m Matrix dimension
 * @param[in]     tp Transpose L?
 */
int trisolveright(const float* L, float* A,
                  const int n, const int m, const char* tp)
{

}

/******************************************************************************
 * Local BLAS/LAPACK
 ******************************************************************************/

static int lsame(const char* a, const char* b)
{
    return (tolower(*a) == tolower(*b));
}

int strsm(const char* side, const char* uplo, const char* transa, const char* diag, int* m, int* n,
          float* alpha, float* a, int* lda, float* b, int* ldb)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    int   i__, j, k, info;
    float temp;
    int   lside;
    int   nrowa;
    int   upper;
    int   nounit;

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  STRSM  solves one of the matrix equations */

    /*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, */

    /*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
    /*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */

    /*     op( A ) = A   or   op( A ) = A'. */

    /*  The matrix X is overwritten on B. */

    /*  Arguments */
    /*  ========== */

    /*  SIDE   - CHARACTER*1. */
    /*           On entry, SIDE specifies whether op( A ) appears on the left */
    /*           or right of X as follows: */
    /*              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */
    /*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */

    /*           Unchanged on exit. */

    /*  UPLO   - CHARACTER*1. */
    /*           On entry, UPLO specifies whether the matrix A is an upper or */
    /*           lower triangular matrix as follows: */
    /*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
    /*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
    /*           Unchanged on exit. */

    /*  TRANSA - CHARACTER*1. */
    /*           On entry, TRANSA specifies the form of op( A ) to be used in */
    /*           the matrix multiplication as follows: */
    /*              TRANSA = 'N' or 'n'   op( A ) = A. */
    /*              TRANSA = 'T' or 't'   op( A ) = A'. */
    /*              TRANSA = 'C' or 'c'   op( A ) = A'. */

    /*           Unchanged on exit. */

    /*  DIAG   - CHARACTER*1. */
    /*           On entry, DIAG specifies whether or not A is unit triangular */
    /*           as follows: */
    /*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
    /*              DIAG = 'N' or 'n'   A is not assumed to be unit */
    /*                                  triangular. */

    /*           Unchanged on exit. */

    /*  M      - INTEGER. */
    /*           On entry, M specifies the number of rows of B. M must be at */
    /*           least zero. */
    /*           Unchanged on exit. */

    /*  N      - INTEGER. */
    /*           On entry, N specifies the number of columns of B.  N must be */
    /*           at least zero. */
    /*           Unchanged on exit. */

    /*  ALPHA  - REAL            . */
    /*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
    /*           zero then  A is not referenced and  B need not be set before */
    /*           entry. */
    /*           Unchanged on exit. */

    /*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m */
    /*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
    /*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
    /*           upper triangular part of the array  A must contain the upper */
    /*           triangular matrix  and the strictly lower triangular part of */
    /*           A is not referenced. */
    /*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
    /*           lower triangular part of the array  A must contain the lower */
    /*           triangular matrix  and the strictly upper triangular part of */
    /*           A is not referenced. */
    /*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
    /*           A  are not referenced either,  but are assumed to be  unity. */
    /*           Unchanged on exit. */

    /*  LDA    - INTEGER. */
    /*           On entry, LDA specifies the first dimension of A as declared */
    /*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
    /*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' */
    /*           then LDA must be at least max( 1, n ). */
    /*           Unchanged on exit. */

    /*  B      - REAL             array of DIMENSION ( LDB, n ). */
    /*           Before entry,  the leading  m by n part of the array  B must */
    /*           contain  the  right-hand  side  matrix  B,  and  on exit  is */
    /*           overwritten by the solution matrix  X. */

    /*  LDB    - INTEGER. */
    /*           On entry, LDB specifies the first dimension of B as declared */
    /*           in  the  calling  (sub)  program.   LDB  must  be  at  least */
    /*           max( 1, m ). */
    /*           Unchanged on exit. */

    /*  Level 3 Blas routine. */

    /*  -- Written on 8-February-1989. */
    /*     Jack Dongarra, Argonne National Laboratory. */
    /*     Iain Duff, AERE Harwell. */
    /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
    /*     Sven Hammarling, Numerical Algorithms Group Ltd. */

    /* Parameter adjustments */
    a_dim1   = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1   = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = lsame(side, "L");
    if (lside)
    {
        nrowa = *m;
    }
    else
    {
        nrowa = *n;
    }
    nounit = lsame(diag, "N");
    upper  = lsame(uplo, "U");

    info = 0;
    if (!lside && !lsame(side, "R"))
    {
        info = 1;
    }
    else if (!upper && !lsame(uplo, "L"))
    {
        info = 2;
    }
    else if (!lsame(transa, "N") && !lsame(transa, "T") && !lsame(transa, "C"))
    {
        info = 3;
    }
    else if (!lsame(diag, "U") && !lsame(diag, "N"))
    {
        info = 4;
    }
    else if (*m < 0)
    {
        info = 5;
    }
    else if (*n < 0)
    {
        info = 6;
    }
    else if (*lda < max(1, nrowa))
    {
        info = 9;
    }
    else if (*ldb < max(1, *m))
    {
        info = 11;
    }
    if (info != 0)
    {
        return 0;
    }
    if (*m == 0 || *n == 0)
    { /*     Quick return if possible. */
        return 0;
    }
    if (*alpha == 0.f)
    { /*     And when  alpha.eq.zero. */
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                b[i__ + j * b_dim1] = 0.f;
            }
        }
        return 0;
    }

    /*     Start the operations. */
    if (lside)
    {
        if (lsame(transa, "N"))
        {
            /*           Form  B := alpha*inv( A )*B. */
            if (upper)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    for (k = *m; k >= 1; --k)
                    {
                        if (b[k + j * b_dim1] != 0.f)
                        {
                            if (nounit)
                            {
                                b[k + j * b_dim1] /= a[k + k * a_dim1];
                            }
                            i__2 = k - 1;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[i__ + k * a_dim1];
                            }
                        }
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    i__2 = *m;
                    for (k = 1; k <= i__2; ++k)
                    {
                        if (b[k + j * b_dim1] != 0.f)
                        {
                            if (nounit)
                            {
                                b[k + j * b_dim1] /= a[k + k * a_dim1];
                            }
                            i__3 = *m;
                            for (i__ = k + 1; i__ <= i__3; ++i__)
                            {
                                b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[i__ + k * a_dim1];
                            }
                        }
                    }
                }
            }
        }
        else
        {

            /*           Form  B := alpha*inv( A' )*B. */

            if (upper)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        temp = *alpha * b[i__ + j * b_dim1];
                        i__3 = i__ - 1;
                        for (k = 1; k <= i__3; ++k)
                        {
                            temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
                        }
                        if (nounit)
                        {
                            temp /= a[i__ + i__ * a_dim1];
                        }
                        b[i__ + j * b_dim1] = temp;
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    for (i__ = *m; i__ >= 1; --i__)
                    {
                        temp = *alpha * b[i__ + j * b_dim1];
                        i__2 = *m;
                        for (k = i__ + 1; k <= i__2; ++k)
                        {
                            temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
                        }
                        if (nounit)
                        {
                            temp /= a[i__ + i__ * a_dim1];
                        }
                        b[i__ + j * b_dim1] = temp;
                    }
                }
            }
        }
    }
    else
    {
        if (lsame(transa, "N"))
        {
            /*           Form  B := alpha*B*inv( A ). */
            if (upper)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    i__2 = j - 1;
                    for (k = 1; k <= i__2; ++k)
                    {
                        if (a[k + j * a_dim1] != 0.f)
                        {
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__)
                            {
                                b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (nounit)
                    {
                        temp = 1.f / a[j + j * a_dim1];
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        }
                    }
                }
            }
            else
            {
                for (j = *n; j >= 1; --j)
                {
                    if (*alpha != 1.f)
                    {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    i__1 = *n;
                    for (k = j + 1; k <= i__1; ++k)
                    {
                        if (a[k + j * a_dim1] != 0.f)
                        {
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (nounit)
                    {
                        temp = 1.f / a[j + j * a_dim1];
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        }
                    }
                }
            }
        }
        else /* Form  B := alpha*B*inv( A' ). */
        {
            if (upper)
            {
                for (k = *n; k >= 1; --k)
                {
                    if (nounit)
                    {
                        temp = 1.f / a[k + k * a_dim1];
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                        }
                    }
                    i__1 = k - 1;
                    for (j = 1; j <= i__1; ++j)
                    {
                        if (a[j + k * a_dim1] != 0.f)
                        {
                            temp = a[j + k * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                b[i__ + j * b_dim1] -= temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (*alpha != 1.f)
                    {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1];
                        }
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (k = 1; k <= i__1; ++k)
                {
                    if (nounit)
                    {
                        temp = 1.f / a[k + k * a_dim1];
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                        }
                    }
                    i__2 = *n;
                    for (j = k + 1; j <= i__2; ++j)
                    {
                        if (a[j + k * a_dim1] != 0.f)
                        {
                            temp = a[j + k * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__)
                            {
                                b[i__ + j * b_dim1] -= temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1];
                        }
                    }
                }
            }
        }
    }

    return 0;
}
