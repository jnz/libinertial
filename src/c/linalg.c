/** @file linalg.c
 * libinertial, Jan Zwiener (jan@zwiener.org)
 *
 * @brief Math functions
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <string.h> /* memset */
#include <ctype.h> /* tolower() */

#include <blasfeo.h>

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "linalg.h"
// #include "blasmini.h"

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

static int ssyrk_(char* uplo, char* trans, int* n, int* k, float* alpha, float* a, int* lda, float* beta,
          float* c__, int* ldc);

inline int max(int a, int b)
{
    return a > b ? a : b;
}

/******************************************************************************
 * FUNCTION BODIES
 ******************************************************************************/

static int lsame(const char* a, const char* b)
{
    return (tolower(*a) == tolower(*b));
}

void matmul(const char* ta, const char* tb, int n, int k, int m, float alpha, const float* A,
            const float* B, float beta, float* C)
{
    int lda    = lsame(ta, "T") ? m : n;
    int ldb    = lsame(tb, "T") ? k : m;
    sgemm_((char*)ta, (char*)tb, &n, &k, &m, &alpha, (float*)A, &lda, (float*)B, &ldb, &beta, C, &n);
}

void matmulsym(const float* A_sym, const float* B, int n, int m, float* C)
{
    float alpha  = 1.0f;
    float beta   = 0.0f;
#if 0
    int   result = ssymm_("L" /* calculate C = A*B not C = B*A */,
                         "U" /* reference upper triangular part of A */, &n, /* rows of B/C */
                         &m,                                                 /* cols of B / C */
                         &alpha, (float*)A_sym, &n, (float*)B, &n, &beta, C, &n);
#endif
    sgemm_("N", "N", &n, &m, &n, &alpha, (float*)A_sym, &n, (float*)B, &n, &beta, C, &n);
}

void mateye(float* A, int n)
{
    memset(A, 0, sizeof(float) * n * n);
    for (int i = 0; i < n; i++)
    {
        MAT_ELEM(A, i, i, n, n) = 1.0f;
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
        if (Ajj <= 0.0f || !isfinite(Ajj))
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

void trisolve(const float* A, float* B, int n, int m, const char* tp)
{
    float alpha = 1.0f;
    strsm_("L" /* left hand*/, "L" /* lower triangular matrix */, (char*)tp /* transpose L? */,
           "N" /* L is not unit triangular */, &n, &m, &alpha, (float*)A, &n, B, &n);
    /* strsm basically just checks for proper matrix dimensions, handle via assert */
}

void trisolveright(const float* L, float* A, int n, int m, const char* tp)
{
    float     alpha = 1.0f;
    strsm_("R" /* right hand*/, "L" /* lower triangular matrix */, (char*)tp /* transpose L? */,
           "N" /* L is not unit triangular */, &m, &n, &alpha, (float*)L, &n, A, &m);
    /* strsm basically just checks for proper matrix dimensions, handle via assert */
}

void symmetricrankupdate(float* P, const float* E, int n, int m)
{
    float alpha = -1.0f;
    float beta  = 1.0f;

    const int result = ssyrk_("U", "N", &n, &m, &alpha, (float*)E, &n, &beta, P, &n);
    assert(result == 0);
}

int udu(const float* A, float* U, float* d, const int m)
{
    /*    A = U*diag(d)*U' decomposition
     *    Source:
     *      1. Golub, Gene H., and Charles F. Van Loan. "Matrix Computations." 4rd ed.,
     *         Johns Hopkins University Press, 2013.
     *      2. Grewal, Weill, Andrews. "Global positioning systems, inertial
     *         navigation, and integration". 1st ed. John Wiley & Sons, New York, 2001.
     *
     *    function [U, d] = udu(M)
     *      [m, ~] = size(M);
     *      U = zeros(m, m); d = zeros(m, 1);
     *
     *      for j = m:-1:1
     *        for i = j:-1:1
     *          sigma = M(i, j);
     *          for k = j + 1:m
     *              sigma = sigma - U(i, k) * d(k) * U(j, k);
     *          end
     *          if i == j
     *              d(j) = sigma;
     *              U(j, j) = 1; % U is a unit triangular matrix
     *          else
     *              U(i, j) = sigma / d(j); % off-diagonal elements of U
     *          end
     *        end
     *      end
     *    end
     */
    int   i, j, k;
    float sigma;

    memset(U, 0, sizeof(U[0]) * m * m);
    memset(d, 0, sizeof(d[0]) * m);

    for (j = m - 1; j >= 0; j--) /* UDU decomposition */
    {
        for (i = j; i >= 0; i--)
        {
            sigma = MAT_ELEM(A, i, j, m, m);
            for (k = j + 1; k < m; k++)
            {
                sigma -= MAT_ELEM(U, i, k, m, m) * d[k] * MAT_ELEM(U, j, k, m, m);
            }
            if (i == j)
            {
                d[j]                    = sigma;
                MAT_ELEM(U, j, j, m, m) = 1.0f;
            }
            else
            {
                if ((d[j] <= 0.0f) || !isfinite(d[j]))
                {
                    /* matrix is not positive definite if d < 0 */
                    return -1;
                }
                MAT_ELEM(U, i, j, m, m) = sigma / d[j];
            }
        }
    }
    return 0;
}

static int ssyrk_(char* uplo, char* trans, int* n, int* k, float* alpha, float* a, int* lda, float* beta,
          float* c__, int* ldc)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    int   i__, j, l, info;
    float temp;
    int   nrowa;
    int   upper;

    /*  Purpose */
    /*  ======= */

    /*  SSYRK  performs one of the symmetric rank k operations */

    /*     C := alpha*A*A' + beta*C, */

    /*  or */

    /*     C := alpha*A'*A + beta*C, */

    /*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
    /*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix */
    /*  in the second case. */

    /*  Arguments */
    /*  ========== */

    /*  UPLO   - CHARACTER*1. */
    /*           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
    /*           triangular  part  of the  array  C  is to be  referenced  as */
    /*           follows: */

    /*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
    /*                                  is to be referenced. */

    /*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
    /*                                  is to be referenced. */

    /*           Unchanged on exit. */

    /*  TRANS  - CHARACTER*1. */
    /*           On entry,  TRANS  specifies the operation to be performed as */
    /*           follows: */

    /*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C. */

    /*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C. */

    /*              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C. */

    /*           Unchanged on exit. */

    /*  N      - INTEGER. */
    /*           On entry,  N specifies the order of the matrix C.  N must be */
    /*           at least zero. */
    /*           Unchanged on exit. */

    /*  K      - INTEGER. */
    /*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
    /*           of  columns   of  the   matrix   A,   and  on   entry   with */
    /*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
    /*           of rows of the matrix  A.  K must be at least zero. */
    /*           Unchanged on exit. */

    /*  ALPHA  - REAL            . */
    /*           On entry, ALPHA specifies the scalar alpha. */
    /*           Unchanged on exit. */

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is */
    /*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
    /*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
    /*           part of the array  A  must contain the matrix  A,  otherwise */
    /*           the leading  k by n  part of the array  A  must contain  the */
    /*           matrix A. */
    /*           Unchanged on exit. */

    /*  LDA    - INTEGER. */
    /*           On entry, LDA specifies the first dimension of A as declared */
    /*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
    /*           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
    /*           be at least  max( 1, k ). */
    /*           Unchanged on exit. */

    /*  BETA   - REAL            . */
    /*           On entry, BETA specifies the scalar beta. */
    /*           Unchanged on exit. */

    /*  C      - REAL             array of DIMENSION ( LDC, n ). */
    /*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
    /*           upper triangular part of the array C must contain the upper */
    /*           triangular part  of the  symmetric matrix  and the strictly */
    /*           lower triangular part of C is not referenced.  On exit, the */
    /*           upper triangular part of the array  C is overwritten by the */
    /*           upper triangular part of the updated matrix. */
    /*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
    /*           lower triangular part of the array C must contain the lower */
    /*           triangular part  of the  symmetric matrix  and the strictly */
    /*           upper triangular part of C is not referenced.  On exit, the */
    /*           lower triangular part of the array  C is overwritten by the */
    /*           lower triangular part of the updated matrix. */

    /*  LDC    - INTEGER. */
    /*           On entry, LDC specifies the first dimension of C as declared */
    /*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
    /*           max( 1, n ). */
    /*           Unchanged on exit. */

    /*  Level 3 Blas routine. */

    /*  -- Written on 8-February-1989. */
    /*     Jack Dongarra, Argonne National Laboratory. */
    /*     Iain Duff, AERE Harwell. */
    /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
    /*     Sven Hammarling, Numerical Algorithms Group Ltd. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1   = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1   = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame(trans, "N"))
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    upper = lsame(uplo, "U");

    info = 0;
    if (!upper && !lsame(uplo, "L"))
    {
        info = 1;
    }
    else if (!lsame(trans, "N") && !lsame(trans, "T") && !lsame(trans, "C"))
    {
        info = 2;
    }
    else if (*n < 0)
    {
        info = 3;
    }
    else if (*k < 0)
    {
        info = 4;
    }
    else if (*lda < max(1, nrowa))
    {
        info = 7;
    }
    else if (*ldc < max(1, *n))
    {
        info = 10;
    }
    if (info != 0)
    {
        return info;
    }

    /*     Quick return if possible. */

    if (*n == 0 || ((*alpha == 0.f || *k == 0) && *beta == 1.f))
    {
        return 0;
    }

    /*     And when  alpha.eq.zero. */

    if (*alpha == 0.f)
    {
        if (upper)
        {
            if (*beta == 0.f)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
        else
        {
            if (*beta == 0.f)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
        return 0;
    }

    /*     Start the operations. */

    if (lsame(trans, "N"))
    {

        /*        Form  C := alpha*A*A' + beta*C. */

        if (upper)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (*beta == 0.f)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    if (a[j + l * a_dim1] != 0.f)
                    {
                        temp = *alpha * a[j + l * a_dim1];
                        i__3 = j;
                        for (i__ = 1; i__ <= i__3; ++i__)
                        {
                            c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
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
                if (*beta == 0.f)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    if (a[j + l * a_dim1] != 0.f)
                    {
                        temp = *alpha * a[j + l * a_dim1];
                        i__3 = *n;
                        for (i__ = j; i__ <= i__3; ++i__)
                        {
                            c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                        }
                    }
                }
            }
        }
    }
    else
    {

        /*        Form  C := alpha*A'*A + beta*C. */

        if (upper)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l)
                    {
                        temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
                    }
                    if (*beta == 0.f)
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    }
                    else
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = *n;
                for (i__ = j; i__ <= i__2; ++i__)
                {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l)
                    {
                        temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
                    }
                    if (*beta == 0.f)
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    }
                    else
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
    }

    return 0;
}
