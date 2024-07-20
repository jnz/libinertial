/** @file linalg.h
 * @brief Embedded linear algebra math library
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

/******************************************************************************
 * DEFINES
 ******************************************************************************/

#define MAT_ELEM(M, row, col, numrows, numcols) (M[row + col * numrows])
#define SQRTF(x) sqrtf(x)

/******************************************************************************
 * TYPEDEFS
 ******************************************************************************/

/******************************************************************************
 * LOCAL DATA DEFINITIONS
 ******************************************************************************/

/******************************************************************************
 * LOCAL FUNCTION PROTOTYPES
 ******************************************************************************/

/******************************************************************************
 * FUNCTION PROTOTYPES
 ******************************************************************************/

/*** @brief matrix multiply C = alpha*A*B + beta*C
 * BLAS equivalent: ?gemm
 *
 * @param[in] ta Supply "T" (transpose A) or "N" (don't transpose A)
 * @param[in] tb Supply "T" (transpose B) or "N" (don't transpose B)
 * @param[in] n Dimension n (rows of A)
 * @param[in] k Dimension k (cols of B)
 * @param[in] m Dimension m (rows of B)
 * @param[in] alpha Factor alpha
 * @param[in] A Input matrix A (n x m)
 * @param[in] B Input matrix B (m x k)
 * @param[in] beta Factor beta
 * @param[in/out] C Output matrix (n x k)
 */
void matmul(const char* ta, const char* tb, int n, int k, int m, float alpha, const float* A,
            const float* B, float beta, float* C);

/** @brief Calculate the lower triangular matrix L, so
 * that L*L' = A.
 * Operation count: n^3/6 with n square roots.
 * BLAS equivalent: ?potrf
 * Optmized for
 *
 * @param[in,out] A Symmetric, positive definite (n x n) matrix. Only upper
 *                  triangular part needs to be given.
 *                  Lower triangular part is overwritten with L.
 * @param[in] n Dimension of A
 * @param[in] onlyWriteLowerPart If set to 0, overwrite the upper
 *                               triangular part with zeros.
 *                               Set to e.g. -1 to leave the upper
 *                               triangular part untouched.
 * @return 0 if successful, -1 if matrix is not positive definite.
 */
int cholesky(float* A, const int n, int onlyWriteLowerPart);

int strsm(const char* side, const char* uplo, const char* transa, const char* diag, int* m, int* n,
          float* alpha, float* a, int* lda, float* b, int* ldb);

/* @} */
