/** @file linalg.h
 * @brief Embedded linear algebra math library
 *
 * Note: all matrices are stored in column-major order.
 *
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
 * BLAS: ?gemm
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
 * that L*L' = A. Operation count: n^3/6 with n square roots.
 * BLAS equivalent: ?potrf
 *
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

/**
 * @brief Triangular solve (right hand side).
 * BLAS: ?trsm
 *
 * Solve matrix equation: X*L = A
 * @param[in]     L Given lower triangular matrix (dimension n x n)
 * @param[in,out] A Matrix being overwritten by X (dimension m x n)
 * @param[in]     n Matrix dimension (rows / columns of L)
 * @param[in]     m Matrix dimension (rows of A)
 * @param[in]     tp Transpose L?
 *
 * @return 0 = successful
 */
int trisolveright(const float* L, float* A, int n, int m, const char* tp);

/* @} */
