/** @file test.c
 * @brief Unit Test File
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <stdio.h>
#include <assert.h> // assert
#include <math.h>   // fabsf

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "linalg.h"
#include "navtoolbox.h"

/******************************************************************************
 * DEFINES
 ******************************************************************************/

#define TEST_FLOAT_WITHIN(delta, expected, actual, message)                                        \
    assert((fabsf((expected) - (actual)) <= delta) && message)

/******************************************************************************
 * TYPEDEFS
 ******************************************************************************/

/******************************************************************************
 * LOCAL DATA DEFINITIONS
 ******************************************************************************/

/******************************************************************************
 * LOCAL FUNCTION PROTOTYPES
 ******************************************************************************/

/** @brief Fill array with a Hilbert matrix.
 * @param[out] H Output Hilbert matrix (n x n).
 * @param[in] n Dimension of H. */
static void hilbert(float* H, int n);

/** @brief Print a matrix to stdout
 * @param[in] R column-major n x m matrix
 * @param[in] n rows
 * @param[in] m cols
 * @param[in] fmt printf format string, e.g. "%.3f"
 * @param[in] name Pretty print with the name of the matrix, can be NULL */
static void matprint(const float* R, const int n, const int m, const char* fmt, const char* name);

/******************************************************************************
 * FUNCTION BODIES
 ******************************************************************************/

/* @} */

static void testlinalg(void)
{
    printf("Running linalg (linear algebra) tests...\n");
    /* Note: all matrices in column-major order */

    // Test Matrix Multiplication
    {
        const float A[4]    = { 1, 4, 3, 2 };             // A = 2 rows x 2 columns
        const float B[6]    = { 2, 3, 5, 6, 3, 9 };       // B= 2 rows x 3 columns
        float       C[6]    = { 1, 1, 1, 1, 1, 1 };       // Output: A*B --> 2 rows x 3 columns
        const float Cexp[6] = { 35, 44, 71, 98, 92, 92 }; // C = 3*A*B + 2*C
        const float alpha   = 3.0f;
        const float beta    = 2.0f;
        matmul("N", "N", 2, 3, 2, alpha, A, B, beta, C);
        for (int i = 0; i < 6; i++)
        {
            TEST_FLOAT_WITHIN(1.0e-08f, Cexp[i], C[i], "Error in matrix multiplication");
        }
        printf("[x] test C = alpha*A*B + beta*C (matmul)\n");
    }
    {
        const float A[6]    = { 9, 6, -5, 10, -3, 9 }; // A = 2 rows x 3 columns
        float       B[4]    = { -1, -1, -1, -1 };      // Output: A*A' --> 2 rows x 2 columns
        const float Bexp[4] = { 115, -23, -23, 217 };  // B = 1*A*A' + 0*B
        const float alpha   = 1.0f;
        const float beta    = 0.0f;
        matmul("N", "T", 2, 2, 3, alpha, A, A, beta, B);
        for (int i = 0; i < 4; i++)
        {
            TEST_FLOAT_WITHIN(1.0e-08f, Bexp[i], B[i], "Error in matrix multiplication");
        }
        printf("[x] test C = A*A' (matmul)\n");
    }
    {
        const float A[4]    = { 1, 4, 3, 2 };       // A = 2 rows x 2 columns
        const float B[6]    = { 2, 3, 5, 6, 3, 9 }; // B= 2 rows x 3 columns
        float       C[6]    = { 2, 2, 2, 2, 2, 2 }; // Output: A*B --> 2 rows x 3 columns
        const float Cexp[6] = { 21, 18, 43.5f, 40.5f, 58.5f, 40.5f }; // C = 1.5*A'*B + 0*C
        const float alpha   = 1.5f;
        const float beta    = 0.0f;
        matmul("T", "N", 2, 3, 2, alpha, A, B, beta, C);
        for (int i = 0; i < 6; i++)
        {
            TEST_FLOAT_WITHIN(1.0e-08f, Cexp[i], C[i], "Error in matrix multiplication");
        }
        printf("[x] test C = alpha*A'*B (matmul)\n");
    }
    {
        const float A[]    = { 1, -10, 5, 3, -20, 7 }; // A = 3 rows x 2 columns
        const float B[]    = { -1, 4, -2, 5, -3, 6 };  // B= 2 rows x 3 columns
        float       C[]    = { 3, 3, 3, 3 };
        const float Cexp[] = { 4, 16, -16, -46 }; // C = 1*A'*B'
        const float alpha  = 1.0f;
        const float beta   = 0.0f;

        matmul("T", "T", 2, 2, 3, alpha, A, B, beta, C);
        for (int i = 0; i < 4; i++)
        {
            TEST_FLOAT_WITHIN(1.0e-08f, Cexp[i], C[i], "Error in matrix multiplication");
        }
        printf("[x] test C = A'*B' (matmul)\n");
    }
    // Test cholesky decomposition
    {
        #define CHOLESKY_TEST_N     8
        const int n = CHOLESKY_TEST_N;
        float     L[CHOLESKY_TEST_N*CHOLESKY_TEST_N];
        hilbert(L, n); // L = hilbert(n) put hilbert matrix into L
        // matprint(L, n, n, "%10.8f", "H");
        const int result = cholesky(L, n, 0); // L = chol(L) inplace calc.
        assert(result == 0 && "Cholesky calculation test failed");
        // matprint(L, n, n, "%10.8f", "L (L*L'=H)");
        // test if L*L' actually is equal to H:
        float LLt[CHOLESKY_TEST_N*CHOLESKY_TEST_N];
        matmul("N", "T", n, n, n, 1.0f, L, L, 0.0, LLt); // LLt = L*L'
        float H[CHOLESKY_TEST_N*CHOLESKY_TEST_N];
        hilbert(H, n); // recreate expected result
        matprint(H, n, n, "%8.6f", "H");
        const float threshold = 1.5e-08f; // comparable error of independent MATLAB test
        for (int i = 0; i < n * n; i++)
        {
            TEST_FLOAT_WITHIN(threshold, H[i], LLt[i], "cholesky decomp. of Hilbert matrix failed");
        }
        printf("[x] Cholesky decomposition on close to singular matrix "
               "(cholesky)\n");
    }
    // Right-hand side triangular solve
    {
        const float L[]    = { 2, 3, 0, 1 }; // 2 x 2 matrix
        float       B[]    = { 8, 18, 28, 2, 4, 6 }; // 3 x 2 matrix
        float       Xexp[] = { 1, 3, 5, 2, 4, 6 }; // X*L = B
        const int result = trisolveright(L, B, 2, 3, "N");
        assert(result == 0);
        const float threshold = 1.0e-08f;
        for (int i = 0; i < 2 * 3; i++)
        {
            TEST_FLOAT_WITHIN(threshold, B[i], Xexp[i], "trisolveright failed");
        }
        printf("[x] Right-hand side triangular solve (trisolveright)\n");
    }
    {
        const float L[]    = { 2, 3, 0, 1 }; // 2 x 2 matrix
        float       B[]    = { 12, 10, 8, 21, 17, 13 }; // 3 x 2 matrix
        float       Xexp[] = { 6, 5, 4, 3, 2, 1 }; // X*L' = B
        const int result = trisolveright(L, B, 2, 3, "T");
        assert(result == 0);
        const float threshold = 1.0e-08f;
        for (int i = 0; i < 2 * 3; i++)
        {
            TEST_FLOAT_WITHIN(threshold, B[i], Xexp[i], "trisolveright with transpose failed");
        }
        printf("[x] Right-hand side triangular solve with transpose (trisolveright)\n");
    }
    {
        float L[3*3];
        float Xexp[3*3];
        float B[3*3];
        hilbert(L, 3);
        hilbert(Xexp, 3);
        cholesky(L, 3, 0);
        matmul("N", "N", 3, 3, 3, 1.0f, Xexp, L, 0.0f, B);
        // matprint(Xexp, 3, 3, "%8.6f", "X");
        // matprint(L, 3, 3, "%8.6f", "L");
        // matprint(B, 3, 3, "%8.6f", "B");
        const int result = trisolveright(L, B, 3, 3, "N");
        assert(result == 0);
        const float threshold = 1.0e-07f;
        for (int i = 0; i < 3 * 3; i++)
        {
            TEST_FLOAT_WITHIN(threshold, B[i], Xexp[i], "trisolveright with transpose failed");
        }
        printf("[x] Right-hand side triangular solve test case #2 (trisolveright)\n");
    }
}

static void testnavtoolbox(void)
{
    printf("Running navtoolbox tests...\n");

    // Test Roll Pitch From Accelerometer, body2nav
    {
        const float f_body[3] = { 0.0f, 0.0f, -0.01f }; /* close to free fall */
        float       roll_rad, pitch_rad;
        nav_roll_pitch_from_accelerometer(f_body, &roll_rad, &pitch_rad);
        TEST_FLOAT_WITHIN(DEG2RAD(1.0e-06f), DEG2RAD(0.0f), roll_rad,
                          "Roll angle calculation incorrect");
        TEST_FLOAT_WITHIN(DEG2RAD(1.0e-06f), DEG2RAD(0.0f), pitch_rad,
                          "Pitch angle calculation incorrect");
    }
    {
        const float f_nav[3] = { 0.0f, 0.0f, -GRAVITY };
        float       R[9];
        float       f_body[3];
        nav_matrix_body2nav(DEG2RAD(10.0f), DEG2RAD(20.0f), 0.0f, R);
        // matprint(R, 3, 3, "%6.3f", "R");
        matmul("T", "N", 3, 1, 3, 1.0f, R, f_nav, 0.0f, f_body);
        // printf("f_body = %.2f %.2f %.2f\n", f_body[0], f_body[1], f_body[2]);
        float roll_rad, pitch_rad;
        nav_roll_pitch_from_accelerometer(f_body, &roll_rad, &pitch_rad);
        // printf("%.1f %.1f\n", RAD2DEG(roll_rad), RAD2DEG(pitch_rad));
        TEST_FLOAT_WITHIN(DEG2RAD(1.0e-06f), DEG2RAD(10.0f), roll_rad,
                          "Roll angle calculation incorrect");
        TEST_FLOAT_WITHIN(DEG2RAD(1.0e-06f), DEG2RAD(20.0f), pitch_rad,
                          "Pitch angle calculation incorrect");
        printf("[x] Body to navigation frame transformation "
               "(nav_matrix_body2nav)\n");
        printf("[x] Initial alignment from accelerometer "
               "(nav_roll_pitch_from_accelerometer)\n");
    }
}

int main(int argc, char** argv)
{
    testlinalg();
    testnavtoolbox();

    printf("\n[OK] All tests completed.\n");

    return 0;
}

static void hilbert(float* H, int n)
{
    /*
     *   Hilbert test matrix, lines are _almost_ linear dependent
     *
     *   [  1       1/2     1/3     ...  1/nA      ]
     *   [  1/2     1/3     1/4     ...  1/(n+1)   ]
     *   [  1/3     1/4     1/5     ...  1/(n+2)   ]
     *   [               ...                       ]
     *   [  1/n     1/(n+1) 1/(n+2) ...  1/(2*n-1) ]
     */

    int start = 1;
    for (int i = 0; i < n; i++) /* row */
    {
        int rowstart = start;
        for (int j = 0; j < n; j++) /* col */
        {
            MAT_ELEM(H, i, j, n, n) = 1.0f / rowstart;
            rowstart++;
        }
        start++;
    }
}

static void matprint(const float* R, const int n, const int m, const char* fmt, const char* name)
{
    if (name)
    {
        printf(" %s =\n", name);
        printf("\t");
    }
    for (int i = 0; i < n; i++) /* row */
    {
        for (int j = 0; j < m; j++) /* col */
        {
            printf(fmt, (double)MAT_ELEM(R, i, j, n, m));
            printf(" ");
        }
        printf("\n");
        if (name && i < (n - 1))
        {
            printf("\t");
        }
    }
}
