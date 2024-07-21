/** @file navtoolbox.c
 * @brief Navigation Toolbox Helper Functions
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <math.h>
#include <assert.h>

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "linalg.h"
#include "navtoolbox.h"

/******************************************************************************
 * DEFINES
 ******************************************************************************/

#define NAV_KALMAN_MAX_STATE_SIZE           32
#define NAV_KALMAN_MAX_MEASUREMENTS         3

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
 * FUNCTION BODIES
 ******************************************************************************/

void nav_roll_pitch_from_accelerometer(const float f[3], float* roll_rad, float* pitch_rad)
{
    assert(f);
    if (roll_rad)
    {
        *roll_rad = atan2f(-f[1], -f[2]); /* eq. 5.89 a */
    }
    if (pitch_rad)
    {
        *pitch_rad = atan2f(f[0], SQRTF(f[1] * f[1] + f[2] * f[2])); /* eq. 5.89 b */
    }
}

void nav_matrix_body2nav(const float roll_rad, const float pitch_rad, const float yaw_rad,
                         float R_output[9])
{
    const float sinr = sinf(roll_rad);
    const float sinp = sinf(pitch_rad);
    const float siny = sinf(yaw_rad);
    const float cosr = cosf(roll_rad);
    const float cosp = cosf(pitch_rad);
    const float cosy = cosf(yaw_rad);
    R_output[0]      = cosp * cosy;
    R_output[3]      = sinr * sinp * cosy - cosr * siny;
    R_output[6]      = cosr * sinp * cosy + sinr * siny;
    R_output[1]      = cosp * siny;
    R_output[4]      = sinr * sinp * siny + cosr * cosy;
    R_output[7]      = cosr * sinp * siny - sinr * cosy;
    R_output[2]      = -sinp;
    R_output[5]      = sinr * cosp;
    R_output[8]      = cosr * cosp;
}

int nav_kalman(float* x, float* P, const float* dl, const float* R, const float* H, int n, int m)
{
    int result;
    float D[NAV_KALMAN_MAX_STATE_SIZE   * NAV_KALMAN_MAX_MEASUREMENTS];
    float L[NAV_KALMAN_MAX_MEASUREMENTS * NAV_KALMAN_MAX_MEASUREMENTS];
    float E[NAV_KALMAN_MAX_STATE_SIZE   * NAV_KALMAN_MAX_MEASUREMENTS];
    float K[NAV_KALMAN_MAX_STATE_SIZE   * NAV_KALMAN_MAX_MEASUREMENTS];

    assert(n <= NAV_KALMAN_MAX_STATE_SIZE);
    assert(m <= NAV_KALMAN_MAX_MEASUREMENTS);

    /*  (1) D = P * A'                  gemm or symm
     *  (2) S = A * D + R               gemm
     *  (3) L = chol(S) (L*L'=S)        potrf
     *  (4) E = D * U ^ -1              trsm
     *  (5) P = P - E * E '             syrk
     *  (6) K = E * U ^ -T              trsm
     *  (7) x = x + K * v               gemm or gemv
     *
     *  n = state variables
     *  m = measurements
     *
     *  Matrix dimensions:
     *  D = n x m
     *  S = m x m
     *  U = m x m
     *  E = n x m
     *  K = n x m */

    /*
    matmul("N", "T", n, m, n, 1.0f, P, H, 0.0f, D);  // (1)
    memcpy(R, size(float)*m*m, U);                   // U = S = R
    matmul("N", "N", m, m, n, 1.0f, A, D, 1.0f, U);  // (2)
    result = cholesky(
    */

}

#if 0
    matmul("NN",n,m,n,1.0,P,H,0.0,F);
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
#endif

/* @} */
