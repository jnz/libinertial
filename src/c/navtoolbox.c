/** @file navtoolbox.c
 * libinertial, Jan Zwiener (jan@zwiener.org)
 *
 * @brief Navigation Toolbox Helper Functions
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <math.h>
#include <assert.h>
#include <string.h> /* memcpy */

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "linalg.h"
#include "navtoolbox.h"

/******************************************************************************
 * DEFINES
 ******************************************************************************/

#ifndef NAV_KALMAN_MAX_STATE_SIZE
#define NAV_KALMAN_MAX_STATE_SIZE 32 /* kalman filter scratchpad buf size */
#endif
#ifndef NAV_KALMAN_MAX_MEASUREMENTS
#define NAV_KALMAN_MAX_MEASUREMENTS 3 /* kalman filter scratchpad buf size */
#endif

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

int nav_kalman(float* x, float* P, const float* dz, const float* R, const float* Ht, int n, int m)
{
    float D[NAV_KALMAN_MAX_STATE_SIZE * NAV_KALMAN_MAX_MEASUREMENTS];
    float L[NAV_KALMAN_MAX_MEASUREMENTS * NAV_KALMAN_MAX_MEASUREMENTS];
    assert(n > 0 && n <= NAV_KALMAN_MAX_STATE_SIZE);
    assert(m > 0 && m <= NAV_KALMAN_MAX_MEASUREMENTS);

    /*  (1) D = P * H'              symm           |   Matrix dimensions:
     *  (2) S = H * D + R           gemm           |   D = n x m
     *  (3) L = chol(S) (L*L'=S)    potrf          |   S = m x m
     *  (4) E = D * L^-T            trsm           |   L = m x m
     *  (5) P = P - E*E'            syrk           |   E = n x m
     *  (6) K = E * L^-1            trsm           |   K = n x m
     *  (7) x = x + K*dz            gemm
     *
     *  n = state variables
     *  m = measurements */

    /* Inplace cholesky decomposition */
    /* Only update the required triangular parts (save instructions and memory access) */
    /* keep symmetry */
    /* numerically stable */

    matmulsym(P, Ht, n, m, D); // (1) D = P * H' (using upper triangular part of P)
    memcpy(L /*dst*/, R /*src*/, sizeof(float) * m * m); // Use L as temp. matrix, preload R
    matmul("T", "N", m, m, n, 1.0f, Ht, D, 1.0f, L);     // (2) L += H*D
    int result =
        cholesky(L, m, 1 /*don't fill upper triangular part of L*/); // (3) L = chol(H*D + R)
                                                                     // (inplace calculation of L)
    if (result != 0)
    {
        return -1;
    }                               // Cholesky fails: bail out (*)
    trisolveright(L, D, m, n, "T"); // (4) given L' and D, solve E*L' = D, for E, overwrite D with E
    symmetricrankupdate(P, D /*E*/, n, m); // (5) P = P - E*E'
    trisolveright(L, D /*E*/, m, n, "N");  // (6) solve K*L = E, for K, overwrite D with K
    matmul("N", "N", n, 1, m, 1.0f, D /*K*/, dz, 1.0f, x); // (7) x = x + K * dz (K is stored in D)

    /* FIXME check for P positive definite (symmetric is automatic)*/
    /* FIXME check for isfinite() in state vector */
    /* FIXME check for chi2 */

    /* (*) If a Cholesky decomposition is found the trsm operations will succeed. */

    return 0;
}

int nav_kalman_udu_scalar(float* x, float* U, float* d, const float dz, const
                          float R, const float* H_line, int n)
{
    float a[NAV_KALMAN_MAX_STATE_SIZE];
    float b[NAV_KALMAN_MAX_STATE_SIZE];
    float alpha = R;
    float gamma = 1.0f/alpha;

    matmul("T", "N", n, 1, n, 1.0f, U, H_line, 0.0f, a); // a = U'*H'
    for (int j=0;j<n;j++)
    {
        b[j] = d[j]*a[j]; // b = D*a = diag(d)*a
    }

    for (int j=0;j<n;j++)
    {
        float beta=alpha;
        alpha += a[j]*b[j];
        float lambda = -a[j]*gamma;
        gamma = 1.0f/alpha; // FIXME test if this is possible and return -1 if not
        d[j] *= beta*gamma;
        for (int i=0;i<j;i++)
        {
            beta = MAT_ELEM(U, i, j, n, n);
            MAT_ELEM(U, i, j, n, n) = beta + b[i]*lambda;
            b[i] += b[j]*beta;
        }
    }

    for (int j=0;j<n;j++)
    {
        x[j] += gamma*dz*b[j];
    }

    return 0;
}

int nav_kalman_udu(float* x, float* U, float* d, const float* z, const float* R,
                   const float* Ht, int n, int m, float chi2_threshold, int downweight_outlier)
{
    int retcode = 0;

    for (int i=0;i<m;i++,Ht+=n) /* iterate over each measurement,
                                   goto next line of H after each iteration */
    {
        float Rv = MAT_ELEM(R, i, i, m, m); /// get scalar measurement variance
        float dz = z[i]; // calculate residual for current scalar measurement
        matmul("N", "N", 1, 1, n, -1.0f, Ht, x, 1.0f, &dz); // dz = z - H(i,:)*x

        // <robust>
        if (chi2_threshold > 0.0f)
        {
            float tmp[NAV_KALMAN_MAX_STATE_SIZE];
            float s; // for chi2 test: s = H*U*diag(d)*U'*H' + R
                     // Chang, G. (2014). Robust Kalman filtering based on
                     // Mahalanobis distance as outlier judging criterion.
                     // Journal of Geodesy, 88(4), 391-401.

            float HPHT = 0.0f; // calc. scalar result of H_line*U*diag(d)*U'*H_line'
            matmul("N", "N", 1, n, n, 1.0f, Ht, U, 0.0f, tmp); // tmp = H(i,:) * U
            for (int j = 0; j < n; j++)
            {
                HPHT += tmp[j] * tmp[j] * d[j];
            }
            s = HPHT + Rv;
            const float mahalanobis_dist_sq = dz * dz / s;
            if (mahalanobis_dist_sq > chi2_threshold) // potential outlier?
            {
                if (!downweight_outlier)
                {
                    continue; // just skip this measurement
                }
                // process this measurement, but reduce the measurement precision
                const float f = mahalanobis_dist_sq / chi2_threshold;
                Rv            = (f - 1.0f) * HPHT + f * Rv;
            }
        }
        // </robust>

        int status = nav_kalman_udu_scalar(x, U, d, dz, Rv, Ht, n);
        if (status != 0)
        {
            retcode = -1; // still process rest of the measurement vector
        }
    }
    return retcode;
}

int decorrelate(float* z, float* Ht, float* R, int n, int m)
{
    /* Basic decorrelation in MATLAB
    [G] = chol(R); % G'*G = R
    zdecorr = (G')\z;
    Hdecorr = (G')\H;
    Rdecorr = eye(length(z)); */

    // in-place cholesky so that L*L' = R:
    int result = cholesky(R, m, 0/* 0 means: fill upper part with zeros */);
    if (result != 0) { return -1; }
    // L*H_decorr = H
    // (L*H_decorr)' = H'
    // H_decorr'*L' = H' solve for H_decorr
    trisolveright(R/*L*/, Ht, m, n, "T");
    trisolve(R/*L*/, z, m, 1, "N");

    return 0;
}

/* @} */
