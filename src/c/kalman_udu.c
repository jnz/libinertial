/** @file kalman_udu.c
 * @author Jan Zwiener (jan@zwiener.org)
 *
 * @brief UDU Kalman Filter
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
#include "blasmini.h" /* strmm_ */

/******************************************************************************
 * DEFINES
 ******************************************************************************/

#ifndef KALMAN_MAX_STATE_SIZE
#define KALMAN_MAX_STATE_SIZE 32 /* kalman filter scratchpad buf size */
#endif
#ifndef KALMAN_MAX_MEASUREMENTS
#define KALMAN_MAX_MEASUREMENTS 3 /* kalman filter scratchpad buf size */
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

int kalman_udu_scalar(float* x, float* U, float* d, const float dz, const float R,
                      const float* H_line, int n)
{
    float a[KALMAN_MAX_STATE_SIZE];
    float b[KALMAN_MAX_STATE_SIZE];
    float alpha = R;
    float gamma = 1.0f / alpha;

    {
        // calculate: a = U'*H'
        int   tmpone   = 1;
        float tmpalpha = 1.0f;
        memcpy(a, H_line, sizeof(a[0]) * n); // preload with H_line
        strmm_("L", "U", "T", "U", &n, &tmpone, &tmpalpha, U, &n, a, &n);
    }

    for (int j = 0; j < n; j++)
    {
        b[j] = d[j] * a[j]; // b = D*a = diag(d)*a
    }

    for (int j = 0; j < n; j++)
    {
        float beta = alpha;
        alpha += a[j] * b[j];
        float lambda = -a[j] * gamma;
        gamma        = 1.0f / alpha; // FIXME test if this is possible and return -1 if not
        d[j] *= beta * gamma;
        for (int i = 0; i < j; i++)
        {
            beta                    = MAT_ELEM(U, i, j, n, n);
            MAT_ELEM(U, i, j, n, n) = beta + b[i] * lambda;
            b[i] += b[j] * beta;
        }
    }

    for (int j = 0; j < n; j++)
    {
        x[j] += gamma * dz * b[j];
    }

    return 0;
}

int kalman_udu(float* x, float* U, float* d, const float* z, const float* R, const float* Ht,
               int n, int m, float chi2_threshold, int downweight_outlier)
{
    int retcode = 0;

    for (int i = 0; i < m; i++, Ht += n) /* iterate over each measurement,
                                            goto next line of H after each iteration */
    {
        float Rv = MAT_ELEM(R, i, i, m, m); /// get scalar measurement variance
        float dz = z[i];                    // calculate residual for current scalar measurement
        matmul("N", "N", 1, 1, n, -1.0f, Ht, x, 1.0f, &dz); // dz = z - H(i,:)*x

        // <robust>
        if (chi2_threshold > 0.0f)
        {
            float tmp[KALMAN_MAX_STATE_SIZE];
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

        int status = kalman_udu_scalar(x, U, d, dz, Rv, Ht, n);
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
    int result = cholesky(R, m, 0 /* 0 means: fill upper part with zeros */);
    if (result != 0)
    {
        return -1;
    }
    // L*H_decorr = H
    // (L*H_decorr)' = H'
    // H_decorr'*L' = H' solve for H_decorr
    trisolveright(R /*L*/, Ht, m, n, "T");
    trisolve(R /*L*/, z, m, 1, "N");

    return 0;
}

#if 0

/** @brief UDU' (Thornton) Filter Temporal / Prediction Step
 *
 *  Catherine Thornton's modified weighted Gram-Schmidt orthogonalization
 *  method for the predictor update of the U-D factors of the covariance matrix
 *  of estimation uncertainty in Kalman filtering. Source: [1].
 *
 *  P = U*D*U' = Uin * diag(din) * Uin'
 *
 *  @param[in,out] x   (optional *) state vector with size (n x 1)
 *  @param[in,out] U   unit upper triangular factor (U) of the modified Cholesky
 *                     factors (U-D factors) of the covariance matrix of
 *                     corrected state estimation uncertainty P^{+} (n x n).
 *                     Updated in-place to the modified factors (U-D)
 *                     of the covariance matrix of predicted state
 *                     estimation uncertainty P^{-}, so that
 *                     U*diag(d)*U' = P^{-} after this function.
 *  @param[in,out] d   diagonal factor (d) vector (n x 1) of the U-D factors
 *                     of the covariance matrix of corrected estimation
 *                     uncertainty P^{+}, so that diag(d) = D.
 *                     Updated in-place so that P^{-} = U*diag(d)*U'
 *  @param[in] Phi     state transition matrix (n x n)
 *  @param[in,out] G   process noise distribution matrix (modified, if necessary to
 *                     make the associated process noise covariance diagonal) (n x r)
 *  @param[in] Q       diagonal covariance matrix of process noise
 *                     in the stochastic system model (r x r)
 *
 * (*) Optional, as a non-linear filter will do the prediction of the state vector
 * with a dedicated (non-linear) function.
 * This will basically just predict the state vector with:
 *
 *      x^{-} = Phi*x^{+}
 *      P^{+} = Phi*P^{-}*Phi' + G*Q*G'
 *
 * References:
 *  [1] Grewal, Weill, Andrews. "Global positioning systems, inertial
 *      navigation, and integration". 1st ed. John Wiley & Sons, New York, 2001.
*/
int kalman_udu_predict(float* x, float* U, float* d, const float* Phi, const float* Gin, const float* Q,
                       int n, int r)
{
    assert(n <=
    assert(r <= KALMAN_MAX_STATE_SIZE);

    if (x)
    {
        float tmp[KALMAN_MAX_STATE_SIZE];
        memcpy(tmp, x, sizeof(x[0])*n);
        matmul("N", "N", n, 1, n, 1.0f, Phi, tmp, 0.0f, x);
    }

    float G[KALMAN_MAX_STATE_SIZE*KALMAN_MAX_STATE_SIZE];
    memcpy(G, Gin, sizeof()*n*r);

    // PhiU  = Phi*Uin;  rows of [PhiU,G] are to be orthogonalized
    float PhiU[KALMAN_MAX_STATE_SIZE*KALMAN_MAX_STATE_SIZE];
    float tmpalpha = 1.0f;
    memcpy(PhiU, Phi, sizeof(Phi[0])*n*n);
    strmm_("R", "U", "N", "U", &n, &n, &tmpalpha, Uin, &n, PhiU, &n);

    // U = eye(n)
    mateye(U, n);


    for (int i = n-1; i >= 0; i--)
    {
        float sigma = 0.0f;
        for (int j=0;j<n;j++)
        {
            sigma += MAT_ELEM(PhiU, i, j, n, n) *
                     MAT_ELEM(PhiU, i, j, n, n) * d

        }

    }

/*
x     = Phi*x;
[n,r] = size(Gin);
G     = Gin;       % move to internal array for destructive updates
PhiU  = Phi*Uin;   % rows of [PhiU,G] are to be orthogonalized
U     = eye(n);    % initialize lower triangular part of U
for i=n:-1:1
    sigma = 0;
    for j=1:n
        sigma = sigma + PhiU(i,j)^2 *din(j);
        if (j <= r)
            sigma = sigma + G(i,j)^2 *Q(j,j);
        end;
    end;
    D(i,i) = sigma
    for j=1:i-1
        sigma = 0;
        for k=1:n
            sigma = sigma + PhiU(i,k)*din(k)*PhiU(j,k);
        end;
        for k=1:r
            sigma = sigma + G(i,k)*Q(k,k)*G(j,k);
        end
        U(j,i) = sigma/D(i,i);
        for k=1:n
            PhiU(j,k) = PhiU(j,k) - U(j,i)*PhiU(i,k);
        end
        for k=1:r
            G(j,k) = G(j,k) - U(j,i)*G(i,k);
        end
    end
end
*/

}
#endif

/* @} */
