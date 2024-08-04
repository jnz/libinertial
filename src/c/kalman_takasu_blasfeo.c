/** @file kalman_takasu_blasfeo.c
 * @author Jan Zwiener (jan@zwiener.org)
 *
 * @brief BLASFEO High Performance Version of the Takasu filter
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <math.h>
#include <assert.h>
#include <string.h> /* memcpy */
#include <stdio.h> /* remove me */

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "kalman_takasu_blasfeo.h"

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

/******************************************************************************
 * FUNCTION BODIES
 ******************************************************************************/

int kalman_takasu_blasfeo(struct blasfeo_smat* x, struct blasfeo_smat* P,
                          const struct blasfeo_smat* dz, const struct blasfeo_smat* R,
                          const struct blasfeo_smat* Ht,
                          struct blasfeo_smat* D,
                          struct blasfeo_smat* L,
                          struct blasfeo_smat* T,
                          struct blasfeo_smat* E,
                          struct blasfeo_smat* xnew,
                          struct blasfeo_smat* Pnew,
                          float chi2_threshold, float* chi2)
{

    /*  (1) D = P * H'              symm           |   Matrix dimensions:
     *  (2) S = H * D + R           gemm           |   D = n x m
     *  (3) L = chol(S) (L*L'=S)    potrf          |   S = m x m
     *  (4) E = D * L^-T            trsm           |   L = m x m
     *  (5) P = P - E*E'            syrk           |   E = n x m
     *  (6) K = E * L^-1            trsm           |   K = n x m
     *  (7) x = x + K*dz            gemm               T = m x m
     *
     *  n = state variables
     *  m = measurements */

    const int n = x->m;
    const int m = dz->m;

    // (1) D = P * H'
    blasfeo_sgemm_nn(n, m, n, 1.0f, P, 0, 0, (struct blasfeo_smat*)Ht, 0, 0, 0.0f, D, 0, 0, D, 0, 0);
    /*
     *printf("D = \n");
     *blasfeo_print_smat(n, m, D, 0, 0);
     *printf("--- \n");
     */
    // (2) T = H*D + R
    blasfeo_sgemm_tn(m, m, n, 1.0f, (struct blasfeo_smat*)Ht, 0, 0, D, 0, 0, 1.0f, (struct blasfeo_smat*)R, 0, 0, T, 0, 0);
    /*
     *printf("S = \n");
     *blasfeo_print_smat(m, m, T, 0, 0);
     *printf("--- \n");
     */
    // (3) L = chol(T) = chol(H*P*H' + R)
    blasfeo_spotrf_l(m, T, 0, 0, L, 0, 0);
    /*
     *printf("L = \n");
     *blasfeo_print_smat(m, m, L, 0, 0);
     *printf("--- \n");
     */
    // (4) given L' and D, solve E*L' = D, for E
    blasfeo_strsm_rltn(n, m, 1.0f, L, 0, 0, D, 0, 0, E, 0, 0);
    /*
     *printf("E = \n");
     *blasfeo_print_smat(n, m, E, 0, 0);
     *printf("--- \n");
     */
    // (5) Pnew = P - E*E'
    blasfeo_ssyrk_un(n, m, -1.0f, E, 0, 0, E, 0, 0, 1.0f, P, 0, 0, Pnew, 0, 0);
    /*
     *printf("Pnew = \n");
     *blasfeo_print_smat(n, n, Pnew, 0, 0);
     *printf("--- \n");
     */
    // (6) solve K*L = E, for K
    blasfeo_strsm_rlnn(n, m, 1.0f, L, 0, 0, E, 0, 0, D/*=K*/, 0, 0);
    /*
     *printf("K = \n");
     *blasfeo_print_smat(n, m, D, 0, 0);
     *printf("--- \n");
     */
    // (7) x = x + K * dz (K is stored in D)
    blasfeo_sgemm_nn(n, 1, m, 1.0f, D, 0, 0, (struct blasfeo_smat*)dz, 0, 0, 1.0f, x, 0, 0, xnew, 0, 0);
    /*
     *printf("xnew = \n");
     *blasfeo_print_smat(n, 1, xnew, 0, 0);
     *printf("--- \n");
     */

    // Copy upper triangular part of P to lower triangular part of P
    for (int k = 0; k < n; k++) {
        for (int l = k + 1; l < n; l++) {
            BLASFEO_SMATEL(Pnew, l, k) = BLASFEO_SMATEL(Pnew, k, l);
        }
    }

    return 0;
}

/* @} */
