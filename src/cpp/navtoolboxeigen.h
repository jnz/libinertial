/** @file navtoolboxeigen.h
 * libinertial, Jan Zwiener (jan@zwiener.org)
 *
 * @brief Navigation Toolbox Helper Functions for Eigen C++ Mathlib
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <Eigen/Dense>

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "../c/navtoolbox.h"

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
 * FUNCTION PROTOTYPES
 ******************************************************************************/

/** @brief Kalman Filter Routine.
 *
 * @param[in,out] x System state (dimension n)
 * @param[in,out] P Upper triangular Covariance matrix of state estimation uncertainty
 * @param[in] dz Difference between measurement and expected measurement: z - H*x (dimension m)
 * @param[in] R Covariance matrix of measurement uncertainty (dimension m x m)
 * @param[in] H Measurement sensitivity matrix (m x n)
 *
 * n = StateDim, m = MeasDim
 *
 * Note (!): only the upper triangular part of P is referenced and updated.
 * @return 0 on success, -1 on error.
 */

template <typename Scalar, int StateDim, int MeasDim>
int nav_kalman_eigen(
    Eigen::Matrix<Scalar, StateDim, 1>& x,
    Eigen::Matrix<Scalar, StateDim, StateDim>& P,
    const Eigen::Matrix<Scalar, MeasDim, 1>& dz,
    const Eigen::Matrix<Scalar, MeasDim, MeasDim>& R,
    const Eigen::Matrix<Scalar, MeasDim, StateDim>& H)
{
    // Takasu formulation:
    Eigen::Matrix<Scalar, StateDim, MeasDim> D;
    Eigen::Matrix<Scalar, MeasDim, MeasDim> L = R; // L is used as a temp matrix and preloaded with R

    // (1) D = P * H'
    D.noalias() = P.template selfadjointView<Eigen::Upper>() * H.transpose();
    // (2) L = H * D + R
    L.noalias() += H * D;
    // (3) L = chol(L)
    Eigen::LLT<Eigen::Matrix<Scalar, MeasDim, MeasDim>> lltOfL(L);
    if (lltOfL.info() != Eigen::Success)
    {
        return -1; // Cholesky decomposition failed
    }
    L = lltOfL.matrixL();

    // (4) E = D * (L')^-1
    Eigen::Matrix<Scalar, StateDim, MeasDim> E = L.template triangularView<Eigen::Lower>().solve(D.transpose()).transpose();

    // (5) P = P - E * E'
    P.template selfadjointView<Eigen::Upper>().rankUpdate(E, -1);

    // (6) K = E * L^-1
    Eigen::Matrix<Scalar, StateDim, MeasDim> K =
        L.transpose().template triangularView<Eigen::Upper>().solve(E.transpose()).transpose();

    x.noalias() += K * dz;

    return 0;
}

/* @} */
