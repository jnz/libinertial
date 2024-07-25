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
 * Note (!): only the upper triangular part of P is referenced and updated.
 * @return 0 on success, -1 on error.
 */
int nav_kalman_eigen(
    Eigen::Matrix<float, Eigen::Dynamic, 1>&                    x,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>&       P,
    const Eigen::Matrix<float, Eigen::Dynamic, 1>&              dz,
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& R,
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& H);

/* @} */
