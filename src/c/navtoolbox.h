/** @file navtoolbox.h
 * libinertial, Jan Zwiener (jan@zwiener.org)
 *
 * @brief Navigation Toolbox Helper Functions
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

#define PI_FLOAT (3.141592653589793f)
#define RAD2DEG(x) ((x) * (180.0f / PI_FLOAT))
#define DEG2RAD(x) ((x) * (PI_FLOAT / 180.0f))
#define CLIGHT (299792458.0)    /* speed of light (m/s) */
#define OMGE (7.2921151467E-5f) /* Earth rotation rate 15deg/h */
#define GRAVITY (9.81f)         /* Gravity */

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

#ifdef __cplusplus
extern "C"
{
#endif

    /** @brief Calculate an approximate orientation from accelerometer data,
     * assuming that the accelerometer measurement is mainly gravity.
     *
     * @param[in] f Specific force measurement x,y,z component (m/s^2)
     * @param[out] roll_rad Output roll angle (rad)
     * @param[out] pitch_rad Output pitch angle (rad) */
    void nav_roll_pitch_from_accelerometer(const float f[3], float* roll_rad, float* pitch_rad);

    /** @brief nav_matrix_body2nav Calculate a matrix R that transforms from
     * the body-frame (b) to the navigation-frame (n): R^n_b.
     * @param[in] roll_rad Roll angle in (rad)
     * @param[in] pitch_rad Pitch angle in (rad)
     * @param[in] yaw_rad Yaw angle in (rad)
     * @param[out] R_output Output 3x3 matrix in column-major format */
    void nav_matrix_body2nav(const float roll_rad, const float pitch_rad, const float yaw_rad,
                             float R_output[9]);

    /** @brief Kalman Filter Routine.
     *
     * @param[in,out] x System state (n x 1)
     * @param[in,out] P Upper triangular Covariance matrix of state estimation uncertainty (n x n)
     * @param[in] dz Measurement residual vector: measurement vs. expected measurement: z - H*x (m x
     * 1)
     * @param[in] R Full covariance matrix of measurement uncertainty (m x m)
     * @param[in] Ht Transposed (!) measurement sensitivity matrix (n x m) (H would be m x n)
     * @param[in] n Number of state variables
     * @param[in] m Number of measurements
     *
     * Note (!): only the upper triangular part of P is referenced and updated.
     *
     * @return 0 on success, -1 on error.
     */
    int nav_kalman(float* x, float* P, const float* dz, const float* R, const float* Ht, int n,
                   int m);

    /** @brief (Robust) Square Root Kalman Filter routine for linear system.
     *
     * @param[in,out] x System state (n x 1)
     * @param[in,out] U Unit upper triangular factor of covariance matrix of a priori state
     * uncertainty (n x n)
     * @param[in,out] d Unit upper triangular factor of covariance matrix of a priori state
     * uncertainty (n x 1)
     * @param[in] z Measurement vector: z = H*x (m x 1)
     * @param[in] R Full covariance matrix of measurement uncertainty (m x m)
     * @param[in] Ht Transposed (!) measurement sensitivity matrix (n x m) (H would be m x n)
     * @param[in] n Number of state variables
     * @param[in] m Number of measurements
     * @param[in] chi2_threshold Scalar threshold for outlier classification. Set to 0.0f to
     * disable.
     * @param[in] downweight_outlier If set to 0, measurements classified as outliers are skipped.
     *
     * @return 0 on success, -1 on error.
     */
    int nav_kalman_udu(float* x, float* U, float* d, const float* z, const float* R,
                       const float* Ht, int n, int m, float chi2_threshold, int downweight_outlier);

    /** @brief Square Root Kalman Filter Routine for scalar measurement.
     *
     * @param[in,out] x System state (n x 1)
     * @param[in,out] U Unit upper triangular factor of covariance matrix of a priori state
     * uncertainty (n x n)
     * @param[in,out] d Unit upper triangular factor of covariance matrix of a priori state
     * uncertainty (n x 1)
     * @param[in] dz Scalar measurement residual: dz = z - H*x (1 x 1)
     * @param[in] R Scalar covariance of measurement uncertainty (1 x 1)
     * @param[in] H_line Row of measurement sensitivity matrix (n x 1)
     * @param[in] n Number of state variables
     */
    int nav_kalman_udu_scalar(float* x, float* U, float* d, const float dz, const float R,
                              const float* H_line, int n);

    /** @brief Decorrelate measurements. For a given covariance matrix R of correlated measurements,
     * calculate a vector of decorrelated measurements (and the matching H-matrix) so that
     * the new covariance is an identity matrix.
     *
     * Note: If R only has diagonal elements, a call to this function does not help.
     *
     * @param[in,out] z Vector of correlated measurements (m x 1).
     * @param[in,out] Ht Transposed measurement sensitivity matrix / design matrix so that z = Ht'*x
     * (n x m).
     * @param[in,out] R Measurement covariance matrix, replaced in-place by chol(R) (m x m).
     *                  So the input is R, overwritten with L such that L*L'=R.
     * @param[in] n Number of columns in H (for a Kalman filter: length of state vector x).
     * @param[in] m Number of measurements in z.
     *
     * If Ht and R do not change, subsequently only measurements z need to be decorrelated.
     * As the input R is replaced by L (such that L*L' = R), L can be reused to
     * decorrelate further measurements:
     *
     *     trisolve(L, z, m, 1, "N");
     *
     * @return 0 if successful, if -1 state of z and H is not guaranteed to be
     *           consistent and must be discarded.
     */
    int decorrelate(float* z, float* Ht, float* R, int n, int m);

#ifdef __cplusplus
}
#endif

/* @} */
