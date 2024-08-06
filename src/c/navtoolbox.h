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

#ifdef __cplusplus
}
#endif

/* @} */
