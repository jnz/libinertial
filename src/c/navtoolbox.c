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

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

#include "navtoolbox.h"
#include "linalg.h"

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

/* @} */
