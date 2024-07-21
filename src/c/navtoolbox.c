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

#if 0
void nav_filter(void)
{
    /* matmul, matinv */
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
}
#endif

/* @} */