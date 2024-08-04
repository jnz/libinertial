/** @file kalman_takasu_blasfeo.h
 * @author Jan Zwiener (jan@zwiener.org)
 *
 * @brief BLASFEO Kalman Takasu Implementation
 *
 * @{ */

/******************************************************************************
 * SYSTEM INCLUDE FILES
 ******************************************************************************/

#include <blasfeo.h>

/******************************************************************************
 * PROJECT INCLUDE FILES
 ******************************************************************************/

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

#ifdef __cplusplus
extern "C"
{
#endif

int kalman_takasu_blasfeo(struct blasfeo_smat* x, struct blasfeo_smat* P,
                          const struct blasfeo_smat* dz, const struct blasfeo_smat* R,
                          const struct blasfeo_smat* Ht,
                          struct blasfeo_smat* D,
                          struct blasfeo_smat* L,
                          struct blasfeo_smat* T,
                          struct blasfeo_smat* E,
                          struct blasfeo_smat* xnew,
                          struct blasfeo_smat* Pnew,
                          float chi2_threshold, float* chi2);

#ifdef __cplusplus
}
#endif

/* @} */
