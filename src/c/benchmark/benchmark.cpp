#include <iostream>
#include <chrono>
#include <random>

#include "benchmark.h"
#include "../navtoolbox.h"
#include "../linalg.h"

#if 0
#include "Eigen/Dense"
using namespace Eigen;

/** Kalman Filter Fusion step.  Expects fixed size matrices 3x3.
 * Max. measurements: 3. static function.
 * @param dl Measurement dl vector (dl = measurement - f(x))
 * @param Qll Covariance of measurements
 * @param A Design/Jacobi matrix
 * @param dx Output state correction vector
 * @param Qxx Output updated state covariance matrix */
void KalmanFusion3x3(const Matrix<float, 3, 1>& dz,
                     const Matrix<float, 3, 3>& R,
                     const Matrix<float, 3, 4>& A,
                     Matrix<float, 4, 1>& x,
                     Matrix<float, 4, 4>& P)
{
    // Using T. Takasu's Kalman filter Cholesky-based method:
    //
    // D = Qxx*A'
    // S = A*D + Qll
    // U = chol(S) (with U'*U = S)
    // E = D*inv(U)
    // K = E*inv(U)'
    // x = x + K*dl
    // P = P - E*E'

    Matrix<float, 4, 3> K;
    Matrix<float, 4, 3> D;
    Matrix<float, 3, 3> S;
    Matrix<float, 3, 3> Uinv;
    Matrix<float, 4, 3> E;

    D.noalias() = P.selfadjointView<Eigen::Upper>()*(A.transpose());
    S.triangularView<Eigen::Upper>() = A*D + R;

    Uinv.setIdentity();
    S.selfadjointView<Eigen::Upper>().llt().matrixU().solveInPlace(Uinv); // calculate Uinv

    E.noalias() = D * Uinv;
    K.noalias() = E * Uinv.transpose();

    P.selfadjointView<Eigen::Upper>().rankUpdate(E, -1); // P = P - E*E'
    x += K*dz;
}
#endif

static void kalman_test1(void)
{
    const float sigma = 0.25f;
    const float bias = 0.66f;

    std::default_random_engine generator;
    std::normal_distribution<float> distribution(0.0, sigma);

    // <kalman filter>
    const float R[3 * 3]  = { sigma, 0, 0, 0, sigma, 0, 0, 0, sigma };
    const float Ht[4 * 3] = { 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1 };

    float x[4]     = { 0, 0, 0, 0 };
    float P[4 * 4] = { 0.24f, 0, 0, 0, 0, 0.24f, 0, 0, 0, 0, 0.24f, 0, 0, 0, 0, 0.84f };

    float z[3];
    float dz[3];
    // </kalman filter>

    volatile float xr[4] = { 0, 0, 0, 0 }; /* volatile: make sure code is not optimized away */

    for (int i=0;i<999999;i++)
    {
        z[0] = distribution(generator) + bias;
        z[1] = distribution(generator) + bias;
        z[2] = distribution(generator) + bias;

        // FIXME calc dz
        dz[0] = z[0];
        dz[1] = z[1];
        dz[2] = z[2];

        nav_kalman(x, P, dz, R, Ht, 4, 3);

        xr[0] = x[0];
        xr[1] = x[1];
        xr[2] = x[2];
        xr[3] = x[3];
    }

    (void)xr;
}

static void kalman_test2(void)
{

}

void benchmark(void)
{
    auto start = std::chrono::high_resolution_clock::now();

    kalman_test1();
    kalman_test2();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    std::cout << "Duration: " << duration.count() << " s" << std::endl;
}
