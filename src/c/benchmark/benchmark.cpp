#include <iostream>
#include <chrono>
#include <random>

#include "benchmark.h"
#include "../navtoolbox.h"
#include "../linalg.h"
#include "../../cpp/navtoolboxeigen.h"
using namespace Eigen;

static int kalman_test1(void)
{
    const float sigma = 0.05f;
    const float bias  = 0.66f;

    std::default_random_engine      generator(42);
    std::normal_distribution<float> distribution(0.0, sigma);

    // <kalman filter>
    const float R[3 * 3]  = { sigma, 0, 0, 0, sigma, 0, 0, 0, sigma };
    const float Ht[4 * 3] = { 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1 };

    float x[4]     = { 0, 0, 0, 0 };
    float P[4 * 4] = { 0.94f, 0, 0, 0, 0, 0.94f, 0, 0, 0, 0, 0.94f, 0, 0, 0, 0, 9.84f };

    float z[3];
    float dz[3];
    // </kalman filter>

    volatile float xr[4] = { 0, 0, 0,
                             0 }; /* make sure the kalman filter results are not optimized away */

    int i;
    for (i = 0; i < 999999; i++)
    {
        z[0] = distribution(generator) + bias;
        z[1] = distribution(generator) + bias;
        z[2] = distribution(generator) + bias;

        // FIXME calc dz = z - H*x
        dz[0] = z[0] - (x[0] + x[3]);
        dz[1] = z[1] - (x[1] + x[3]);
        dz[2] = z[2] - (x[2] + x[3]);

        nav_kalman(x, P, dz, R, Ht, 4, 3);

        xr[0] = x[0];
        xr[1] = x[1];
        xr[2] = x[2];
        xr[3] = x[3];
    }
    printf("Position: %.4f %.4f %.4f\n", (double)xr[0], (double)xr[1], (double)xr[2]);
    printf("Final bias estimation: %.4f (actual: %.4f)\n", (double)xr[3], (double)bias);
    return i;
}

// ----------------------------------------------------------------------------------

static int kalman_test_eigen(void)
{
    const float sigma = 0.05f;
    const float bias  = 0.66f;

    std::default_random_engine      generator(42);
    std::normal_distribution<float> distribution(0.0, sigma);

    // <kalman filter>
    const float R_data[3 * 3]  = { sigma, 0, 0, 0, sigma, 0, 0, 0, sigma };
    const float Ht_data[4 * 3] = { 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1 };

    float x_data[4]     = { 0, 0, 0, 0 };
    float P_data[4 * 4] = { 0.94f, 0, 0, 0, 0, 0.94f, 0, 0, 0, 0, 0.94f, 0, 0, 0, 0, 9.84f };
    // </kalman filter>
    // Convert the arrays to Eigen variables with dynamic sizes
    MatrixXf            R  = MatrixXf::Map(R_data, 3, 3);
    MatrixXf            Ht = MatrixXf::Map(Ht_data, 4, 3);
    MatrixXf            H  = Ht.transpose(); // Transpose to get H
    VectorXf            x  = VectorXf::Map(x_data, 4);
    MatrixXf            P  = MatrixXf::Map(P_data, 4, 4);
    Matrix<float, 3, 1> z;
    Matrix<float, 3, 1> dz;

    volatile float xr[4] = { 0, 0, 0, 0 }; /* volatile: make sure code is not optimized away */

    int i;
    for (i = 0; i < 999999; i++)
    {
        z(0) = distribution(generator) + bias;
        z(1) = distribution(generator) + bias;
        z(2) = distribution(generator) + bias;

        dz = z - H * x;
        nav_kalman_eigen(x, P, dz, R, H);

        xr[0] = x(0);
        xr[1] = x(1);
        xr[2] = x(2);
        xr[3] = x(3);
    }
    printf("Eigen Kalman Position: %.4f %.4f %.4f\n", (double)xr[0], (double)xr[1], (double)xr[2]);
    printf("Eigen Final bias estimation: %.4f (actual: %.4f)\n", (double)xr[3], (double)bias);
    return i;
}
// ----------------------------------------------------------------------------------

void benchmark(void)
{
    {
        auto start = std::chrono::high_resolution_clock::now();

        int loops = kalman_test1();

        auto                          end      = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Duration Test: " << duration.count() << " s. "
                  << "sec/loop: " << (duration.count() / loops) << std::endl;
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        int loops = kalman_test_eigen();

        auto                          end      = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Duration Eigen Test: " << duration.count() << " s. "
                  << "sec/loop: " << (duration.count() / loops) << std::endl;
    }
}
