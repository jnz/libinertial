#include <iostream>
#include <chrono>
#include <random>

#include "benchmark.h"
#include "../navtoolbox.h"
#include "../linalg.h"
#include "../../cpp/navtoolboxeigen.h"
#include "../kalman_takasu_blasfeo.h"

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

        nav_kalman(x, P, dz, R, Ht, 4, 3, 0.0f, NULL);

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

static int kalman_test1_blasfeo(void)
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

    const int n = 4;
    const int m = 3;

    struct blasfeo_smat sR;
    int R_size = blasfeo_memsize_smat(m, m);
    void *R_mem_align;
    v_zeros_align(&R_mem_align, R_size);
    blasfeo_create_smat(m, m, &sR, R_mem_align);
    blasfeo_pack_smat(m, m, (float*)R, m, &sR, 0, 0);  // convert from column-major to BLASFEO smat

    struct blasfeo_smat sdz;
    int dz_size = blasfeo_memsize_smat(m, 1);
    void *dz_mem_align;
    v_zeros_align(&dz_mem_align, dz_size);
    blasfeo_create_smat(m, 1, &sdz, dz_mem_align);
    blasfeo_pack_smat(m, 1, (float*)dz, m, &sdz, 0, 0);  // convert from column-major to BLASFEO smat

    struct blasfeo_smat sHt;
    int ht_size = blasfeo_memsize_smat(n, m);
    void *ht_mem_align;
    v_zeros_align(&ht_mem_align, ht_size);
    blasfeo_create_smat(n, m, &sHt, ht_mem_align);
    blasfeo_pack_smat(n, m, (float*)Ht, n, &sHt, 0, 0);  // convert from column-major to BLASFEO smat

    struct blasfeo_smat sx;
    int x_size = blasfeo_memsize_smat(n, 1);
    void *x_mem_align;
    v_zeros_align(&x_mem_align, x_size);
    blasfeo_create_smat(n, 1, &sx, x_mem_align);
    blasfeo_pack_smat(n, 1, (float*)x, n, &sx, 0, 0);  // convert from column-major to BLASFEO smat

    struct blasfeo_smat sP;
    int p_size = blasfeo_memsize_smat(n, n);
    void *p_mem_align;
    v_zeros_align(&p_mem_align, p_size);
    blasfeo_create_smat(n, n, &sP, p_mem_align);
    blasfeo_pack_smat(n, n, (float*)P, n, &sP, 0, 0);  // convert from column-major to BLASFEO smat

    struct blasfeo_smat sD;
    int d_size = blasfeo_memsize_smat(n, m);
    void *d_mem_align;
    v_zeros_align(&d_mem_align, d_size);
    blasfeo_create_smat(n, m, &sD, d_mem_align);

    struct blasfeo_smat sL;
    int l_size = blasfeo_memsize_smat(m, m);
    void *l_mem_align;
    v_zeros_align(&l_mem_align, l_size);
    blasfeo_create_smat(m, m, &sL, l_mem_align);

    struct blasfeo_smat sT;
    int t_size = blasfeo_memsize_smat(m, m);
    void *t_mem_align;
    v_zeros_align(&t_mem_align, t_size);
    blasfeo_create_smat(m, m, &sT, t_mem_align);

    struct blasfeo_smat sE;
    int e_size = blasfeo_memsize_smat(n, m);
    void *e_mem_align;
    v_zeros_align(&e_mem_align, e_size);
    blasfeo_create_smat(n, m, &sE, e_mem_align);

    struct blasfeo_smat sxnew;
    int xnew_size = blasfeo_memsize_smat(n, 1);
    void *xnew_mem_align;
    v_zeros_align(&xnew_mem_align, xnew_size);
    blasfeo_create_smat(n, 1, &sxnew, xnew_mem_align);

    struct blasfeo_smat sPnew;
    int pnew_size = blasfeo_memsize_smat(n, n);
    void *pnew_mem_align;
    v_zeros_align(&pnew_mem_align, pnew_size);
    blasfeo_create_smat(n, n, &sPnew, pnew_mem_align);

    int i;
    for (i = 0; i < 999999; i++)
    {
        z[0] = distribution(generator) + bias;
        z[1] = distribution(generator) + bias;
        z[2] = distribution(generator) + bias;

        /* dz[0] = z[0] - (x[0] + x[3]);
           dz[1] = z[1] - (x[1] + x[3]);
           dz[2] = z[2] - (x[2] + x[3]); */

        BLASFEO_SMATEL(&sdz, 0, 0) = z[0] - (BLASFEO_SMATEL(&sx, 0, 0) + BLASFEO_SMATEL(&sx, 3, 0));
        BLASFEO_SMATEL(&sdz, 1, 0) = z[1] - (BLASFEO_SMATEL(&sx, 1, 0) + BLASFEO_SMATEL(&sx, 3, 0));
        BLASFEO_SMATEL(&sdz, 2, 0) = z[2] - (BLASFEO_SMATEL(&sx, 2, 0) + BLASFEO_SMATEL(&sx, 3, 0));

        kalman_takasu_blasfeo(&sx, &sP, &sdz, &sR, &sHt, &sD, &sL, &sT, &sE, &sxnew, &sPnew, 0.0f, NULL);

        blasfeo_sgecp(n, 1, &sxnew, 0, 0, &sx, 0, 0); // sx = sxnew
        blasfeo_sgecp(n, n, &sPnew, 0, 0, &sP, 0, 0); // sP = sPnew

        blasfeo_unpack_smat(4, 1, &sxnew, 0, 0, x, 4);

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
    const int StateDim = 4;
    const int MeasDim  = 3;

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
    Matrix<float, MeasDim, MeasDim>   R  = Matrix<float, MeasDim, MeasDim>::Map(R_data);
    Matrix<float, StateDim, 1>        x  = Matrix<float, StateDim, 1>::Map(x_data);
    Matrix<float, StateDim, StateDim> P  = Matrix<float, StateDim, StateDim>::Map(P_data);
    Matrix<float, StateDim, MeasDim>  Ht = Matrix<float, StateDim, MeasDim>::Map(Ht_data);
    Matrix<float, MeasDim, StateDim>  H  = Ht.transpose();

    Matrix<float, MeasDim, 1> z;
    Matrix<float, MeasDim, 1> dz;

    volatile float xr[4] = { 0, 0, 0, 0 }; /* volatile: make sure code is not optimized away */

    int i;
    for (i = 0; i < 999999; i++)
    {
        z(0) = distribution(generator) + bias;
        z(1) = distribution(generator) + bias;
        z(2) = distribution(generator) + bias;

        dz = z - H * x;
        nav_kalman_eigen<float, StateDim, MeasDim>(x, P, dz, R, H);

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

static int kalman_test_bierman(void)
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
    float U[4 * 4] = { 0 };
    float d[4];
    udu(P, U, d, 4);

    float z[3];
    // </kalman filter>

    volatile float xr[4] = { 0, 0, 0,
                             0 }; /* make sure the kalman filter results are not optimized away */

    int i;
    for (i = 0; i < 999999; i++)
    {
        z[0] = distribution(generator) + bias;
        z[1] = distribution(generator) + bias;
        z[2] = distribution(generator) + bias;

        nav_kalman_udu(x, U, d, z, R, Ht, 4, 3, 0.0, 0);

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

        int loops = kalman_test1_blasfeo();

        auto                          end      = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "BLASFEO Duration Test: " << duration.count() << " s. "
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

    {
        auto start = std::chrono::high_resolution_clock::now();

        int loops = kalman_test_bierman();

        auto                          end      = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Duration Bierman Test: " << duration.count() << " s. "
                  << "sec/loop: " << (duration.count() / loops) << std::endl;
    }
}
