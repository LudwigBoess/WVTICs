#include "kernel.h"
#include "globals.h"

//We omit checking for r < h here and do this in the respective code using the functions
float sph_kernel(const float r, const float h)
{
#ifdef SPH_CUBIC_SPLINE
    return sph_kernel_M4(r, h);
#else
#ifdef SPH_WC2
    return sph_kernel_WC2(r, h);
#else
#ifdef SPH_WC8
    return sph_kernel_WC8(r, h);
#else
#ifdef SPH_WC10
    return sph_kernel_WC10(r, h);
#else
#ifdef SPH_WC12
    return sph_kernel_WC12(r, h);
#else
    // default
    return sph_kernel_WC6(r, h);
#endif // SPH_WC12
#endif // SPH_WC10
#endif // SPH_WC8
#endif // SPH_WC2
#endif // SPH_CUBIC_SPLINE
}

float sph_kernel_derivative(const float r, const float h)
{
#ifdef SPH_CUBIC_SPLINE
    return sph_kernel_derivative_M4(r, h);
#else
#ifdef SPH_WC2
    return sph_kernel_derivative_WC2(r, h);
#else
#ifdef SPH_WC8
    return sph_kernel_derivative_WC8(r, h);
#else
#ifdef SPH_WC10
    return sph_kernel_derivative_WC10(r, h);
#else
#ifdef SPH_WC12
    return sph_kernel_derivative_WC12(r, h);
#else
    // default
    return sph_kernel_derivative_WC6(r, h);
#endif // SPH_WC12
#endif // SPH_WC10
#endif // SPH_WC8
#endif // SPH_WC2
#endif // SPH_CUBIC_SPLINE
}

double bias_correction(const float h)
{
#ifndef SPH_CUBIC_SPLINE
#ifdef SPH_WC2
    return bias_correction_WC2(h);
#else
    return bias_correction_WC6(h);
#endif // SPH_WC2
#else
    return 0.0;
#endif // SPH_CUBIC_SPLINE
}

// See [Steinwandel, Böss, Dolag 2022]

// Wendland C12

float sph_kernel_WC12(const float r, const float h)
{
    const double u = r / h;
    const double t = 1 - u;
    const double t13 = t * t * t * t * t * t * t * t * t * t * t * t * t;

#ifdef TWO_DIM
    double norm = 406980.0 / (2162579 * pi) / p2(h);
#else
    double norm = 445740.0 / (947489 * pi) / p3(h);
#endif //TWO_DIM

    return norm * t13 * t * (77 + 7078 * u + 6755 * u * u + 24500 * u * u * u + 54915 * u * u * u * u + 73206 * u * u * u * u * u + 46189 * u * u * u * u * u * u);
}

float sph_kernel_derivative_WC12(const float r, const float h)
{
    const float u = r / h;
    const double t = 1 - u;
    const double t13 = t * t * t * t * t * t * t * t * t * t * t * t * t;

#ifdef TWO_DIM
    double norm = 406980.0 / (2162579 * pi) / p2(h);
#else
    double norm = 445740.0 / (947489 * pi) / p3(h);
#endif //TWO_DIM

    return norm / h * -20 * t13 * u * (46189 * u * u * u * u * u + 55689 * u * u * u * u + 31122 * u * u * u + 374 * u + 34);
}

// Wendland C10
float sph_kernel_WC10(const float r, const float h)
{
    const double u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 19.0 / (11 * pi) / p2(h);
#else
    double norm = 1995.0 / (512 * pi) / p3(h);
#endif //TWO_DIM

    return norm * t * t * t * t * t * t * t * t * t * t * t * t * (9 + 108 * u + 566 * u * u + 1644 * u * u * u + 2697 * u * u * u * u + 2048 * u * u * u * u * u);
}

float sph_kernel_derivative_WC10(const float r, const float h)
{
    const float u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 19.0 / (11 * pi) / p2(h);
#else
    double norm = 1995.0 / (512 * pi) / p3(h);
#endif //TWO_DIM

    return norm / h * -8 * t * t * t * t * t * t * t * t * t * t * t * u * (4352 * u * u * u * u + 4078 * u * u * u + 1743 * u * u + 374 * u + 34);
}

// Wendland C8

float sph_kernel_WC8(const float r, const float h)
{
    const double u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 8.0 / (3 * pi) / p2(h);
#else
    double norm = 357.0 / (64 * pi) / p3(h);
#endif //TWO_DIM

    return norm * t * t * t * t * t * t * t * t * t * t * (5 + 50 * u + 210 * u * u + 450 * u * u * u + 429 * u * u * u * u);
}

float sph_kernel_derivative_WC8(const float r, const float h)
{
    const float u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 8.0 / (3 * pi) / p2(h);
#else
    double norm = 357.0 / (64 * pi) / p3(h);
#endif //TWO_DIM

    return norm / h * -26 * t * t * t * t * t * t * t * t * t * u * (231 * u * u * u + 159 * u * u + 45 * u + 5);
}

// See [Dehnen & Ali 12]

float sph_kernel_WC6(const float r, const float h)
{
    const double u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 78.0 / (7 * pi) / p2(h);
#else
    double norm = 1365.0 / (64 * pi) / p3(h);
#endif //TWO_DIM

    return norm * t * t * t * t * t * t * t * t * (1 + 8 * u + 25 * u * u + 32 * u * u * u);
}

float sph_kernel_derivative_WC6(const float r, const float h)
{
    const float u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 78.0 / (7 * pi) / p2(h);
#else
    double norm = 1365.0 / (64 * pi) / p3(h);
#endif //TWO_DIM

    return norm / h * -22.0 * t * t * t * t * t * t * t * u * (16 * u * u + 7 * u + 1);
}

double bias_correction_WC6(const float h)
{
#ifdef TWO_DIM
    return 0.0;
#else
    return -0.0116 * pow(DESNNGB * 0.01, -2.236) * Problem.Mpart * sph_kernel_WC6(0, h);
#endif //TWO_DIM
}

float sph_kernel_WC2(const float r, const float h)
{
    const double u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 7.0 / pi / p2(h);
#else
    double norm = 21.0 / (2 * pi) / p3(h);
#endif //TWO_DIM

    return norm * t * t * t * t * (1 + 4 * u);
}

float sph_kernel_derivative_WC2(const float r, const float h)
{
    const float u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 7.0 / pi / p2(h);
#else
    double norm = 21.0 / (2 * pi) / p3(h);
#endif //TWO_DIM

    return norm / h * -20.0 * t * t * t * u;
}

double bias_correction_WC2(const float h)
{
#ifdef TWO_DIM
    return 0.0;
#else
    return -0.0294 * pow(DESNNGB * 0.01, -0.977) * Problem.Mpart * sph_kernel_WC2(0, h);
#endif //TWO_DIM
}

float sph_kernel_M4(const float r, const float h) // cubic spline
{
    const double u = r / h;
    const double t = 1 - u;
    const double v = (u > 0.5 ? 0.0 : 0.5 - u);

#ifdef TWO_DIM
    double norm = 80.0 / (7.0 * pi) / p2(h);
#else
    double norm = 16.0 / pi / p3(h);
#endif //TWO_DIM

    return norm * (p3(t) - 4.0 * p3(v));
}

float sph_kernel_derivative_M4(const float r, const float h)
{
    const double u = r / h;
    const double t = 1 - u;
    const double v = (u > 0.5 ? 0.0 : 0.5 - u);

#ifdef TWO_DIM
    double norm = 80.0 / (7.0 * pi) / p2(h);
#else
    double norm = 16.0 / pi / p3(h);
#endif //TWO_DIM

    return norm / h * (-3.0 * p2(t) + 12.0 * p2(v));
}
