#ifndef WVTICS_KERNEL_H
#define WVTICS_KERNEL_H



float sph_kernel ( const float r, const float h );
float sph_kernel_derivative ( const float r, const float h );
double bias_correction ( const float h );

float sph_kernel_M4 ( const float r, const float h );
float sph_kernel_derivative_M4 ( const float r, const float h );

float sph_kernel_WC2 ( const float r, const float h );
float sph_kernel_derivative_WC2 ( const float r, const float h );
double bias_correction_WC2 ( const float h );

float sph_kernel_WC4 ( const float r, const float h );
float sph_kernel_derivative_WC4 ( const float r, const float h );
double bias_correction_WC4 ( const float h );

float sph_kernel_WC6 ( const float r, const float h );
float sph_kernel_derivative_WC6 ( const float r, const float h );
double bias_correction_WC6 ( const float h );

float sph_kernel_WC8(const float r, const float h);
float sph_kernel_derivative_WC8(const float r, const float h);

float sph_kernel_WC10(const float r, const float h);
float sph_kernel_derivative_WC10(const float r, const float h);

float sph_kernel_WC12 ( const float r, const float h );
float sph_kernel_derivative_WC12 ( const float r, const float h );


#endif //WVTICS_KERNEL_H
