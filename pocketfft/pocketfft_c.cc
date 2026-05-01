/* C wrapper implementation around pocketfft ( header-only C++ template lib ).
 * Instantiates float and double versions and exposes a stable C ABI. */

#include "pocketfft_hdronly.h"
#include "pocketfft_c.h"

#include <complex>
#include <vector>
#include <exception>

using std::complex;
using std::vector;
using pocketfft::shape_t;
using pocketfft::stride_t;

namespace {

inline shape_t to_shape(size_t n, const size_t *p) {
    return shape_t(p, p + n);
}
inline stride_t to_stride(size_t n, const ptrdiff_t *p) {
    return stride_t(p, p + n);
}

} /* anonymous namespace */

#define WRAP_BEGIN try {
#define WRAP_END } catch (const std::exception &) { return 1; } catch (...) { return 2; } return 0;

extern "C" {

/* ---- c2c ---- */
int pocketfft_c2c_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, void *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::c2c(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes), forward != 0,
        static_cast<const complex<float>*>(data_in),
        static_cast<complex<float>*>(data_out), fct, nthreads);
WRAP_END
}

int pocketfft_c2c_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, void *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::c2c(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes), forward != 0,
        static_cast<const complex<double>*>(data_in),
        static_cast<complex<double>*>(data_out), fct, nthreads);
WRAP_END
}

/* ---- r2c (single axis) ---- */
int pocketfft_r2c_f32(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const float *data_in, void *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2c(to_shape(ndim, shape_in), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), axis, forward != 0,
        data_in, static_cast<complex<float>*>(data_out), fct, nthreads);
WRAP_END
}

int pocketfft_r2c_f64(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const double *data_in, void *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2c(to_shape(ndim, shape_in), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), axis, forward != 0,
        data_in, static_cast<complex<double>*>(data_out), fct, nthreads);
WRAP_END
}

/* ---- r2c (multi axis) ---- */
int pocketfft_r2c_axes_f32(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const float *data_in, void *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2c(to_shape(ndim, shape_in), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes), forward != 0,
        data_in, static_cast<complex<float>*>(data_out), fct, nthreads);
WRAP_END
}

int pocketfft_r2c_axes_f64(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const double *data_in, void *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2c(to_shape(ndim, shape_in), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes), forward != 0,
        data_in, static_cast<complex<double>*>(data_out), fct, nthreads);
WRAP_END
}

/* ---- c2r (single axis) ---- */
int pocketfft_c2r_f32(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const void *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::c2r(to_shape(ndim, shape_out), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), axis, forward != 0,
        static_cast<const complex<float>*>(data_in), data_out, fct, nthreads);
WRAP_END
}

int pocketfft_c2r_f64(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const void *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::c2r(to_shape(ndim, shape_out), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), axis, forward != 0,
        static_cast<const complex<double>*>(data_in), data_out, fct, nthreads);
WRAP_END
}

/* ---- c2r (multi axis) ---- */
int pocketfft_c2r_axes_f32(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::c2r(to_shape(ndim, shape_out), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes), forward != 0,
        static_cast<const complex<float>*>(data_in), data_out, fct, nthreads);
WRAP_END
}

int pocketfft_c2r_axes_f64(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::c2r(to_shape(ndim, shape_out), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes), forward != 0,
        static_cast<const complex<double>*>(data_in), data_out, fct, nthreads);
WRAP_END
}

/* ---- r2r_fftpack ---- */
int pocketfft_r2r_fftpack_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int real2hermitian, int forward,
    const float *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_fftpack(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes),
        real2hermitian != 0, forward != 0,
        data_in, data_out, fct, nthreads);
WRAP_END
}

int pocketfft_r2r_fftpack_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int real2hermitian, int forward,
    const double *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_fftpack(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes),
        real2hermitian != 0, forward != 0,
        data_in, data_out, fct, nthreads);
WRAP_END
}

/* ---- r2r_separable_hartley ---- */
int pocketfft_r2r_separable_hartley_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_separable_hartley(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

int pocketfft_r2r_separable_hartley_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_separable_hartley(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

/* ---- r2r_genuine_hartley ---- */
int pocketfft_r2r_genuine_hartley_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_genuine_hartley(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

int pocketfft_r2r_genuine_hartley_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_genuine_hartley(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

/* ---- r2r_separable_fht ---- */
int pocketfft_r2r_separable_fht_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_separable_fht(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

int pocketfft_r2r_separable_fht_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_separable_fht(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

/* ---- r2r_genuine_fht ---- */
int pocketfft_r2r_genuine_fht_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_genuine_fht(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

int pocketfft_r2r_genuine_fht_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::r2r_genuine_fht(to_shape(ndim, shape),
        to_stride(ndim, stride_in), to_stride(ndim, stride_out),
        to_shape(naxes, axes), data_in, data_out, fct, nthreads);
WRAP_END
}

/* ---- dct ---- */
int pocketfft_dct_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const float *data_in, float *data_out,
    float fct, int ortho, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::dct(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes),
        type, data_in, data_out, fct, ortho != 0, nthreads);
WRAP_END
}

int pocketfft_dct_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const double *data_in, double *data_out,
    double fct, int ortho, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::dct(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes),
        type, data_in, data_out, fct, ortho != 0, nthreads);
WRAP_END
}

/* ---- dst ---- */
int pocketfft_dst_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const float *data_in, float *data_out,
    float fct, int ortho, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::dst(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes),
        type, data_in, data_out, fct, ortho != 0, nthreads);
WRAP_END
}

int pocketfft_dst_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const double *data_in, double *data_out,
    double fct, int ortho, size_t nthreads)
{
WRAP_BEGIN
    pocketfft::dst(to_shape(ndim, shape), to_stride(ndim, stride_in),
        to_stride(ndim, stride_out), to_shape(naxes, axes),
        type, data_in, data_out, fct, ortho != 0, nthreads);
WRAP_END
}

} /* extern "C" */
