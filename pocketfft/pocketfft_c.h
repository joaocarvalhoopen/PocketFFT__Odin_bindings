/* Plain-C interface to pocketfft (subset: float and double).
 * Each function returns 0 on success, non-zero on failure (exception caught). */

#ifndef POCKETFFT_C_H
#define POCKETFFT_C_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ndim       : number of dimensions (length of shape, stride_in, stride_out)
 * shape      : array of dimension lengths (size_t[ndim])
 * stride_in  : input strides in BYTES (ptrdiff_t[ndim])
 * stride_out : output strides in BYTES (ptrdiff_t[ndim])
 * naxes      : number of axes to transform (length of axes)
 * axes       : axes to transform (size_t[naxes])
 * forward    : non-zero for forward transform
 * data_in    : input buffer pointer
 * data_out   : output buffer pointer
 * fct        : scaling factor applied to result
 * nthreads   : 0 = use all cores, otherwise number of threads
 */

/* complex-to-complex */
int pocketfft_c2c_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, void *data_out, float fct, size_t nthreads);

int pocketfft_c2c_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, void *data_out, double fct, size_t nthreads);

/* real-to-complex (single axis) */
int pocketfft_r2c_f32(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const float *data_in, void *data_out, float fct, size_t nthreads);

int pocketfft_r2c_f64(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const double *data_in, void *data_out, double fct, size_t nthreads);

/* real-to-complex (multi axis) */
int pocketfft_r2c_axes_f32(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const float *data_in, void *data_out, float fct, size_t nthreads);

int pocketfft_r2c_axes_f64(size_t ndim, const size_t *shape_in,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const double *data_in, void *data_out, double fct, size_t nthreads);

/* complex-to-real (single axis) */
int pocketfft_c2r_f32(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const void *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_c2r_f64(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t axis, int forward,
    const void *data_in, double *data_out, double fct, size_t nthreads);

/* complex-to-real (multi axis) */
int pocketfft_c2r_axes_f32(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_c2r_axes_f64(size_t ndim, const size_t *shape_out,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes, int forward,
    const void *data_in, double *data_out, double fct, size_t nthreads);

/* r2r FFTPACK-style half-complex */
int pocketfft_r2r_fftpack_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int real2hermitian, int forward,
    const float *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_r2r_fftpack_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int real2hermitian, int forward,
    const double *data_in, double *data_out, double fct, size_t nthreads);

/* separable Hartley */
int pocketfft_r2r_separable_hartley_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_r2r_separable_hartley_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads);

/* genuine Hartley */
int pocketfft_r2r_genuine_hartley_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_r2r_genuine_hartley_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads);

/* separable FHT */
int pocketfft_r2r_separable_fht_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_r2r_separable_fht_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads);

/* genuine FHT */
int pocketfft_r2r_genuine_fht_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const float *data_in, float *data_out, float fct, size_t nthreads);

int pocketfft_r2r_genuine_fht_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    const double *data_in, double *data_out, double fct, size_t nthreads);

/* discrete cosine transform (type 1..4) */
int pocketfft_dct_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const float *data_in, float *data_out,
    float fct, int ortho, size_t nthreads);

int pocketfft_dct_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const double *data_in, double *data_out,
    double fct, int ortho, size_t nthreads);

/* discrete sine transform (type 1..4) */
int pocketfft_dst_f32(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const float *data_in, float *data_out,
    float fct, int ortho, size_t nthreads);

int pocketfft_dst_f64(size_t ndim, const size_t *shape,
    const ptrdiff_t *stride_in, const ptrdiff_t *stride_out,
    size_t naxes, const size_t *axes,
    int type, const double *data_in, double *data_out,
    double fct, int ortho, size_t nthreads);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* POCKETFFT_C_H */
