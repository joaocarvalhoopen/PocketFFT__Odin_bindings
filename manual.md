# PocketFFT Odin bindings — Manual

A reference for every public symbol exposed by `pocketfft.odin`.

The package wraps the upstream pocketfft C++ library via a thin C ABI
(`libpocketfft_c.a`). Two API layers are available:

1. **Raw foreign procedures** — minimal-overhead, 1-to-1 with the C ABI; take
   bare `[^]` pointers and `c.size_t`/`c.ptrdiff_t`.
2. **Ergonomic wrappers** — accept `[]T` slices, validate buffer sizes, and
   return a typed `Error`. Recommended for application code.

Both layers share the same conventions:

* Strides are in **bytes**, not elements.
* Shapes describe the **logical multidimensional** array layout (row-major
  by default — see `strides_row_major`).
* Multidimensional data is stored as a flat `[]T` slice. Element at index
  `(i0, i1, …, iN-1)` lives at offset `Σ ik * stride[k] / sizeof(T)` for
  contiguous arrays, or at byte offset `Σ ik * stride[k]` for arbitrary
  strides.
* For `r2c`/`c2r` transforms the complex side has length `s/2 + 1` along the
  transformed axis, where `s` is the real axis length.
* `forward = true` ⇒ minus sign in the complex exponent;
  `forward = false` ⇒ plus sign. None of the transforms apply normalization
  automatically; if you need an inverse transform that returns the original
  data, scale by `1/N` (or whatever convention you prefer).

All wrappers come in two overloads (`f32` and `f64`); pass `f32`/`complex64`
or `f64`/`complex128` slices and Odin selects the right one. Long double is
not supported (Odin has no native type).

---

## Constants

| Name       | Value | Description                      |
| ---------- | ----- | -------------------------------- |
| `FORWARD`  | `true`  | Forward transform (e<sup>-iωt</sup>). |
| `BACKWARD` | `false` | Backward transform (e<sup>+iωt</sup>). |

## Error type

```odin
Error :: enum {
    None             = 0,  // success
    Cxx_Exception    = 1,  // pocketfft raised a std::exception (e.g. invalid args)
    Unknown          = 2,  // unknown C++ exception
    Invalid_Buffer   = 3,  // input/output slice too small for the given shape
    Empty_Shape      = 4,  // shape parameter has zero dimensions
}
```

## Helpers

### `strides_row_major(shape, elem_size, allocator) -> []int`

Compute byte-strides for a contiguous row-major array.

| Parameter   | Type             | Description                                                |
| ----------- | ---------------- | ---------------------------------------------------------- |
| `shape`     | `[]uint`         | Logical array dimensions (length = number of axes).        |
| `elem_size` | `int`            | Size in bytes of one element (e.g. `size_of(complex128)`). |
| `allocator` | `mem.Allocator`  | Optional. Defaults to `context.allocator`.                 |
| **returns** | `[]int`          | Stride per axis, in bytes. Caller must `delete` it.        |

Example for a `[8, 16] f64` row-major array:

```odin
str := strides_row_major({8, 16}, size_of(f64))
// str == [128, 8]   (16 elements × 8 bytes per row, then 8 bytes per col)
```

---

# 1. Complex-to-complex transforms

## `c2c_simple(shape, forward, data_in, data_out, fct, nthreads) -> Error`

Easy entry point: contiguous, row-major, transform along **all** axes.

| Parameter   | Type                            | Default | Description |
| ----------- | ------------------------------- | ------- | ----------- |
| `shape`     | `[]uint`                        | —       | Logical dimensions of the array. |
| `forward`   | `bool`                          | —       | `FORWARD` or `BACKWARD`. |
| `data_in`   | `[]complex64` or `[]complex128` | —       | Input array (length = product of `shape`). |
| `data_out`  | same complex type as `data_in`  | —       | Output array (length = product of `shape`). May alias `data_in`. |
| `fct`       | `f32` / `f64`                   | `1`     | Scalar applied to every output element. |
| `nthreads`  | `uint`                          | `1`     | Threads to use. `0` ⇒ all logical CPUs. |
| **returns** | `Error`                         |         | `.None` on success. |

## `c2c_f32_proc(shape, axes, stride_in, stride_out, forward, data_in, data_out, fct, nthreads) -> Error`
## `c2c_f64_proc(shape, axes, stride_in, stride_out, forward, data_in, data_out, fct, nthreads) -> Error`

Full-control variant for non-contiguous data, partial axes, custom strides.

| Parameter    | Type                          | Description |
| ------------ | ----------------------------- | ----------- |
| `shape`      | `[]uint`                      | Array dimensions. |
| `axes`       | `[]uint`                      | Indices of axes to transform (subset of `0..len(shape)-1`). |
| `stride_in`  | `[]int`                       | Input byte strides per axis. |
| `stride_out` | `[]int`                       | Output byte strides per axis. |
| `forward`    | `bool`                        | Direction. |
| `data_in`    | `[]complex64`/`[]complex128`  | Input buffer. |
| `data_out`   | matching complex type         | Output buffer. |
| `fct`        | `f32`/`f64`                   | Output scaling factor. |
| `nthreads`   | `uint`                        | Threads. |

## Raw foreign procs: `c2c_f32`, `c2c_f64`

Direct C ABI; take `[^]c.size_t`, `[^]c.ptrdiff_t`, `rawptr` and an integer
`forward` (0/1). Return `c.int` (0 on success). Use only when you need to
avoid the slice-validation overhead.

---

# 2. Real-to-complex transforms

## `r2c_simple(shape, forward, data_in, data_out, fct, nthreads) -> Error`

Contiguous r2c. Transforms along the **last** axis only. The output array
must be sized for `shape[0] * … * shape[N-2] * (shape[N-1]/2 + 1)` complex
elements.

| Parameter   | Type                                | Default | Description |
| ----------- | ----------------------------------- | ------- | ----------- |
| `shape`     | `[]uint`                            | —       | Shape of the **real** input. |
| `forward`   | `bool`                              | —       | `FORWARD` for the standard real FFT. |
| `data_in`   | `[]f32` / `[]f64`                   | —       | Real input. |
| `data_out`  | `[]complex64` / `[]complex128`      | —       | Complex output (Hermitian half). |
| `fct`       | `f32` / `f64`                       | `1`     | Output scaling. |
| `nthreads`  | `uint`                              | `1`     | Threads. |

## `r2c_single_f32 / r2c_single_f64(shape_in, stride_in, stride_out, axis, forward, data_in, data_out, fct, nthreads) -> Error`

Single-axis variant with explicit strides.

| Parameter    | Type   | Description |
| ------------ | ------ | ----------- |
| `shape_in`   | `[]uint` | Shape of the real input array. |
| `stride_in`  | `[]int`  | Byte strides of the real input. |
| `stride_out` | `[]int`  | Byte strides of the complex output. |
| `axis`       | `uint`   | Axis along which the r2c FFT is performed. |
| `forward`    | `bool`   | Direction (typically `FORWARD`). |
| `data_in`    | `[]f32`/`[]f64`               | Real input. |
| `data_out`   | `[]complex64`/`[]complex128`  | Complex output. |
| `fct`        | `f32`/`f64` | Output scaling. |
| `nthreads`   | `uint`      | Threads. |

## Raw foreign procs

| Name              | Notes |
| ----------------- | ----- |
| `r2c_f32`/`r2c_f64`             | Single-axis raw C ABI. |
| `r2c_axes_f32`/`r2c_axes_f64`   | Multi-axis raw C ABI. The last entry of `axes` is the r2c axis; remaining axes get a follow-up c2c. |

---

# 3. Complex-to-real transforms

## `c2r_simple(shape_out, forward, data_in, data_out, fct, nthreads) -> Error`

Contiguous c2r. Transforms along the **last** axis. Input length must be
`shape[0] * … * shape[N-2] * (shape[N-1]/2 + 1)` complex elements.

| Parameter   | Type                                | Default | Description |
| ----------- | ----------------------------------- | ------- | ----------- |
| `shape_out` | `[]uint`                            | —       | Shape of the **real** output. |
| `forward`   | `bool`                              | —       | Usually `BACKWARD` for the standard inverse real FFT. |
| `data_in`   | `[]complex64` / `[]complex128`      | —       | Hermitian-half complex input. |
| `data_out`  | `[]f32` / `[]f64`                   | —       | Real output. |
| `fct`       | `f32` / `f64`                       | `1`     | Output scaling (commonly `1/N`). |
| `nthreads`  | `uint`                              | `1`     | Threads. |

## `c2r_single_f32 / c2r_single_f64(shape_out, stride_in, stride_out, axis, forward, data_in, data_out, fct, nthreads) -> Error`

Single-axis variant with explicit strides — see the r2c equivalent.

## Raw foreign procs

| Name                          | Notes |
| ----------------------------- | ----- |
| `c2r_f32`/`c2r_f64`           | Single-axis raw C ABI. |
| `c2r_axes_f32`/`c2r_axes_f64` | Multi-axis: c2c is performed along all axes except `axes.back()`, then c2r along that last one. |

---

# 4. Real-to-real transforms

All r2r wrappers take the same parameter shape:

| Parameter         | Type            | Description |
| ----------------- | --------------- | ----------- |
| `shape`           | `[]uint`        | Array dimensions. |
| `axes`            | `[]uint`        | Axes to transform. |
| `stride_in`       | `[]int`         | Input byte strides per axis. |
| `stride_out`      | `[]int`         | Output byte strides per axis. |
| `data_in`         | `[]f32`/`[]f64` | Input buffer. |
| `data_out`        | matching real type | Output buffer. |
| `fct`             | `f32`/`f64`     | Output scaling. Default `1`. |
| `nthreads`        | `uint`          | Threads. Default `1`. |

The proc-group entry points are: `r2r_fftpack`, `r2r_separable_hartley`,
`r2r_genuine_hartley`, `r2r_separable_fht`, `r2r_genuine_fht`. Each
dispatches to a `_f32_proc` / `_f64_proc` overload based on argument type.

## `r2r_fftpack(shape, axes, stride_in, stride_out, real2hermitian, forward, data_in, data_out, fct, nthreads) -> Error`

FFTPACK-style half-complex transform. Has two extra `bool` parameters:

| Parameter         | Type   | Description |
| ----------------- | ------ | ----------- |
| `real2hermitian`  | `bool` | `true` ⇒ encode the Hermitian spectrum into the real output (forward direction). `false` ⇒ decode it back (backward direction). |
| `forward`         | `bool` | Sign of the complex exponent (independent of `real2hermitian`). |

To round-trip a real signal `x`:

```odin
r2r_fftpack(shape, axes, str, str, true,  true,  x,   freq, 1.0)
r2r_fftpack(shape, axes, str, str, false, false, freq, x2,  1.0/f64(N))
// x2 ≈ x
```

## `r2r_separable_hartley(shape, axes, stride_in, stride_out, data_in, data_out, fct, nthreads) -> Error`
## `r2r_genuine_hartley(shape, axes, stride_in, stride_out, data_in, data_out, fct, nthreads) -> Error`

Hartley transforms.
* **Separable**: applies a 1-D Hartley along each requested axis in turn,
  with `re(x) + im(x)` collapsing complex intermediates.
* **Genuine**: performs a full multidimensional FFT first, then computes
  `re + im` from the result. For a single axis the two are equivalent;
  for higher dimensions they differ.

Both are involutive up to scaling: applying twice with `fct=1` returns
`N * input`.

## `r2r_separable_fht(...) -> Error`
## `r2r_genuine_fht(...) -> Error`

Same parameters as the Hartley variants above. These follow FFTW's "FHT"
sign convention: the imaginary part is **subtracted** instead of added
(`re - im`). Use these unless you specifically need the
backwards-compatible Hartley behavior.

---

# 5. Discrete cosine / sine transforms

## `dct(shape, axes, stride_in, stride_out, type, data_in, data_out, fct, ortho, nthreads) -> Error`
## `dst(shape, axes, stride_in, stride_out, type, data_in, data_out, fct, ortho, nthreads) -> Error`

| Parameter    | Type            | Default | Description |
| ------------ | --------------- | ------- | ----------- |
| `shape`      | `[]uint`        | —       | Array dimensions. |
| `axes`       | `[]uint`        | —       | Axes to transform. |
| `stride_in`  | `[]int`         | —       | Input byte strides per axis. |
| `stride_out` | `[]int`         | —       | Output byte strides per axis. |
| `type`       | `int`           | —       | Transform type, **must be `1`, `2`, `3`, or `4`**. |
| `data_in`    | `[]f32`/`[]f64` | —       | Input buffer. |
| `data_out`   | matching real type | —    | Output buffer. |
| `fct`        | `f32`/`f64`     | `1`     | Output scaling. |
| `ortho`      | `bool`          | `false` | When `true`, apply the boundary corrections that make the transform orthogonal (see below). |
| `nthreads`   | `uint`          | `1`     | Threads. |

### Inverse pairs (with `ortho = false`)

| Forward type | Inverse type | Scale on inverse |
| ------------ | ------------ | ---------------- |
| DCT-II       | DCT-III      | `1 / (2N)`       |
| DCT-I        | DCT-I        | `1 / (2(N-1))`   |
| DCT-IV       | DCT-IV       | `1 / (2N)`       |
| DST-II       | DST-III      | `1 / (2N)`       |
| DST-I        | DST-I        | `1 / (2(N+1))`   |
| DST-IV       | DST-IV       | `1 / (2N)`       |

### Orthogonal flag (`ortho = true`)

For each 1-D sub-transform pocketfft additionally:

* **DCT-I**: multiplies the first/last input by `√2`, divides the
  first/last output by `√2`.
* **DCT-II**: divides the first output by `√2`.
* **DCT-III**: multiplies the first input by `√2`.
* **DCT-IV**: nothing.
* **DST-I**: nothing.
* **DST-II**: divides the last output by `√2`.
* **DST-III**: multiplies the last input by `√2`.
* **DST-IV**: nothing.

These corrections do **not** include the `√(2/N)` global scaling that an
"orthonormal" DCT/DST normally carries; pass it in `fct` if you want the
fully orthonormal definition.

---

# 6. Selecting the right call

| Need                                                | Recommended call               |
| --------------------------------------------------- | ------------------------------ |
| 1-D or row-major contiguous c2c                     | `c2c_simple`                   |
| 1-D or row-major contiguous r2c / c2r              | `r2c_simple` / `c2r_simple`    |
| Multi-axis r2c / c2r                                | `r2c_axes_*` / `c2r_axes_*` (raw) |
| Custom strides (interleaved / column-major / non-contiguous) | `c2c_f64_proc`, `r2c_single_f64`, `c2r_single_f64`, etc. |
| Maximum performance, you handle validation yourself | Raw foreign procs (`c2c_f64`, …) |
| Real-to-real (cosine, sine, Hartley)                | `dct`, `dst`, `r2r_*`          |

---

# 7. Threading

Every wrapper accepts an `nthreads: uint` argument. Pass `0` to use all
logical CPUs. The library silently ignores the value if pocketfft was built
with `POCKETFFT_NO_MULTITHREADING`, or for 1-D transforms.

# 8. Error handling

The wrappers translate three classes of failure to `Error` values:

* `.Empty_Shape` — caught locally before calling into C.
* `.Invalid_Buffer` — caught locally for `c2c_simple` family wrappers when
  the slice you supplied is shorter than the product of `shape`.
* `.Cxx_Exception` / `.Unknown` — pocketfft itself rejected the arguments
  (e.g. duplicate axes, mismatched shape size, invalid DCT/DST type).

The raw foreign procedures return `c.int` directly: `0 = success`,
`1 = std::exception`, `2 = unknown exception`.
