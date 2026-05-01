# PocketFFT - Odin bindings
Bindings for the pocketfft fast FFT library. 

## Description
Odin bindings for the [pocketfft](https://gitlab.mpcdf.mpg.de/mtr/pocketfft) [pocketfft_github](https://github.com/mreineck/pocketfft)
header-only C++ FFT library. <br>
This is a fast library. <br>
Because pocketfft is template-only C++, a thin C ABI wrapper
(`pocketfft_c.{h,cc}`) instantiates the templates for `float` and `double`
(`long double` is omitted as Odin has no equivalent). The Odin package then
links against the resulting static library.

## License
-------

3-clause BSD

## How to compile and run

```
make clean
make pocketfft_lib
make test
make run
```


## Files

| File                     | Purpose                                            |
| ------------------------ | -------------------------------------------------- |
| `pocketfft_hdronly.h`    | Upstream pocketfft header-only library             |
| `pocketfft_c.h`          | C ABI declarations                                 |
| `pocketfft_c.cc`         | C wrapper implementation (instantiates templates)  |
| `Makefile`               | Builds `libpocketfft_c.a`                          |
| `pocketfft.odin`         | Odin bindings + ergonomic wrappers                 |
| `pocketfft_test.odin`    | Test program covering every binding                |

## Build & test

```sh
make                                                          # builds libpocketfft_c.a
odin build . -extra-linker-flags:"-lstdc++ -lpthread -lm" \
             -out:pocketfft_test
./pocketfft_test
```

## API at a glance

Two layers are exposed:

1. **Raw foreign procs** mirroring the C ABI 1-to-1
   (`c2c_f32`, `c2c_f64`, `r2c_f32`, …, `dst_f64`).
2. **Ergonomic Odin wrappers** that take `[]T` slices, validate buffer sizes,
   and return a typed `Error`:

   * `c2c_simple(shape, FORWARD/BACKWARD, src, dst, fct)` — contiguous c2c
   * `r2c_simple` / `c2r_simple`                          — contiguous r2c/c2r
   * `r2r_fftpack`                                        — FFTPACK half-complex
   * `r2r_separable_hartley` / `r2r_genuine_hartley`
   * `r2r_separable_fht`     / `r2r_genuine_fht`
   * `dct(shape, axes, …, type, src, dst, fct, ortho)`
   * `dst(…)`

Helpers:

* `strides_row_major(shape, elem_size)` computes byte strides for a contiguous
  row-major array (pocketfft strides are in **bytes**).

### Minimal example

```odin
import "./pocketfft"

main :: proc( ) {
	
    shape := [ ]uint{ 1024 }
    src   := make( [ ]complex128, 1024 );
    defer delete( src )
    dst   := make( [ ]complex128, 1024 );
    defer delete( dst )
    // ... fill src ...
    err := pocketfft.c2c_simple( shape, pocketfft.FORWARD, src, dst, 1.0 )
    assert( err == .None )
}
```

## Coverage

The included test exercises every wrapped function (both `f32` and `f64`
variants where applicable) and verifies:

* round-trip identity for `c2c`, `r2c`/`c2r`, `r2r_fftpack`,
  `r2r_separable_hartley`/`r2r_genuine_hartley`,
  `r2r_separable_fht`/`r2r_genuine_fht`,
  `dct` (types 1/2-3/4), `dst` (types 1/2-3/4)
* that a known pure-tone signal lands on the correct frequency bin
* error reporting for empty shapes and undersized buffers

All 25 test cases currently pass.

All 2 high performance measure cases pass.


## Have fun
Best regards, <br>
Joao Carvalho <br>
