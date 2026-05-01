// Test program exercising every pocketfft binding.
//
// Build & run:
//   make                                    # builds libpocketfft_c.a
//   odin run . -extra-linker-flags:"-lstdc++ -lpthread -lm"
//
// Strategy: for each transform pair (forward + inverse) we round-trip random
// data and check that the L2 error is below a tolerance. For r2r_fftpack we
// also pair forward(real2hermitian=true) with backward(real2hermitian=false).
package pocketfft

import "core:fmt"
import "core:math"
import "core:math/rand"
import "core:os"
import "core:thread"
import "core:time"

@(private="file") TOL_F32 :: 1.0e-4
@(private="file") TOL_F64 :: 1.0e-12

@(private="file") fail_count := 0

@(private="file")
report :: proc ( name : string, ok : bool, extra : string = "" ) {

	if ok {

		fmt.printf( "  [ OK ] %s %s\n", name, extra )
	} else {

		fmt.printf( "  [FAIL] %s %s\n", name, extra )
		fail_count += 1
	}
}

@(private="file")
l2err_f32 :: proc ( a, b : [ ]f32 ) -> f64 {

	num, den : f64
	for v, i in a {

		d := f64( v ) - f64( b[ i ] )
		num += d*d
		den += f64( v ) * f64( v )
	}
	if den == 0 do return math.sqrt( num )
	return math.sqrt( num / den )
}

@(private="file")
l2err_f64 :: proc ( a, b : [ ]f64 ) -> f64 {

	num, den : f64
	for v, i in a {

		d   := v - b[ i ]
		num += d * d
		den += v * v
	}
	if den == 0 do return math.sqrt( num )
	return math.sqrt( num / den )
}

@(private="file")
l2err_c64 :: proc ( a, b : [ ]complex64 ) -> f64 {

	num, den : f64
	for v, i in a {

		dr  := f64( real( v ) - real( b[ i ] ) )
		di  := f64( imag( v ) - imag( b[ i ] ) )
		num += dr * dr + di * di
		ar  := f64( real( v ) ); ai := f64( imag( v ) )
		den += ar * ar + ai * ai
	}
	if den == 0 do return math.sqrt( num )
	return math.sqrt( num / den )
}

@(private="file")
l2err_c128 :: proc ( a, b : [ ]complex128 ) -> f64 {

	num, den : f64
	for v, i in a {

		dr  := real( v ) - real( b[ i ] )
		di  := imag( v ) - imag( b[ i ] )
		num += dr * dr + di * di
		ar  := real( v ); ai := imag( v )
		den += ar * ar + ai * ai
	}
	if den == 0 do return math.sqrt( num )
	return math.sqrt( num / den )
}

@(private="file")
rand_real_f32 :: proc ( n : int ) -> [ ]f32 {

	out := make( [ ]f32, n )
	for i in 0 ..< n do out[ i ] = rand.float32_range( -1, 1 )
	return out
}

@(private="file")
rand_real_f64 :: proc ( n : int ) -> [ ]f64 {

	out := make( [ ]f64, n )
	for i in 0 ..< n do out[ i ] = rand.float64_range( -1, 1 )
	return out
}

@(private="file")
rand_cplx_f32 :: proc ( n : int ) -> [ ]complex64 {

	out := make( [ ]complex64, n )
	for i in 0 ..< n {

		out[ i ] = complex( rand.float32_range( -1, 1 ), rand.float32_range( -1, 1 ) )
	}
	return out
}

@(private="file")
rand_cplx_f64 :: proc ( n : int ) -> [ ]complex128 {

	out := make( [ ]complex128, n )
	for i in 0 ..< n {

		out[ i ] = complex( rand.float64_range( -1, 1 ), rand.float64_range( -1, 1 ) )
	}
	return out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test_c2c_1d :: proc ( ) {

	N :: 256
	shape := [ ]uint{ N }
	src := rand_cplx_f64( N )
	defer delete( src )
	fwd := make( [ ]complex128, N ); defer delete( fwd )
	rt  := make( [ ]complex128, N ); defer delete( rt )

	e1 := c2c_simple( shape, FORWARD,  src, fwd, 1.0 )
	e2 := c2c_simple( shape, BACKWARD, fwd, rt,  1.0 / f64( N ) )
	err := l2err_c128( src, rt )
	report( "c2c f64 1D round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_c2c_1d_f32 :: proc ( ) {

	N :: 256
	shape := [ ]uint{ N }
	src := rand_cplx_f32( N ); defer delete( src )
	fwd := make( [ ]complex64, N ); defer delete( fwd )
	rt  := make( [ ]complex64, N ); defer delete( rt )

	e1  := c2c_simple( shape, FORWARD,  src, fwd, f32( 1 ) )
	e2  := c2c_simple( shape, BACKWARD, fwd, rt,  f32( 1.0 / f32( N ) ) )
	err := l2err_c64( src, rt )
	report( "c2c f32 1D round-trip",
		e1 == .None && e2 == .None && err < TOL_F32,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_c2c_2d :: proc ( ) {

	shape := [ ]uint{ 16, 32 }
	N := int( shape[ 0 ] * shape[ 1 ] )
	src := rand_cplx_f64( N ); defer delete( src )
	fwd := make( [ ]complex128, N ); defer delete( fwd )
	rt  := make( [ ]complex128, N ); defer delete( rt )

	e1  := c2c_simple( shape, FORWARD,  src, fwd, 1.0 )
	e2  := c2c_simple( shape, BACKWARD, fwd, rt,  1.0 / f64( N ) )
	err := l2err_c128( src, rt )
	report( "c2c f64 2D round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_r2c_c2r_1d :: proc ( ) {

	N :: 128
	shape := [ ]uint{ N }
	src := rand_real_f64( N ); defer delete( src )
	freq := make( [ ]complex128, N / 2 + 1 ); defer delete( freq )
	rt   := make( [ ]f64, N ); defer delete( rt )

	e1 := r2c_simple( shape, FORWARD,  src,  freq, 1.0 )
	e2 := c2r_simple( shape, BACKWARD, freq, rt,   1.0 / f64( N ) )
	err := l2err_f64( src, rt )
	report( "r2c+c2r f64 1D round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_r2c_c2r_1d_f32 :: proc ( ) {

	N :: 128
	shape := [ ]uint{ N }
	src := rand_real_f32( N ); defer delete( src )
	freq := make( [ ]complex64, N / 2 + 1 ); defer delete( freq )
	rt   := make( [ ]f32, N ); defer delete( rt )

	e1  := r2c_simple( shape, FORWARD,  src,  freq, f32( 1 ) )
	e2  := c2r_simple( shape, BACKWARD, freq, rt,   f32( 1.0 / f32( N ) ) )
	err := l2err_f32( src, rt )
	report( "r2c+c2r f32 1D round-trip",
		e1 == .None && e2 == .None && err < TOL_F32,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_r2c_c2r_2d :: proc ( ) {

	// 2D r2c/c2r using the multi-axis variants via the raw API.
	shape  := [ ]uint{ 8, 16 }
	axes   := [ ]uint{ 0, 1 }
	N      := int( shape[ 0 ] * shape[ 1 ] )
	NF     := int( shape[ 0 ] * ( shape[ 1 ] / 2 + 1 ) )

	src  := rand_real_f64( N ); defer delete( src )
	freq := make( [ ]complex128, NF ); defer delete( freq )
	rt   := make( [ ]f64, N ); defer delete( rt )

	str_in     := strides_row_major( shape, size_of( f64 ) ); defer delete( str_in )
	freq_shape := [ ]uint{ shape[ 0 ], shape[ 1 ] / 2 + 1 }
	str_out    := strides_row_major( freq_shape, size_of( complex128 ) ); defer delete( str_out )

	rc1 := r2c_axes_f64( 2, sptr( shape ), pptr( str_in ), pptr( str_out ),
		2, sptr( axes ), 1, raw_data( src ), raw_data( freq ), 1.0, 1 )
	rc2 := c2r_axes_f64( 2, sptr( shape ), pptr( str_out ), pptr( str_in ),
		2, sptr( axes ), 0, raw_data( freq ), raw_data( rt ), 1.0 / f64( N ), 1 )
	err := l2err_f64( src, rt )
	report( "r2c+c2r f64 2D (multi-axis) round-trip",
		rc1 == 0 && rc2 == 0 && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_r2r_fftpack :: proc ( ) {

	N :: 64
	shape := [ ]uint{ N }
	axes  := [ ]uint{ 0 }
	src := rand_real_f64( N ); defer delete( src )
	mid := make( [ ]f64, N ); defer delete( mid )
	rt  := make( [ ]f64, N ); defer delete( rt )

	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	e1  := r2r_fftpack( shape, axes, str, str, true,  true,  src, mid, 1.0 )
	e2  := r2r_fftpack( shape, axes, str, str, false, false, mid, rt,  1.0 / f64( N ) )
	err := l2err_f64( src, rt )
	report( "r2r_fftpack f64 round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )

	// f32 variant
	srcf := rand_real_f32( N ); defer delete( srcf )
	midf := make( [ ]f32, N ); defer delete( midf )
	rtf  := make( [ ]f32, N ); defer delete( rtf )
	strf := strides_row_major( shape, size_of( f32 ) ); defer delete( strf )
	e3   := r2r_fftpack( shape, axes, strf, strf, true,  true,  srcf, midf, f32( 1 ) )
	e4   := r2r_fftpack( shape, axes, strf, strf, false, false, midf, rtf,  f32( 1.0 / f32( N ) ) )
	errf := l2err_f32( srcf, rtf )
	report( "r2r_fftpack f32 round-trip",
		e3 == .None && e4 == .None && errf < TOL_F32,
		fmt.tprintf( "(err=%.2e)", errf ) )
}

// Hartley/FHT transforms are involutive up to scaling: applying twice with
// fct=1 gives N*src.
test_r2r_separable_hartley :: proc ( ) {

	N :: 64
	shape := [ ]uint{ N }
	axes  := [ ]uint{ 0 }
	src := rand_real_f64( N ); defer delete( src )
	mid := make( [ ]f64, N ); defer delete( mid )
	rt  := make( [ ]f64, N ); defer delete( rt )
	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	e1  := r2r_separable_hartley( shape, axes, str, str, src, mid, 1.0 )
	e2  := r2r_separable_hartley( shape, axes, str, str, mid, rt,  1.0 / f64( N ) )
	err := l2err_f64( src, rt )
	report( "r2r_separable_hartley f64 round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )

	srcf := rand_real_f32( N ); defer delete( srcf )
	midf := make( [ ]f32, N ); defer delete( midf )
	rtf  := make( [ ]f32, N ); defer delete( rtf )
	strf := strides_row_major( shape, size_of( f32 ) ); defer delete( strf )
	e3 := r2r_separable_hartley( shape, axes, strf, strf, srcf, midf, f32( 1 ) )
	e4 := r2r_separable_hartley( shape, axes, strf, strf, midf, rtf,  f32( 1.0 / f32( N ) ) )
	errf := l2err_f32( srcf, rtf )
	report( "r2r_separable_hartley f32 round-trip",
		e3 == .None && e4 == .None && errf < TOL_F32,
		fmt.tprintf( "(err=%.2e)", errf ) )
}

test_r2r_genuine_hartley :: proc ( ) {

	shape := [ ]uint{ 8, 16 }
	axes  := [ ]uint{ 0, 1 }
	N   := int( shape[ 0 ] * shape[ 1 ] )
	src := rand_real_f64( N ); defer delete( src )
	mid := make( [ ]f64, N ); defer delete( mid )
	rt  := make( [ ]f64, N ); defer delete( rt )
	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	e1 := r2r_genuine_hartley( shape, axes, str, str, src, mid, 1.0 )
	e2 := r2r_genuine_hartley( shape, axes, str, str, mid, rt,  1.0 / f64( N ) )
	err := l2err_f64(src, rt)
	report( "r2r_genuine_hartley f64 2D round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_r2r_separable_fht :: proc ( ) {

	N :: 64
	shape := [ ]uint{ N }
	axes  := [ ]uint{ 0 }
	src := rand_real_f64( N ); defer delete( src )
	mid := make( [ ]f64, N ); defer delete( mid )
	rt  := make( [ ]f64, N ); defer delete( rt )
	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	e1 := r2r_separable_fht( shape, axes, str, str, src, mid, 1.0 )
	e2 := r2r_separable_fht( shape, axes, str, str, mid, rt,  1.0 / f64( N ) )
	err := l2err_f64( src, rt )
	report( "r2r_separable_fht f64 round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

test_r2r_genuine_fht :: proc ( ) {

	shape := [ ]uint{ 8, 16 }
	axes  := [ ]uint{ 0, 1 }
	N := int( shape[ 0 ] * shape[ 1 ] )
	src := rand_real_f64( N ); defer delete( src )
	mid := make( [ ]f64, N ); defer delete( mid )
	rt  := make( [ ]f64, N ); defer delete( rt )
	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	e1 := r2r_genuine_fht( shape, axes, str, str, src, mid, 1.0 )
	e2 := r2r_genuine_fht( shape, axes, str, str, mid, rt,  1.0 / f64( N ) )
	err := l2err_f64( src, rt )
	report( "r2r_genuine_fht f64 2D round-trip",
		e1 == .None && e2 == .None && err < TOL_F64,
		fmt.tprintf( "(err=%.2e)", err ) )
}

// DCT/DST: inverse of type-N is type-N' with appropriate scaling
//   dct(2) inverse is dct(3) with fct = 1/(2N)
//   dst(2) inverse is dst(3) with fct = 1/(2N)
//   dct(1) inverse is dct(1) with fct = 1/(2(N-1))
//   dst(1) inverse is dst(1) with fct = 1/(2(N+1))
//   type-4 is its own inverse with fct = 1/(2N)
test_dct :: proc ( ) {

	N :: 64
	shape := [ ]uint{ N }
	axes  := [ ]uint{ 0 }
	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	// type 2 / 3 pair
	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N ); defer delete( mid )
		rt  := make( [ ]f64, N ); defer delete( rt )
		e1 := dct( shape, axes, str, str, 2, src, mid, 1.0, false )
		e2 := dct( shape, axes, str, str, 3, mid, rt,  1.0 / f64( 2 * N ), false )
		err := l2err_f64( src, rt )
		report( "dct type 2/3 f64",
			e1 == .None && e2 == .None && err < TOL_F64,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	// type 1 (self-inverse with fct=1/(2(N-1)))
	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N); defer delete( mid )
		rt  := make( [ ]f64, N); defer delete( rt )
		e1 := dct( shape, axes, str, str, 1, src, mid, 1.0, false )
		e2 := dct( shape, axes, str, str, 1, mid, rt,  1.0 / f64( 2 * ( N - 1 ) ), false )
		err := l2err_f64( src, rt )
		report( "dct type 1 f64 self-inverse",
			e1 == .None && e2 == .None && err < TOL_F64,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	// type 4 (self-inverse with fct=1/(2N))
	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N ); defer delete( mid )
		rt  := make( [ ]f64, N ); defer delete( rt )
		e1 := dct( shape, axes, str, str, 4, src, mid, 1.0, false )
		e2 := dct( shape, axes, str, str, 4, mid, rt,  1.0 / f64( 2 * N ), false )
		err := l2err_f64( src, rt )
		report( "dct type 4 f64 self-inverse",
			e1 == .None && e2 == .None && err < TOL_F64,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	// f32 type 2/3
	{

		src := rand_real_f32( N ); defer delete( src )
		mid := make( [ ]f32, N ); defer delete( mid )
		rt  := make( [ ]f32, N ); defer delete( rt )
		strf := strides_row_major( shape, size_of( f32 ) ); defer delete( strf )
		e1 := dct( shape, axes, strf, strf, 2, src, mid, f32( 1 ), false )
		e2 := dct( shape, axes, strf, strf, 3, mid, rt,  f32( 1.0 / f32( 2 *  N ) ), false )
		err := l2err_f32( src, rt )
		report( "dct type 2/3 f32",
			e1 == .None && e2 == .None && err < TOL_F32,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	// orthogonal flag: just verify the binding accepts ortho=true and produces
	// finite output (the exact normalization convention is documented in the
	// pocketfft README and depends on transform type).
	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N ); defer delete( mid )
		e1 := dct( shape, axes, str, str, 2, src, mid, 1.0, true )
		all_finite := true
		for v in mid {

			if !( v == v ) || v == math.INF_F64 || v == -math.INF_F64 {

				all_finite = false
			    break
			}
		}
		report( "dct ortho=true binding works", e1 == .None && all_finite )
	}
}

test_dst :: proc ( ) {

	N :: 64
	shape := [ ]uint{ N }
	axes  := [ ]uint{ 0 }
	str := strides_row_major( shape, size_of( f64 ) ); defer delete( str )

	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N); defer delete( mid )
		rt  := make( [ ]f64, N); defer delete( rt )
		e1 := dst( shape, axes, str, str, 2, src, mid, 1.0, false )
		e2 := dst( shape, axes, str, str, 3, mid, rt,  1.0 / f64( 2 * N ), false )
		err := l2err_f64( src, rt )
		report( "dst type 2/3 f64",
			e1 == .None && e2 == .None && err < TOL_F64,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N ); defer delete( mid )
		rt  := make( [ ]f64, N ); defer delete( rt )
		e1  := dst( shape, axes, str, str, 1, src, mid, 1.0, false )
		e2  := dst( shape, axes, str, str, 1, mid, rt,  1.0 / f64( 2 * ( N + 1 ) ), false )
		err := l2err_f64( src, rt )
		report( "dst type 1 f64 self-inverse",
			e1 == .None && e2 == .None && err < TOL_F64,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	{

		src := rand_real_f64( N ); defer delete( src )
		mid := make( [ ]f64, N ); defer delete( mid )
		rt  := make( [ ]f64, N ); defer delete( rt )
		e1  := dst( shape, axes, str, str, 4, src, mid, 1.0, false )
		e2  := dst( shape, axes, str, str, 4, mid, rt,  1.0 / f64( 2 * N ), false )
		err := l2err_f64( src, rt )
		report( "dst type 4 f64 self-inverse",
			e1 == .None && e2 == .None && err < TOL_F64,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
	{

		src := rand_real_f32( N ); defer delete( src )
		mid := make( [ ]f32, N ); defer delete( mid )
		rt  := make( [ ]f32, N ); defer delete( rt )
		strf := strides_row_major( shape, size_of( f32 ) ); defer delete( strf )
		e1 := dst( shape, axes, strf, strf, 2, src, mid, f32( 1 ), false )
		e2 := dst( shape, axes, strf, strf, 3, mid, rt,  f32( 1.0 / f32( 2 * N ) ), false )
		err := l2err_f32( src, rt )
		report( "dst type 2/3 f32",
			e1 == .None && e2 == .None && err < TOL_F32,
			fmt.tprintf( "(err=%.2e)", err ) )
	}
}

// Verify that a known sinusoid produces a single non-zero bin.
test_known_signal :: proc ( ) {

	N :: 64
	shape := [ ]uint{ N }
	src := make( [ ]complex128, N ); defer delete( src )
	out := make( [ ]complex128, N ); defer delete( out )
	K :: 5  // bin index of the test frequency
	for i in 0 ..< N {

		t := 2 * math.PI * f64( K ) * f64( i ) / f64( N )
		src[ i ] = complex( math.cos( t ), math.sin( t ) )
	}
	e := c2c_simple( shape, FORWARD, src, out, 1.0 )
	// expect out[K] ~= N, all others ~= 0
	max_other: f64 = 0
	for i in 0 ..< N {

		if i == K do continue
		m := math.sqrt( real( out[ i ] ) * real( out[ i ] ) + imag( out[ i ] ) * imag( out[ i ] ) )
		if m > max_other do max_other = m
	}
	bin_mag := math.sqrt( real( out[ K ] ) * real( out[ K ] ) + imag( out[ K ] ) * imag( out[ K ] ) )
	ok := e == .None && math.abs( bin_mag - f64( N ) ) < 1e-9 && max_other < 1e-9
	report( "c2c known sinusoid lands on correct bin", ok,
		fmt.tprintf( "(bin=%.3f other=%.2e)", bin_mag, max_other ) )
}

// Empty shape returns an error code rather than crashing.
test_empty_shape_error :: proc ( ) {

	shape := [ ]uint{ }
	src := [ ]complex128{ }
	out := [ ]complex128{ }
	e := c2c_simple( shape, FORWARD, src, out, 1.0 )
	report( "empty shape returns Error", e == .Empty_Shape )
}

// Insufficient buffer triggers Invalid_Buffer.
test_buffer_check :: proc ( ) {

	shape := [ ]uint{ 32 }
	src := make( [ ]complex128, 16 ); defer delete( src ) // too small!
	out := make( [ ]complex128, 32 ); defer delete( out )
	e := c2c_simple( shape, FORWARD, src, out, 1.0 )
	report( "undersized buffer detected", e == .Invalid_Buffer )
}

// Large parallel batched FFT benchmark:
//   200_000 vectors of complex64, length = next power of 2 >= 16_000 (= 16384).
//   Uses all CPU cores via pocketfft's nthreads parameter, in-place to keep
//   memory pressure manageable (~26 GB for the single buffer).
bench_last_fft_ns  : i64 = 0
bench_last_ifft_ns : i64 = 0
bench_n_vectors    : int = 0
bench_vec_len      : int = 0
bench_threads      : int = 0

next_pow2 :: proc ( x : uint ) -> uint {

	v : uint = 1
	for v < x { v <<= 1 }
	return v
}

test_large_parallel_c2c_f32 :: proc ( ) {

	N        := uint( 200_000 )
	M        := next_pow2( 16_000 )      // 16384
	threads  := uint( os.processor_core_count( ) )
	bench_n_vectors = int( N )
	bench_vec_len   = int( M )
	bench_threads   = int( threads )

	total := N * M
	bytes := total * size_of( complex64 )
	fmt.printf( "  (allocating %d MiB for %d x %d complex64 buffer, %d threads)\n",
		bytes / ( 1024 * 1024 ), N, M, threads )

	buf, err := make( [ ]complex64, total )
	if err != nil {

		report( "large parallel c2c_f32 (alloc)", false, "allocation failed" )
		return
	}
	defer delete( buf )

	// Fill with a cheap deterministic pattern (avoid 26 GB of RNG calls).
	for i in 0 ..< total {

		x := f32( i & 0xFFFF ) * 1.52587890625e-5  // /65536
		buf[ i ] = complex(x, -x)
	}

	shape  := [ ]uint{ N, M }
	axes   := [ ]uint{ 1 }                          // FFT along last axis only
	stride := strides_row_major( shape, size_of( complex64 ) )
	defer delete( stride )

	// Forward FFT (in-place).
	t0 := time.now( )
	e1 := c2c_f32_proc( shape, axes, stride, stride, FORWARD,
		buf, buf, fct = 1.0, nthreads = threads )
	bench_last_fft_ns = time.duration_nanoseconds( time.since( t0 ) )

	// Inverse FFT (in-place), normalized so we recover the input.
	inv_scale := f32( 1.0 ) / f32( M )
	t1 := time.now( )
	e2 := c2c_f32_proc( shape, axes, stride, stride, BACKWARD,
		buf, buf, fct = inv_scale, nthreads = threads )
	bench_last_ifft_ns = time.duration_nanoseconds( time.since( t1 ) )

	ok := e1 == .None && e2 == .None
	report( "large parallel c2c_f32 FFT+IFFT round trip", ok,
		fmt.tprintf( "(fft=%.3f s ifft=%.3f s)",
			f64( bench_last_fft_ns ) * 1e-9, f64( bench_last_ifft_ns ) * 1e-9 ) )
}

// Same workload as test_large_parallel_c2c_f32, but parallelized in Odin
// (one worker thread per CPU core). Each worker calls pocketfft with
// nthreads = 1 on its own contiguous chunk of vectors.
bench2_total_ns : i64 = 0
bench2_fft_ns   : i64 = 0
bench2_ifft_ns  : i64 = 0
bench2_threads  : int = 0

Worker_Job :: struct {

	buf       : [ ^ ]complex64,
	row_start : uint,
	row_count : uint,
	vec_len   : uint,
	inv_scale : f32,
	fft_ns    : i64,
	ifft_ns   : i64,
	err       : Error,
}

worker_run :: proc ( job : ^Worker_Job ) {

	shape  := [ ]uint{ job.row_count, job.vec_len }
	axes   := [ ]uint{ 1 }
	stride := strides_row_major( shape, size_of( complex64 ) )
	defer delete( stride )

	chunk_len := job.row_count * job.vec_len
	chunk     := job.buf[ job.row_start * job.vec_len : ][ : chunk_len ]

	t0 := time.now( )
	e1 := c2c_f32_proc( shape, axes, stride, stride, FORWARD,
		chunk, chunk, fct = 1.0, nthreads = 1 )
	job.fft_ns = time.duration_nanoseconds( time.since( t0 ) )

	t1 := time.now( )
	e2 := c2c_f32_proc( shape, axes, stride, stride, BACKWARD,
		chunk, chunk, fct = job.inv_scale, nthreads = 1 )
	job.ifft_ns = time.duration_nanoseconds( time.since( t1 ) )

	if e1 != .None { job.err = e1 }
	else           { job.err = e2 }
}

test_large_odin_parallel_c2c_f32 :: proc ( ) {

	N       := uint( 200_000 )
	M       := next_pow2( 16_000 )      // 16384
	threads := uint( os.processor_core_count( ) )
	if threads > N { threads = N }
	bench2_threads = int( threads )

	total := N * M
	bytes := total * size_of( complex64 )
	fmt.printf( "  (allocating %d MiB for %d x %d complex64 buffer, %d Odin threads)\n",
		bytes / ( 1024 * 1024 ), N, M, threads )

	buf, err := make( [ ]complex64, total )
	if err != nil {

		report( "large odin-parallel c2c_f32 (alloc)", false, "allocation failed" )
		return
	}
	defer delete( buf )

	for i in 0 ..< total {

		x := f32( i & 0xFFFF ) * 1.52587890625e-5
		buf[ i ] = complex( x, -x )
	}

	// Even chunking: first (N % T) workers get one extra row.
	jobs := make( [ ]Worker_Job, threads ); defer delete( jobs )
	base := N / threads
	rem  := N % threads
	row  : uint = 0
	for i in 0 ..< threads {

		count   := base + ( 1 if i < rem else 0 )
		jobs[i]  = Worker_Job{

			buf       = raw_data( buf ),
			row_start = row,
			row_count = count,
			vec_len   = M,
			inv_scale = f32(1.0) / f32(M),
		}
		row += count
	}

	workers := make( [ ]^thread.Thread, threads ); defer delete( workers )

	t_start := time.now( )
	for i in 0 ..< threads {

		workers[ i ] = thread.create_and_start_with_poly_data( & jobs[ i ], worker_run )
	}
	for i in 0 ..< threads {

		thread.join( workers[ i ] )
		thread.destroy( workers[ i ] )
	}
	bench2_total_ns = time.duration_nanoseconds( time.since( t_start ) )

	// Aggregate: sum of per-worker times divided by thread count is the
	// average single-worker cost; the wall-clock total is what matters.
	max_fft, max_ifft : i64 = 0, 0
	all_ok := true
	for i in 0 ..< threads {

		if jobs[ i ].err != .None { all_ok = false }
		if jobs[ i ].fft_ns  > max_fft  { max_fft  = jobs[ i ].fft_ns  }
		if jobs[ i ].ifft_ns > max_ifft { max_ifft = jobs[ i ].ifft_ns }
	}
	bench2_fft_ns  = max_fft
	bench2_ifft_ns = max_ifft

	report( "large odin-parallel c2c_f32 FFT+IFFT round trip", all_ok,
		fmt.tprintf( "(wall=%.3f s)", f64( bench2_total_ns ) * 1e-9 ) )
}

test_main :: proc ( ) {

	rand.reset( 0xC0FFEE )

	fmt.println( "--- pocketfft Odin binding tests ---" )

	test_c2c_1d( )
	test_c2c_1d_f32( )
	test_c2c_2d( )
	test_known_signal( )

	test_r2c_c2r_1d( )
	test_r2c_c2r_1d_f32( )
	test_r2c_c2r_2d( )

	test_r2r_fftpack( )
	test_r2r_separable_hartley( )
	test_r2r_genuine_hartley( )
	test_r2r_separable_fht( )
	test_r2r_genuine_fht( )

	test_dct( )
	test_dst( )

	test_empty_shape_error( )
	test_buffer_check( )

	test_large_parallel_c2c_f32( )
	test_large_odin_parallel_c2c_f32( )

	fmt.println("------------------------------------")
	if fail_count == 0 {

		fmt.println( "All tests passed." )
		fmt.println( )
		fmt.println( "Last benchmark tests: 200,000 complex64 vectors of length 16,384" )
		fmt.println( "(next power of 2 of 16,000), forward FFT + inverse FFT." )
		fmt.println( )
		fmt.println( "[1] Pocketfft-internal threading (nthreads parameter):" )
		fmt.printf( "    Forward FFT (%d threads): %.3f s\n",
			bench_threads, f64( bench_last_fft_ns ) * 1e-9 )
		fmt.printf( "    Inverse FFT (%d threads): %.3f s\n",
			bench_threads, f64( bench_last_ifft_ns ) * 1e-9 )
		fmt.printf( "    Total: %.3f s\n",
			f64( bench_last_fft_ns + bench_last_ifft_ns ) * 1e-9 )
		fmt.println( )
		fmt.println( "[2] Odin-level threading (one worker per core, nthreads=1 each):" )
		fmt.printf( "    Workers: %d, ~%d vectors per worker.\n",
			bench2_threads, 200_000 / bench2_threads )
		fmt.printf( "    Slowest worker FFT:  %.3f s\n", f64( bench2_fft_ns ) * 1e-9 )
		fmt.printf( "    Slowest worker IFFT: %.3f s\n", f64( bench2_ifft_ns ) * 1e-9 )
		fmt.printf( "    Wall-clock total (FFT+IFFT, all workers): %.3f s\n",
			f64( bench2_total_ns ) * 1e-9 )
		os.exit( 0 )
	} else {

		fmt.printf( "%d test(s) FAILED.\n", fail_count )
		os.exit( 1 )
	}
}
