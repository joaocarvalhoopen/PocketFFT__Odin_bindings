// Odin bindings for pocketfft (https://gitlab.mpcdf.mpg.de/mtr/pocketfft).
//
// The bindings link against `libpocketfft_c.a`, a thin C wrapper around the
// header-only C++ template library. Build it with `make` first.
//
// Two layers are provided:
//   * raw `foreign` bindings ( lower-case prefix `c_*` ) that mirror the C ABI
//   * ergonomic Odin wrappers that take [ ]T slices, compute strides for
//     contiguous row-major arrays automatically, and return a typed Error.
//
// Only `f32` and `f64` are supported ( Odin has no native long double ).
package pocketfft

import "core:c"

// ---------------------------------------------------------------------------
// Foreign declarations ( raw C ABI ).
// ---------------------------------------------------------------------------

when ODIN_OS == .Windows {

	foreign import pocketfft_c "../bin/libpocketfft_c.lib"
} else {

	foreign import pocketfft_c "../bin/libpocketfft_c.a"
}

@(default_calling_convention = "c", link_prefix = "pocketfft_")
foreign pocketfft_c {

	c2c_f32 :: proc ( ndim : c.size_t, shape: [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t, forward : c.int,
		data_in : rawptr, data_out : rawptr, fct : f32, nthreads : c.size_t ) -> c.int ---

	c2c_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t, forward : c.int,
		data_in : rawptr, data_out : rawptr, fct : f64, nthreads : c.size_t ) -> c.int ---

	r2c_f32 :: proc ( ndim : c.size_t, shape_in : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		axis : c.size_t, forward : c.int,
		data_in : [ ^ ]f32, data_out : rawptr, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2c_f64 :: proc ( ndim : c.size_t, shape_in : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		axis : c.size_t, forward : c.int,
		data_in : [ ^ ]f64, data_out : rawptr, fct : f64, nthreads : c.size_t ) -> c.int ---

	r2c_axes_f32 :: proc ( ndim : c.size_t, shape_in : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out: [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t, forward : c.int,
		data_in : [ ^ ]f32, data_out : rawptr, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2c_axes_f64 :: proc ( ndim : c.size_t, shape_in : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t, forward : c.int,
		data_in : [ ^ ]f64, data_out : rawptr, fct : f64, nthreads : c.size_t ) -> c.int ---

	c2r_f32 :: proc ( ndim : c.size_t, shape_out : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		axis : c.size_t, forward : c.int,
		data_in : rawptr, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	c2r_f64 :: proc ( ndim : c.size_t, shape_out : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		axis : c.size_t, forward : c.int,
		data_in : rawptr, data_out : [ ^ ]f64, fct : f64, nthreads : c.size_t ) -> c.int ---

	c2r_axes_f32 :: proc ( ndim : c.size_t, shape_out : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t, forward : c.int,
		data_in : rawptr, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	c2r_axes_f64 :: proc ( ndim : c.size_t, shape_out : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t, forward : c.int,
		data_in : rawptr, data_out : [ ^ ]f64, fct : f64, nthreads : c.size_t ) -> c.int ---

	r2r_fftpack_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		real2hermitian : c.int, forward : c.int,
		data_in : [ ^ ]f32, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2r_fftpack_f64 :: proc ( ndim : c.size_t, shape: [ ^ ]c.size_t,
		stride_in: [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		real2hermitian : c.int, forward : c.int,
		data_in : [ ^ ]f64, data_out : [ ^ ]f64, fct : f64, nthreads : c.size_t ) -> c.int ---

	r2r_separable_hartley_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f32, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2r_separable_hartley_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f64, data_out: [ ^ ]f64, fct : f64, nthreads : c.size_t ) -> c.int ---

	r2r_genuine_hartley_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f32, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2r_genuine_hartley_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f64, data_out : [ ^ ]f64, fct: f64, nthreads : c.size_t ) -> c.int ---

	r2r_separable_fht_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f32, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2r_separable_fht_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f64, data_out : [ ^ ]f64, fct : f64, nthreads : c.size_t ) -> c.int ---

	r2r_genuine_fht_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes: [ ^ ]c.size_t,
		data_in : [ ^ ]f32, data_out : [ ^ ]f32, fct : f32, nthreads : c.size_t ) -> c.int ---

	r2r_genuine_fht_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		data_in : [ ^ ]f64, data_out : [ ^ ]f64, fct : f64, nthreads : c.size_t ) -> c.int ---

	dct_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		type : c.int, data_in : [ ^ ]f32, data_out : [ ^ ]f32,
		fct : f32, ortho : c.int, nthreads : c.size_t ) -> c.int ---

	dct_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		type : c.int, data_in : [ ^ ]f64, data_out : [ ^ ]f64,
		fct : f64, ortho : c.int, nthreads : c.size_t ) -> c.int ---

	dst_f32 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		type : c.int, data_in : [ ^ ]f32, data_out : [ ^ ]f32,
		fct : f32, ortho : c.int, nthreads : c.size_t ) -> c.int ---

	dst_f64 :: proc ( ndim : c.size_t, shape : [ ^ ]c.size_t,
		stride_in : [ ^ ]c.ptrdiff_t, stride_out : [ ^ ]c.ptrdiff_t,
		naxes : c.size_t, axes : [ ^ ]c.size_t,
		type : c.int, data_in : [ ^ ]f64, data_out : [ ^ ]f64,
		fct : f64, ortho : c.int, nthreads : c.size_t ) -> c.int ---
}

// ---------------------------------------------------------------------------
// Ergonomic Odin API.
// ---------------------------------------------------------------------------

FORWARD  :: true
BACKWARD :: false

Error :: enum {

	None             = 0,
	Cxx_Exception    = 1,  // std::exception thrown inside pocketfft
	Unknown          = 2,  // non-std exception
	Invalid_Buffer   = 3,  // input/output slice too small for given shape
	Empty_Shape      = 4,  // shape has 0 dimensions
}

// Compute byte-strides for a contiguous row-major array of `elem_size`-byte
// elements with the given shape. Caller owns the returned slice.
strides_row_major :: proc( shape : [ ]uint, elem_size: int, allocator := context.allocator ) -> [ ]int {

	out := make( [ ]int, len( shape ), allocator )
	if len( shape ) == 0 do return out
	tmp := elem_size
	for i := len( shape ) - 1; i >= 0; i -= 1 {

		out[ i ] = tmp
		tmp *= int( shape[ i ] )
	}

	return out
}

@(private)
prod :: proc( shape : [ ]uint ) -> uint {

	n : uint = 1
	for s in shape do n *= s
	return n
}

@(private)
to_axes :: proc ( shape : [ ]uint, allocator := context.allocator ) -> [ ]uint {

	out := make( [ ]uint, len( shape ), allocator )
	for i in 0 ..< len( shape ) do out[ i ] = uint( i )
	return out
}

// Cast helper: [ ]uint -> [ ^ ]c.size_t ( works because both are usize on common targets ).
@(private)
sptr :: proc ( s : [ ]uint ) -> [ ^ ]c.size_t {

	return cast( [ ^ ]c.size_t )raw_data( s ) if len( s ) > 0 else nil
}

@(private)
pptr :: proc ( s : [ ]int ) -> [ ^ ]c.ptrdiff_t {

	return cast( [ ^ ]c.ptrdiff_t )raw_data( s ) if len( s ) > 0 else nil
}

// -- complex-to-complex ------------------------------------------------------

c2c :: proc{ c2c_f32, c2c_f64 }

c2c_f32_proc :: proc ( shape, axes : [ ]uint,
	stride_in, stride_out : [ ]int,
	forward : bool,
	data_in, data_out : [ ]complex64,
	fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	if uint( len( data_in ) ) < prod( shape ) || uint( len( data_out ) ) < prod( shape ) do return .Invalid_Buffer
	rc := c2c_f32( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ), 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

c2c_f64_proc :: proc ( shape, axes : [ ]uint,
	stride_in, stride_out : [ ]int,
	forward : bool,
	data_in, data_out : [ ]complex128,
	fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	if uint( len( data_in ) ) < prod( shape ) || uint(len( data_out ) ) < prod( shape ) do return .Invalid_Buffer
	rc := c2c_f64( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ), 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// Convenience: contiguous c2c, all axes, computes strides.
c2c_simple :: proc{ c2c_simple_f32, c2c_simple_f64 }

c2c_simple_f32 :: proc ( shape : [ ]uint, forward : bool,
	data_in, data_out : [ ]complex64, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	stride := strides_row_major( shape, size_of( complex64 ) )
	defer delete( stride )
	axes := to_axes( shape )
	defer delete( axes )
	return c2c_f32_proc( shape, axes, stride, stride, forward, data_in, data_out, fct, nthreads )
}

c2c_simple_f64 :: proc ( shape : [ ]uint, forward : bool,
	data_in, data_out : [ ]complex128, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	stride := strides_row_major( shape, size_of( complex128 ) )
	defer delete( stride )
	axes := to_axes( shape )
	defer delete( axes )
	return c2c_f64_proc( shape, axes, stride, stride, forward, data_in, data_out, fct, nthreads )
}

// -- real-to-complex (single axis) ------------------------------------------

r2c :: proc{ r2c_single_f32, r2c_single_f64 }

r2c_single_f32 :: proc ( shape_in : [ ]uint, stride_in, stride_out : [ ]int,
	axis: uint, forward: bool,
	data_in : [ ]f32, data_out : [ ]complex64, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{
	if len( shape_in ) == 0 do return .Empty_Shape
	rc := r2c_f32( c.size_t( len( shape_in ) ), sptr( shape_in ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( axis ), 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

r2c_single_f64 :: proc ( shape_in : [ ]uint, stride_in, stride_out : [ ]int,
	axis : uint, forward : bool,
	data_in : [ ]f64, data_out : [ ]complex128, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape_in ) == 0 do return .Empty_Shape
	rc := r2c_f64( c.size_t( len( shape_in ) ), sptr( shape_in ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( axis ), 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// Convenience : contiguous 1D r2c, output length must be n/2+1.
r2c_simple :: proc{ r2c_simple_f32, r2c_simple_f64 }

r2c_simple_f32 :: proc ( shape: [ ]uint, forward : bool,
	data_in : [ ]f32, data_out : [ ]complex64, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	stride_in := strides_row_major( shape, size_of( f32 ) )
	defer delete( stride_in )
	out_shape := make( [ ]uint, len( shape ) )
	defer delete( out_shape )
	copy( out_shape, shape )
	out_shape[ len( out_shape ) - 1 ] = shape[ len( shape ) - 1 ] / 2 + 1
	stride_out := strides_row_major( out_shape, size_of( complex64 ) )
	defer delete( stride_out )
	return r2c_single_f32( shape, stride_in, stride_out, uint( len( shape ) - 1 ), forward,
		data_in, data_out, fct, nthreads )
}

r2c_simple_f64 :: proc ( shape : [ ]uint, forward : bool,
	data_in : [ ]f64, data_out : [ ]complex128, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	stride_in := strides_row_major( shape, size_of( f64 ) )
	defer delete( stride_in )
	out_shape := make( [ ]uint, len( shape ) )
	defer delete( out_shape )
	copy( out_shape, shape )
	out_shape[ len( out_shape ) - 1 ] = shape[ len( shape ) - 1 ] / 2 + 1
	stride_out := strides_row_major( out_shape, size_of( complex128 ) )
	defer delete( stride_out )
	return r2c_single_f64( shape, stride_in, stride_out, uint( len( shape ) - 1 ), forward,
		data_in, data_out, fct, nthreads )
}

// -- complex-to-real ( single axis ) ------------------------------------------

c2r :: proc{ c2r_single_f32, c2r_single_f64 }

c2r_single_f32 :: proc ( shape_out : [ ]uint, stride_in, stride_out : [ ]int,
	axis : uint, forward : bool,
	data_in : [ ]complex64, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{
	if len( shape_out ) == 0 do return .Empty_Shape
	rc := c2r_f32( c.size_t( len( shape_out ) ), sptr( shape_out ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( axis ), 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

c2r_single_f64 :: proc ( shape_out : [ ]uint, stride_in, stride_out : [ ]int,
	axis : uint, forward : bool,
	data_in : [ ]complex128, data_out : [ ]f64, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{
	if len( shape_out ) == 0 do return .Empty_Shape
	rc := c2r_f64( c.size_t( len( shape_out ) ), sptr( shape_out ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( axis ), 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

c2r_simple :: proc{ c2r_simple_f32, c2r_simple_f64 }

c2r_simple_f32 :: proc ( shape_out : [ ]uint, forward : bool,
	data_in : [ ]complex64, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	stride_out := strides_row_major( shape_out, size_of( f32 ) )
	defer delete( stride_out )
	in_shape := make( [ ]uint, len( shape_out ) )
	defer delete( in_shape )
	copy( in_shape, shape_out )
	in_shape[ len( in_shape ) - 1 ] = shape_out[ len( shape_out ) - 1 ] / 2 + 1
	stride_in := strides_row_major( in_shape, size_of( complex64 ) )
	defer delete( stride_in )
	return c2r_single_f32( shape_out, stride_in, stride_out, uint( len( shape_out ) - 1 ), forward,
		data_in, data_out, fct, nthreads )
}

c2r_simple_f64 :: proc ( shape_out : [ ]uint, forward : bool,
	data_in : [ ]complex128, data_out: [ ]f64, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	stride_out := strides_row_major( shape_out, size_of( f64 ) )
	defer delete( stride_out )
	in_shape := make( [ ]uint, len( shape_out ) )
	defer delete( in_shape )
	copy( in_shape, shape_out )
	in_shape[ len( in_shape ) - 1 ] = shape_out[ len( shape_out ) - 1 ] / 2 + 1
	stride_in := strides_row_major( in_shape, size_of( complex128 ) )
	defer delete( stride_in )
	return c2r_single_f64( shape_out, stride_in, stride_out, uint( len( shape_out ) - 1 ), forward,
		data_in, data_out, fct, nthreads )
}

// -- r2r FFTPACK half-complex -----------------------------------------------

r2r_fftpack :: proc{ r2r_fftpack_f32_proc, r2r_fftpack_f64_proc }

r2r_fftpack_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	real2hermitian, forward : bool,
	data_in, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_fftpack_f32( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ),
		1 if real2hermitian else 0, 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

r2r_fftpack_f64_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	real2hermitian, forward : bool,
	data_in, data_out : [ ]f64, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_fftpack_f64( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ),
		1 if real2hermitian else 0, 1 if forward else 0,
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// -- r2r separable Hartley --------------------------------------------------

r2r_separable_hartley :: proc{ r2r_separable_hartley_f32_proc, r2r_separable_hartley_f64_proc }

r2r_separable_hartley_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_separable_hartley_f32( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

r2r_separable_hartley_f64_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f64, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_separable_hartley_f64( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// -- r2r genuine Hartley ----------------------------------------------------

r2r_genuine_hartley :: proc{ r2r_genuine_hartley_f32_proc, r2r_genuine_hartley_f64_proc }

r2r_genuine_hartley_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_genuine_hartley_f32( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

r2r_genuine_hartley_f64_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f64, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_genuine_hartley_f64( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// -- r2r separable FHT ------------------------------------------------------

r2r_separable_fht :: proc{ r2r_separable_fht_f32_proc, r2r_separable_fht_f64_proc }

r2r_separable_fht_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_separable_fht_f32( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

r2r_separable_fht_f64_proc :: proc( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f64, fct : f64 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_separable_fht_f64( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// -- r2r genuine FHT --------------------------------------------------------

r2r_genuine_fht :: proc{ r2r_genuine_fht_f32_proc, r2r_genuine_fht_f64_proc }

r2r_genuine_fht_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f32, fct : f32 = 1, nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_genuine_fht_f32( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

r2r_genuine_fht_f64_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	data_in, data_out : [ ]f64, fct : f64 = 1, nthreads: uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := r2r_genuine_fht_f64( c.size_t( len( shape ) ), sptr( shape ),
		pptr( stride_in ), pptr( stride_out ), c.size_t( len( axes ) ), sptr( axes ),
		raw_data( data_in ), raw_data( data_out ), fct, c.size_t( nthreads ) )
	return Error( rc )
}

// -- DCT / DST --------------------------------------------------------------

dct :: proc{ dct_f32_proc, dct_f64_proc }
dst :: proc{ dst_f32_proc, dst_f64_proc }

dct_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	type : int, data_in, data_out : [ ]f32, fct : f32 = 1, ortho : bool = false,
	nthreads : uint = 1 ) -> Error
{
	if len( shape ) == 0 do return .Empty_Shape
	rc := dct_f32( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ),
		c.int( type ), raw_data( data_in ), raw_data( data_out ),
		fct, 1 if ortho else 0, c.size_t( nthreads ) )
	return Error( rc )
}

dct_f64_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	type : int, data_in, data_out : [ ]f64, fct : f64 = 1, ortho : bool = false,
	nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := dct_f64( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ),
		c.int( type ), raw_data( data_in ), raw_data( data_out ),
		fct, 1 if ortho else 0, c.size_t( nthreads ) )
	return Error( rc )
}

dst_f32_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out: [ ]int,
	type : int, data_in, data_out : [ ]f32, fct : f32 = 1, ortho : bool = false,
	nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := dst_f32( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ),
		c.int( type ), raw_data( data_in ), raw_data( data_out ),
		fct, 1 if ortho else 0, c.size_t( nthreads ) )
	return Error( rc )
}

dst_f64_proc :: proc ( shape, axes : [ ]uint, stride_in, stride_out : [ ]int,
	type : int, data_in, data_out : [ ]f64, fct : f64 = 1, ortho : bool = false,
	nthreads : uint = 1 ) -> Error
{

	if len( shape ) == 0 do return .Empty_Shape
	rc := dst_f64( c.size_t( len( shape ) ), sptr( shape ), pptr( stride_in ), pptr( stride_out ),
		c.size_t( len( axes ) ), sptr( axes ),
		c.int( type ), raw_data( data_in ), raw_data( data_out ),
		fct, 1 if ortho else 0, c.size_t( nthreads ) )
	return Error( rc )
}
