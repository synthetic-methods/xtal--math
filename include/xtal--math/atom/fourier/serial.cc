#pragma once
#include "./any.cc"
#include "./serial.hh"// testing...





XTAL_ENV_(push)
namespace xtal::atom::math::fourier::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//atic_assert(_xtd::trivially_initializable<serial_t<float[2]>>);
static_assert(_xtd::trivially_destructible <serial_t<float[2]>>);
static_assert(_xtd::trivially_copyable     <serial_t<float[2]>>);
static_assert(_xtd::trivially_movable      <serial_t<float[2]>>);
//atic_assert(                     atomic_q<serial_t<float[2]>>);


////////////////////////////////////////////////////////////////////////////////

TAG_("atom", "serial")
{
	using _fit = bond::fit<>;
	using T_delta = typename _fit::delta_type;
	using T_sigma = typename _fit::sigma_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using D1 = serial_t<int[1]>;
	using D2 = serial_t<int[2]>;
	using D3 = serial_t<int[3]>;
	using D4 = serial_t<int[4]>;
	
	TRY_("partial construction")
	{
		D4 d4{1000};

		TRUE_(d4 == D4{1000, 0, 0, 0});

	}
	TRY_("multiplication")
	{
		TRUE_(D2 {          10, 1} * D2 {          20, 2} == D2 {                  200,   40});
		TRUE_(D3 {     100, 10, 1} * D3 {     200, 20, 2} == D3 {         20000,  4000,  600});
		TRUE_(D4{1000, 100, 10, 1} * D4{2000, 200, 20, 2} == D4{2000000, 400000, 60000, 8000});

	}
	TRY_("integration")
	{
		D4 d4{1000, 100, 10, 1};

		TRUE_(++d4 == D4{1100, 110, 11, 1});
		TRUE_(++d4 == D4{1210, 121, 12, 1});
		TRUE_(++d4 == D4{1331, 133, 13, 1});
		TRUE_(--d4 == D4{1210, 121, 12, 1});
		TRUE_(--d4 == D4{1100, 110, 11, 1});

	}
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
