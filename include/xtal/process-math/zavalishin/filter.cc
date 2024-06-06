#pragma once
#include "./any.cc"
#include "./filter.hh"// testing...






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("whatever")
{
	using op = bond::operating;

	using T_sigma = typename op::sigma_t;
	using T_delta = typename op::delta_t;
	using T_alpha = typename op::alpha_t;
	using T_aphex = typename op::aphex_t;
	XTAL_LET_(T_alpha) pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	using U_phi = algebra::d_::circular_t<T_alpha[2]>;

	auto mt19937_f = typename op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(true);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
