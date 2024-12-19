#pragma once
#include "./any.cc"
#include "./sinoid.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::goomtrex::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("sinoid")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	static constexpr T_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("sinoid stuff")
	{
		TRUE_(true);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
