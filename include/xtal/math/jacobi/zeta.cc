#pragma once
#include "./any.cc"
#include "./zeta.hh"// testing...






XTAL_ENV_(push)
namespace xtal::math::jacobi::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using re = bond::realized;

	using T_sigma = typename re::sigma_t;
	using T_delta = typename re::delta_t;
	using T_alpha = typename re::alpha_t;
	using T_aphex = typename re::aphex_t;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	using U_phi = algebra::differential::modular_t<T_alpha[2]>;
	using W_phi = _std::complex<U_phi>;

	auto mt19937_f = typename re::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		//\
		W_phi   constexpr t = { 0.22, 0.11};
		T_aphex constexpr t = { 0.22, 0.11};
		T_aphex constexpr s = { 0.22, 1.11};
		T_aphex constexpr r = {-0.150796, 0.106422};

		TRUE_(bond::computrim_f<16>(r) == bond::computrim_f<16>(zeta_f(t, s)));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
