#pragma once
#include "./any.cc"
#include "./square.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("math")
{
	using _fit = bond::fit<>;
	using U_delta = typename _fit::delta_type;
	using U_sigma = typename _fit::sigma_type;
	using U_alpha = typename _fit::alpha_type;
	using U_aphex = typename _fit::aphex_type;
	
	using W_alpha = atom::couple_t<U_alpha[2]>;
	using W_aphex = atom::couple_t<U_aphex[2]>;

	TRY_("square (singular)")
	{
		TRUE_(square_f(atom::couple_f(2, 3)) == atom::couple_f(4, 9));
		TRUE_(square_f(W_alpha{2, 3}) == W_alpha{4, 9});

		U_aphex u01{0, 1}; U_aphex w01 = u01*u01;
		U_aphex u23{2, 3}; U_aphex w23 = u23*u23;

		TRUE_(square_f(atom::couple_f(u01, u23)) == atom::couple_f(w01, w23));
		TRUE_(square_f(W_aphex{u01, u23}) == W_aphex{w01, w23});

		U_aphex u{ 0.17364817768602775, 0.98480775300884071};
		U_aphex w{-0.93969262077264359, 0.34202014336211378};
		TRUE_(check_f<-2>(w, square_f(u)));

	}
	TRY_("square (plural)")
	{
		TRUE_(square_f(1.0, 2.0, 3.0) == 1.0*1.0 + 2.0*2.0 + 3.0*3.0);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
