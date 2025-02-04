#pragma once
#include "./any.cc"
#include "./term.hh"// testing...





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

	TRY_("term")
	{
		TRUE_(term_f(one, 2.0, 3.0) == 7.0);
		TRUE_(term_f(one, 2.0, 3.0, 5.0, 7.0) == 211.0);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
