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
	using _fix = bond::fixture<>;
	using U_delta = typename _fix::delta_type;
	using U_sigma = typename _fix::sigma_type;
	using U_alpha = typename _fix::alpha_type;
	using U_aphex = typename _fix::aphex_type;
	
	using W_alpha = atom::couple_t<U_alpha[2]>;
	using W_aphex = atom::couple_t<U_aphex[2]>;

	TRY_("term")
	{
		TRUE_(factorial_f<5>() == 120);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
