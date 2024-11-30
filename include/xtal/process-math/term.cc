#pragma once
#include "./any.cc"
#include "./term.hh"// testing...

#include "./pade/unity.hh"




XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("math")
{
	using _op = bond::operating;
	using U_delta = typename _op::delta_type;
	using U_sigma = typename _op::sigma_type;
	using U_alpha = typename _op::alpha_type;
	using U_aphex = typename _op::aphex_type;
	
	using W_alpha = algebra::scalar_t<U_alpha[2]>;
	using W_aphex = algebra::scalar_t<U_aphex[2]>;

	TRY_("term")
	{
		TRUE_(term_f(1.0, 2.0, 3.0) == 7.0);
		TRUE_(term_f(1.0, 2.0, 3.0, 5.0, 7.0) == 211.0);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
