#pragma once
#include "./any.cc"





#include "./monomial.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("monomial")
{
	using _fit = bond::fit<>;
	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	auto mt19937_f = typename _fit::MT19937();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("monomial evaluation")
	{
		TRUE_(monomial_t< 7>{}.method(3.) == 2187.);
		TRUE_(monomial_f< 7>          (3.) == 2187.);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
