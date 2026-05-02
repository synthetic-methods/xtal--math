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
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	auto mt19937_f = typename U_fit::MT19937();
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
