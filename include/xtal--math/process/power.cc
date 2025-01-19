#pragma once
#include "./any.cc"
#include "./power.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("power")
{
	using _fix = bond::fixture<>;
	using T_sigma = typename _fix::sigma_type;
	using T_delta = typename _fix::delta_type;
	using T_alpha = typename _fix::alpha_type;
	using T_aphex = typename _fix::aphex_type;

	auto mt19937_f = typename _fix::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("power evaluation")
	{
		TRUE_(power_t< 7>::function(3.) == 2187.);
		TRUE_(power_f< 7>          (3.) == 2187.);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
