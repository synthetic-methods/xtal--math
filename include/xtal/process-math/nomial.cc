#pragma once
#include "./any.cc"
#include "./nomial.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("nomial")
{
	using _op = bond::operating;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(nomial_t< 7>::function(3.) == 2187.);
		TRUE_(nomial_f< 7>          (3.) == 2187.);

		TRUE_(check_f<-1>(nomial_t<-3>::function(3.), T_alpha{1.442249570307408382321638310780110L}));
		TRUE_(check_f<-1>(nomial_f<-3>          (3.), T_alpha{1.442249570307408382321638310780110L}));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
