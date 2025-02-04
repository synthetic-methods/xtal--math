#pragma once
#include "./any.cc"
#include "./disunity.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("disunity")
{
	using _fit = bond::fit<>;
	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	TRY_("task")
	{
		TRUE_(check_f<22>(disunity_f<-1, -1>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 2.0), disunity_f<(+1), 1>(T_aphex {2, 3})));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 0.5), disunity_f<(-1), 1>(T_aphex {2, 3})));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
