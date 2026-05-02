#pragma once
#include "./any.cc"





#include "./unify.hh"
XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("unify")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	TRY_("task")
	{
		TRUE_(check_f<22>(unify_f<-1, -1>(pow(U_aphex {2, 3}, -2.0)), U_aphex {2, 3}));
		TRUE_(check_f<22>(pow(U_aphex {2, 3}, 2.0), unify_f<(+1), 1>(U_aphex {2, 3})));
		TRUE_(check_f<22>(pow(U_aphex {2, 3}, 0.5), unify_f<(-1), 1>(U_aphex {2, 3})));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
