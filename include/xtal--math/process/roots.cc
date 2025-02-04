#pragma once
#include "./any.cc"
#include "./roots.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("roots")
{
	using _fit = bond::fit<>;
	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	TRY_("evaluation")
	{
		TRUE_(check_f<-1>(get<0>(roots_f< 4>(0.5)), root_f< 4>(0.5)));
		TRUE_(check_f<-1>(get<1>(roots_f< 4>(0.5)), root_f<-4>(0.5)));

		TRUE_(check_f<-1>(get<0>(roots_f<-4>(0.5)), root_f<-4>(0.5)));
		TRUE_(check_f<-1>(get<1>(roots_f<-4>(0.5)), root_f< 4>(0.5)));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
