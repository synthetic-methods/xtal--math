#pragma once
#include "./any.cc"





#include "./roots.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("roots")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

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
