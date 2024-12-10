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
	using _op = bond::operating;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

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
