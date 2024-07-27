#pragma once
#include "./any.cc"
#include "./root.hh"// testing...






XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("root")
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
		TRUE_(check_f<22>(root_f<-2>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 0.5), root_f< 2>(T_aphex {2, 3})));

	}
	TRY_("punctured evaluation")
	{
		TRUE_(1.0 + root_f<-1, 0>(0.0) == 1.0 + root_f<-1, 0>(0.0));
		TRUE_(1.0 + root_f<-2, 0>(0.0) == 1.0 + root_f<-2, 0>(0.0));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
