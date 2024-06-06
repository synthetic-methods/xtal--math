#pragma once
#include "./any.cc"
#include "./square.hh"// testing...






XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("square")
{
	using op = bond::operating;
	using T_sigma = typename op::sigma_t;
	using T_delta = typename op::delta_t;
	using T_alpha = typename op::alpha_t;
	using T_aphex = typename op::aphex_t;
	XTAL_LET_(T_alpha) one =  1;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	TRY_("task")
	{
		using _std::pow;

		TRUE_(check_f<22>(square_f<-1, -1>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 2.0), square_f<(+1), 1>(T_aphex {2, 3})));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 0.5), square_f<(-1), 1>(T_aphex {2, 3})));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
