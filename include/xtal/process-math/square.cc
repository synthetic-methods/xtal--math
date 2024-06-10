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
	using _op = bond::operating;
	using T_sigma = typename _op::sigma_t;
	using T_delta = typename _op::delta_t;
	using T_alpha = typename _op::alpha_t;
	using T_aphex = typename _op::aphex_t;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

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
