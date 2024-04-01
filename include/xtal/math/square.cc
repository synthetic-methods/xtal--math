#pragma once
#include "./any.cc"
#include "./square.ii"// testing...






XTAL_ENV_(push)
namespace xtal::math::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("square")
{
	using re = bond::realized;
	using T_sigma = typename re::sigma_t;
	using T_delta = typename re::delta_t;
	using T_alpha = typename re::alpha_t;
	using T_aphex = typename re::aphex_t;
	XTAL_LET_(T_alpha) one =  1;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	TRY_("task")
	{
		TRUE_(check_f<22>(_std::pow(T_aphex {1, 1}, 2.0), square_f<(+1), 1>(T_aphex {1, 1})));
		TRUE_(check_f<22>(_std::pow(T_aphex {1, 1}, 0.5), square_f<(-1), 1>(T_aphex {1, 1})));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
