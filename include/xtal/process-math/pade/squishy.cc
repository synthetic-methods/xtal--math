#pragma once
#include "./any.cc"
#include "./squishy.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("squishy")
{
	using _op = bond::operating;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;

	TRY_("task")
	{
		TRUE_(check_f<22>(squishy_f<-1, -1>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 2.0), squishy_f<(+1), 1>(T_aphex {2, 3})));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 0.5), squishy_f<(-1), 1>(T_aphex {2, 3})));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
