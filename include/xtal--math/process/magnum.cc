#pragma once
#include "./any.cc"
#include "./magnum.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("magnum")
{
	using T_op = bond::operating;
	using T_sigma = typename T_op::sigma_type;
	using T_delta = typename T_op::delta_type;
	using T_alpha = typename T_op::alpha_type;
	using T_aphex = typename T_op::aphex_type;

	auto mt19937_f = typename T_op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("magnum edition")
	{
		T_alpha alpha{-0.123};
		T_aphex aphex{-0.123, 0.456};

		T_alpha const alpha_sgn{-1    };
		T_alpha const alpha_mgn{ 0.123};

		T_aphex const aphex_sgn{-0.2604290310426322085923800386808580e0L, 0.9654929931336608817105116031598300e0L};
		T_aphex const aphex_mgn{ 0.4722975756871932162539451383054256e0L, 0.0000000000000000000000000000000000e0L};

		TRUE_(check_f<-1>(magnum_t<>::edit(alpha), alpha_mgn));
		TRUE_(check_f<-1>(magnum_t<>::edit(aphex), aphex_mgn));

		TRUE_(check_f<-1>(alpha, alpha_sgn));
		TRUE_(check_f<-1>(aphex, aphex_sgn));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
