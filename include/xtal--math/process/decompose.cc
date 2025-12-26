#pragma once
#include "./any.cc"





#include "./decompose.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("decompose")
{
	using T_fit = bond::fit<>;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;

	auto mt19937_f = typename T_fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("magnum edition")
	{
		T_alpha alpha{-0.123};
		T_aphex aphex{-0.123, 0.456};

		T_alpha const alpha_sgn{-1    };
		T_alpha const alpha_mgn{ 0.123};

		T_aphex const aphex_sgn{-0.2604290310426322085923800386808580e0L, 0.9654929931336608817105116031598300e0L};
		T_aphex const aphex_mgn{ 0.4722975756871932162539451383054256e0L, 0.0000000000000000000000000000000000e0L};

		TRUE_(check_f<-1>(decompose_t<unsigned>::method_e(alpha), alpha_mgn));
		TRUE_(check_f<-1>(decompose_t<unsigned>::method_e(aphex), aphex_mgn));

		TRUE_(check_f<-1>(alpha, alpha_sgn));
		TRUE_(check_f<-1>(aphex, aphex_sgn));

		TRUE_(same_q<float, typename bond::fit<float>::alpha_type>);
		TRUE_(same_q<float, typename bond::fit<float>::alpha_type>);

		TRUE_(same_q<float, decltype(bond::fit<float>::alpha_1)>);
		TRUE_(same_q<float, decltype(bond::fit<float>::alpha_1)>);

		TRUE_(same_q<float, decltype(std::copysign(bond::fit<float>::alpha_1, +1.F))>);
		TRUE_(same_q<float, decltype(std::copysign(bond::fit<float>::alpha_1, -1.F))>);

		TRUE_(same_q<float, decltype(_xtd::copysign(bond::fit<float>::alpha_1, +1.F))>);
		TRUE_(same_q<float, decltype(_xtd::copysign(bond::fit<float>::alpha_1, -1.F))>);

		TRUE_(same_q<float, decltype(decompose_f<unsigned>(+1.F))>);
		TRUE_(same_q<float, decltype(decompose_f<unsigned>(-1.F))>);

	};
	TRY_("signum edition")
	{
		T_alpha alpha{-0.123};
		T_aphex aphex{-0.123, 0.456};

		T_alpha const alpha_sgn{-1    };
		T_alpha const alpha_mgn{ 0.123};

		T_aphex const aphex_sgn{-0.2604290310426322085923800386808580e0L, 0.9654929931336608817105116031598300e0L};
		T_aphex const aphex_mgn{ 0.4722975756871932162539451383054256e0L, 0.0000000000000000000000000000000000e0L};

		TRUE_(check_f<-1>(decompose_t<signed>::method_e(alpha), alpha_sgn));
		TRUE_(check_f<-1>(decompose_t<signed>::method_e(aphex), aphex_sgn));

		TRUE_(check_f<-1>(alpha, alpha_mgn));
		TRUE_(check_f<-1>(aphex, aphex_mgn));

	};
	TRY_("assigned_f( 1)")
	{
		TRUE_( 1. == decompose_f<signed>(T_sigma{ 1}));
		TRUE_( 1. == decompose_f<signed>(T_delta{ 1}));
		TRUE_( 1. == decompose_f<signed>(T_alpha{ 1}));
	//	TRUE_( 1. == decompose_f<signed>(true));

	};
	TRY_("assigned_f( 0)")
	{
		TRUE_(-1. == decompose_f<signed>(T_sigma{ 0}));
		TRUE_(-1. == decompose_f<signed>(T_delta{ 0}));
	//	TRUE_(-1. == decompose_f<signed>(T_alpha{-0}));
	//	TRUE_(-1. == decompose_f<signed>(false));

	};
	TRY_("assigned_f(-1)")
	{
		TRUE_(-1. == decompose_f<signed>(T_delta{-1}));
		TRUE_(-1. == decompose_f<signed>(T_alpha{-1}));

	};
}
////////////////////////////////////////////////////////////////////////////////
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
