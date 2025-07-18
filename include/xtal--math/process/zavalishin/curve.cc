#pragma once
#include "./any.cc"
#include "./curve.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("curve")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("curve cardinal")
	{
		TRUE_(check_f<-32>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-3,-1>::template method_f< 4>(0.5, 0.5)));
		TRUE_(check_f<-35>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-3,-1>::template method_f< 3>(0.5, 0.5)));
		TRUE_(check_f<-39>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-3,-1>::template method_f< 2>(0.5, 0.5)));
		TRUE_(check_f<-44>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-3,-1>::template method_f< 1>(0.5, 0.5)));
		TRUE_(check_f<-41>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-2,-1>::template method_f< 4>(0.5, 0.5)));
		TRUE_(check_f<-41>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-2,-1>::template method_f< 3>(0.5, 0.5)));
		TRUE_(check_f<-45>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-2,-1>::template method_f< 2>(0.5, 0.5)));
		TRUE_(check_f<-45>(curve_t<-3,-1>::template method_f<~0>(0.5, 0.5), curve_t<-2,-1>::template method_f< 1>(0.5, 0.5)));

		TRUE_(check_f<-36>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 3,-1>::template method_f< 4>(0.5, 0.5)));
		TRUE_(check_f<-38>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 3,-1>::template method_f< 3>(0.5, 0.5)));
		TRUE_(check_f<-41>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 3,-1>::template method_f< 2>(0.5, 0.5)));
		TRUE_(check_f<-44>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 3,-1>::template method_f< 1>(0.5, 0.5)));
		TRUE_(check_f<-38>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 2,-1>::template method_f< 4>(0.5, 0.5)));
		TRUE_(check_f<-41>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 2,-1>::template method_f< 3>(0.5, 0.5)));
		TRUE_(check_f<-44>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 2,-1>::template method_f< 2>(0.5, 0.5)));
		TRUE_(check_f<-45>(curve_t< 3,-1>::template method_f<~0>(0.5, 0.5), curve_t< 2,-1>::template method_f< 1>(0.5, 0.5)));

	}
	TRY_("curve isomorphism")
	{
		TRUE_(check_f<-40>(0.5, curve_t<-3>::template method_f< 2>(curve_t< 3>::template method_f< 2>(0.5, 0.5), 0.5)));
		TRUE_(check_f<-40>(0.5, curve_t< 3>::template method_f< 2>(curve_t<-3>::template method_f< 2>(0.5, 0.5), 0.5)));

		TRUE_(check_f<- 7>(0.5, curve_t<-2>::template method_f< 4>(curve_t< 2>::template method_f< 4>(0.5, 0.5), 0.5)));
		TRUE_(check_f<- 7>(0.5, curve_t< 2>::template method_f< 4>(curve_t<-2>::template method_f< 4>(0.5, 0.5), 0.5)));

	}
	/**/
	TRY_("curve evaluation")
	{
		TRUE_(check_f<-2>(0.6180339887498949,  curve_t<-2, -1>::template method_f<1, 1>(-1.0)));
		TRUE_(check_f<-2>(0.4142135623730949,  curve_t<-2, -1>::template method_f<1, 1>(-2.0)));
		TRUE_(check_f<-2>(0.3027756377319948,  curve_t<-2, -1>::template method_f<1, 1>(-3.0)));
		TRUE_(check_f<-2>(0.2360679774997898,  curve_t<-2, -1>::template method_f<1, 1>(-4.0)));
		TRUE_(check_f<-4>(0.1925824035672532,  curve_t<-2, -1>::template method_f<1, 1>(-5.0)));
		TRUE_(check_f<-2>(0.1622776601683800,  curve_t<-2, -1>::template method_f<1, 1>(-6.0)));
		TRUE_(check_f<-2>(0.1400549446402577,  curve_t<-2, -1>::template method_f<1, 1>(-7.0)));
		TRUE_(check_f<-2>(0.1231056256176600,  curve_t<-2, -1>::template method_f<1, 1>(-8.0)));
		TRUE_(check_f<-1>(0.1097722286464435,  curve_t<-2, -1>::template method_f<1, 1>(-9.0)));

		TRUE_(check_f<-5>(1.6180339887498949,  curve_t<-2, -1>::template method_f<1, 1>( 1.0)));
		TRUE_(check_f<-5>(2.4142135623730949,  curve_t<-2, -1>::template method_f<1, 1>( 2.0)));
		TRUE_(check_f<-5>(3.3027756377319948,  curve_t<-2, -1>::template method_f<1, 1>( 3.0)));
		TRUE_(check_f<-5>(4.2360679774997898,  curve_t<-2, -1>::template method_f<1, 1>( 4.0)));
		TRUE_(check_f<-5>(5.1925824035672532,  curve_t<-2, -1>::template method_f<1, 1>( 5.0)));
		TRUE_(check_f<-5>(6.1622776601683800,  curve_t<-2, -1>::template method_f<1, 1>( 6.0)));
		TRUE_(check_f<-5>(7.1400549446402577,  curve_t<-2, -1>::template method_f<1, 1>( 7.0)));
		TRUE_(check_f<-5>(8.1231056256176600,  curve_t<-2, -1>::template method_f<1, 1>( 8.0)));
		TRUE_(check_f<-5>(9.1097722286464435,  curve_t<-2, -1>::template method_f<1, 1>( 9.0)));

		TRUE_(check_f<-9>(1.6180339887498950,  curve_t<-2, -1>::template method_f<1, 1>( 1.0)));
		TRUE_(check_f<-9>(2.4142135623730950,  curve_t<-2, -1>::template method_f<1, 1>( 2.0)));
		TRUE_(check_f<-9>(3.3027756377319952,  curve_t<-2, -1>::template method_f<1, 1>( 3.0)));
		TRUE_(check_f<-9>(4.2360679774997900,  curve_t<-2, -1>::template method_f<1, 1>( 4.0)));
		TRUE_(check_f<-9>(5.1925824035672520,  curve_t<-2, -1>::template method_f<1, 1>( 5.0)));
		TRUE_(check_f<-9>(6.1622776601683790,  curve_t<-2, -1>::template method_f<1, 1>( 6.0)));
		TRUE_(check_f<-9>(7.1400549446402590,  curve_t<-2, -1>::template method_f<1, 1>( 7.0)));
		TRUE_(check_f<-9>(8.1231056256176600,  curve_t<-2, -1>::template method_f<1, 1>( 8.0)));
		TRUE_(check_f<-9>(9.1097722286464420,  curve_t<-2, -1>::template method_f<1, 1>( 9.0)));

		TRUE_(check_f<-9>(0.4438171626239681, -curve_t< 3,  0>::template method_f<~0>(-0.5,-0.5)));
		TRUE_(check_f< 9>(0.4438171626239681, -curve_t< 3,  0>::template method_f< 9>(-0.5,-0.5)));


		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f<~0>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 9>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 8>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 7>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 6>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 5>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 4>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 3>( 0.5, 0.5)));
		TRUE_(check_f< 8>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 2>( 0.5, 0.5)));
		TRUE_(check_f< 7>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 1>( 0.5, 0.5)));
		TRUE_(check_f< 2>(0.4438171626239681,  curve_t< 3,  0>::template method_f< 0>( 0.5, 0.5)));
	
		TRUE_(check_f<20>(0.5739597111931870, -curve_t<-3,  0>::template method_f<~0>(-0.5, -0.5)));
		TRUE_(check_f<20>(0.5739597111931870, -curve_t<-3,  0>::template method_f< 9>(-0.5, -0.5)));
		TRUE_(check_f<13>(0.5739597111931870, -curve_t<-3,  0>::template method_f< 2>(-0.5, -0.5)));
		TRUE_(check_f< 8>(0.5739597111931870, -curve_t<-3,  0>::template method_f< 1>(-0.5, -0.5)));

		TRUE_(check_f<20>(0.57395971119318700, curve_t<-3,  0>::template method_f<~0>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, curve_t<-3,  0>::template method_f< 9>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, curve_t<-3,  0>::template method_f< 8>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, curve_t<-3,  0>::template method_f< 7>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, curve_t<-3,  0>::template method_f< 6>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, curve_t<-3,  0>::template method_f< 5>( 0.5,  0.5)));
		TRUE_(check_f<19>(0.57395971119318700, curve_t<-3,  0>::template method_f< 4>( 0.5,  0.5)));
		TRUE_(check_f<17>(0.57395971119318700, curve_t<-3,  0>::template method_f< 3>( 0.5,  0.5)));
		TRUE_(check_f<13>(0.57395971119318700, curve_t<-3,  0>::template method_f< 2>( 0.5,  0.5)));
		TRUE_(check_f< 8>(0.57395971119318700, curve_t<-3,  0>::template method_f< 1>( 0.5,  0.5)));
		TRUE_(check_f< 2>(0.57395971119318700, curve_t<-3,  0>::template method_f< 0>( 0.5,  0.5)));

		TRUE_(check_f<15>(exp(0.5) - one, curve_t<-3>::template method_f<9>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, curve_t<-3>::template method_f<8>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, curve_t<-3>::template method_f<7>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, curve_t<-3>::template method_f<6>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, curve_t<-3>::template method_f<5>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, curve_t<-3>::template method_f<4>(0.5, 1.0)));
		TRUE_(check_f<13>(exp(0.5) - one, curve_t<-3>::template method_f<3>(0.5, 1.0)));
		TRUE_(check_f<11>(exp(0.5) - one, curve_t<-3>::template method_f<2>(0.5, 1.0)));
		TRUE_(check_f< 5>(exp(0.5) - one, curve_t<-3>::template method_f<1>(0.5, 1.0)));
		TRUE_(check_f< 1>(exp(0.5) - one, curve_t<-3>::template method_f<0>(0.5, 1.0)));

		TRUE_(check_f<36>(sinh(0.5), curve_t<-3>::template method_f<9>(0.5, 0.0)));
		TRUE_(check_f<33>(sinh(0.5), curve_t<-3>::template method_f<8>(0.5, 0.0)));
		TRUE_(check_f<27>(sinh(0.5), curve_t<-3>::template method_f<7>(0.5, 0.0)));
		TRUE_(check_f<24>(sinh(0.5), curve_t<-3>::template method_f<6>(0.5, 0.0)));
		TRUE_(check_f<22>(sinh(0.5), curve_t<-3>::template method_f<5>(0.5, 0.0)));
		TRUE_(check_f<18>(sinh(0.5), curve_t<-3>::template method_f<4>(0.5, 0.0)));
		TRUE_(check_f<15>(sinh(0.5), curve_t<-3>::template method_f<3>(0.5, 0.0)));
		TRUE_(check_f<11>(sinh(0.5), curve_t<-3>::template method_f<2>(0.5, 0.0)));
		TRUE_(check_f< 9>(sinh(0.5), curve_t<-3>::template method_f<1>(0.5, 0.0)));
		TRUE_(check_f< 2>(sinh(0.5), curve_t<-3>::template method_f<0>(0.5, 0.0)));

		TRUE_(check_f<31>(sinh(1.0), curve_t<-3>::template method_f<9>(1.0, 0.0)));
		TRUE_(check_f<28>(sinh(1.0), curve_t<-3>::template method_f<8>(1.0, 0.0)));
		TRUE_(check_f<22>(sinh(1.0), curve_t<-3>::template method_f<7>(1.0, 0.0)));
		TRUE_(check_f<18>(sinh(1.0), curve_t<-3>::template method_f<6>(1.0, 0.0)));
		TRUE_(check_f<16>(sinh(1.0), curve_t<-3>::template method_f<5>(1.0, 0.0)));
		TRUE_(check_f<13>(sinh(1.0), curve_t<-3>::template method_f<4>(1.0, 0.0)));
		TRUE_(check_f< 9>(sinh(1.0), curve_t<-3>::template method_f<3>(1.0, 0.0)));
		TRUE_(check_f< 6>(sinh(1.0), curve_t<-3>::template method_f<2>(1.0, 0.0)));
		TRUE_(check_f< 4>(sinh(1.0), curve_t<-3>::template method_f<1>(1.0, 0.0)));
		TRUE_(check_f< 1>(sinh(1.0), curve_t<-3>::template method_f<0>(1.0, 0.0)));

	};
	/***/

	EST_("curve_t<-3>::method_f<~0, 0> (exp(#) - one)")
	{
		T_alpha o{1};
		for (T_sigma i = 0x100; ~--i;) {
			//\
			o *= curve_t<-3, -0>::template method_f<~0, 0>(_fit::mantissa_f(mt19937_f));
			o *= exp(_fit::mantissa_f(mt19937_f)) - one;
		}
		return o;

	};
	/**/
	EST_("curve_t<-3>::method_f< 8, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= curve_t<-3, -0>::template method_f< 8, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("curve_t<-3>::method_f< 4, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= curve_t<-3, -0>::template method_f< 4, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("curve_t<-3>::method_f< 3, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= curve_t<-3, -0>::template method_f< 3, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("curve_t<-3>::method_f< 2, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= curve_t<-3, -0>::template method_f< 2, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("curve_t<-3>::method_f< 1, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= curve_t<-3, -0>::template method_f< 1, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	/***/

}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
