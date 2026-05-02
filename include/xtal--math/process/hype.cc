#pragma once
#include "./any.cc"





#include "./hype.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("hype")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _fit::MT19937();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("hype cardinal")
	{
		TRUE_(check_f<-32>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-3,-1>{}.template method< 4>(0.5, 0.5)));
		TRUE_(check_f<-35>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-3,-1>{}.template method< 3>(0.5, 0.5)));
		TRUE_(check_f<-39>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-3,-1>{}.template method< 2>(0.5, 0.5)));
		TRUE_(check_f<-44>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-3,-1>{}.template method< 1>(0.5, 0.5)));
		TRUE_(check_f<-41>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-2,-1>{}.template method< 4>(0.5, 0.5)));
		TRUE_(check_f<-41>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-2,-1>{}.template method< 3>(0.5, 0.5)));
		TRUE_(check_f<-45>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-2,-1>{}.template method< 2>(0.5, 0.5)));
		TRUE_(check_f<-45>(hype_t<-3,-1>{}.template method<~0>(0.5, 0.5), hype_t<-2,-1>{}.template method< 1>(0.5, 0.5)));

		TRUE_(check_f<-36>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 3,-1>{}.template method< 4>(0.5, 0.5)));
		TRUE_(check_f<-38>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 3,-1>{}.template method< 3>(0.5, 0.5)));
		TRUE_(check_f<-41>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 3,-1>{}.template method< 2>(0.5, 0.5)));
		TRUE_(check_f<-44>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 3,-1>{}.template method< 1>(0.5, 0.5)));
		TRUE_(check_f<-38>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 2,-1>{}.template method< 4>(0.5, 0.5)));
		TRUE_(check_f<-41>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 2,-1>{}.template method< 3>(0.5, 0.5)));
		TRUE_(check_f<-44>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 2,-1>{}.template method< 2>(0.5, 0.5)));
		TRUE_(check_f<-45>(hype_t< 3,-1>{}.template method<~0>(0.5, 0.5), hype_t< 2,-1>{}.template method< 1>(0.5, 0.5)));

	}
	TRY_("hype isomorphism")
	{
		TRUE_(check_f<-40>(0.5, hype_t<-3>{}.template method< 2>(hype_t< 3>{}.template method< 2>(0.5, 0.5), 0.5)));
		TRUE_(check_f<-40>(0.5, hype_t< 3>{}.template method< 2>(hype_t<-3>{}.template method< 2>(0.5, 0.5), 0.5)));

		TRUE_(check_f<- 7>(0.5, hype_t<-2>{}.template method< 4>(hype_t< 2>{}.template method< 4>(0.5, 0.5), 0.5)));
		TRUE_(check_f<- 7>(0.5, hype_t< 2>{}.template method< 4>(hype_t<-2>{}.template method< 4>(0.5, 0.5), 0.5)));

	}
	/**/
	TRY_("hype evaluation")
	{
		TRUE_(check_f<-2>(0.6180339887498949,  hype_t<-2, -1>{}.template method<1, 1>(-1.0)));
		TRUE_(check_f<-2>(0.4142135623730949,  hype_t<-2, -1>{}.template method<1, 1>(-2.0)));
		TRUE_(check_f<-2>(0.3027756377319948,  hype_t<-2, -1>{}.template method<1, 1>(-3.0)));
		TRUE_(check_f<-2>(0.2360679774997898,  hype_t<-2, -1>{}.template method<1, 1>(-4.0)));
		TRUE_(check_f<-4>(0.1925824035672532,  hype_t<-2, -1>{}.template method<1, 1>(-5.0)));
		TRUE_(check_f<-2>(0.1622776601683800,  hype_t<-2, -1>{}.template method<1, 1>(-6.0)));
		TRUE_(check_f<-2>(0.1400549446402577,  hype_t<-2, -1>{}.template method<1, 1>(-7.0)));
		TRUE_(check_f<-2>(0.1231056256176600,  hype_t<-2, -1>{}.template method<1, 1>(-8.0)));
		TRUE_(check_f<-1>(0.1097722286464435,  hype_t<-2, -1>{}.template method<1, 1>(-9.0)));

		TRUE_(check_f<-5>(1.6180339887498949,  hype_t<-2, -1>{}.template method<1, 1>( 1.0)));
		TRUE_(check_f<-5>(2.4142135623730949,  hype_t<-2, -1>{}.template method<1, 1>( 2.0)));
		TRUE_(check_f<-5>(3.3027756377319948,  hype_t<-2, -1>{}.template method<1, 1>( 3.0)));
		TRUE_(check_f<-5>(4.2360679774997898,  hype_t<-2, -1>{}.template method<1, 1>( 4.0)));
		TRUE_(check_f<-5>(5.1925824035672532,  hype_t<-2, -1>{}.template method<1, 1>( 5.0)));
		TRUE_(check_f<-5>(6.1622776601683800,  hype_t<-2, -1>{}.template method<1, 1>( 6.0)));
		TRUE_(check_f<-5>(7.1400549446402577,  hype_t<-2, -1>{}.template method<1, 1>( 7.0)));
		TRUE_(check_f<-5>(8.1231056256176600,  hype_t<-2, -1>{}.template method<1, 1>( 8.0)));
		TRUE_(check_f<-5>(9.1097722286464435,  hype_t<-2, -1>{}.template method<1, 1>( 9.0)));

		TRUE_(check_f<-9>(1.6180339887498950,  hype_t<-2, -1>{}.template method<1, 1>( 1.0)));
		TRUE_(check_f<-9>(2.4142135623730950,  hype_t<-2, -1>{}.template method<1, 1>( 2.0)));
		TRUE_(check_f<-9>(3.3027756377319952,  hype_t<-2, -1>{}.template method<1, 1>( 3.0)));
		TRUE_(check_f<-9>(4.2360679774997900,  hype_t<-2, -1>{}.template method<1, 1>( 4.0)));
		TRUE_(check_f<-9>(5.1925824035672520,  hype_t<-2, -1>{}.template method<1, 1>( 5.0)));
		TRUE_(check_f<-9>(6.1622776601683790,  hype_t<-2, -1>{}.template method<1, 1>( 6.0)));
		TRUE_(check_f<-9>(7.1400549446402590,  hype_t<-2, -1>{}.template method<1, 1>( 7.0)));
		TRUE_(check_f<-9>(8.1231056256176600,  hype_t<-2, -1>{}.template method<1, 1>( 8.0)));
		TRUE_(check_f<-9>(9.1097722286464420,  hype_t<-2, -1>{}.template method<1, 1>( 9.0)));

		TRUE_(check_f<-9>(0.4438171626239681, -hype_t< 3,  0>{}.template method<~0>(-0.5,-0.5)));
		TRUE_(check_f< 9>(0.4438171626239681, -hype_t< 3,  0>{}.template method< 9>(-0.5,-0.5)));


		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method<~0>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 9>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 8>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 7>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 6>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 5>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 4>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 3>( 0.5, 0.5)));
		TRUE_(check_f< 8>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 2>( 0.5, 0.5)));
		TRUE_(check_f< 7>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 1>( 0.5, 0.5)));
		TRUE_(check_f< 2>(0.4438171626239681,  hype_t< 3,  0>{}.template method< 0>( 0.5, 0.5)));
	
		TRUE_(check_f<20>(0.5739597111931870, -hype_t<-3,  0>{}.template method<~0>(-0.5, -0.5)));
		TRUE_(check_f<20>(0.5739597111931870, -hype_t<-3,  0>{}.template method< 9>(-0.5, -0.5)));
		TRUE_(check_f<13>(0.5739597111931870, -hype_t<-3,  0>{}.template method< 2>(-0.5, -0.5)));
		TRUE_(check_f< 8>(0.5739597111931870, -hype_t<-3,  0>{}.template method< 1>(-0.5, -0.5)));

		TRUE_(check_f<20>(0.57395971119318700, hype_t<-3,  0>{}.template method<~0>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, hype_t<-3,  0>{}.template method< 9>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, hype_t<-3,  0>{}.template method< 8>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, hype_t<-3,  0>{}.template method< 7>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, hype_t<-3,  0>{}.template method< 6>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, hype_t<-3,  0>{}.template method< 5>( 0.5,  0.5)));
		TRUE_(check_f<19>(0.57395971119318700, hype_t<-3,  0>{}.template method< 4>( 0.5,  0.5)));
		TRUE_(check_f<17>(0.57395971119318700, hype_t<-3,  0>{}.template method< 3>( 0.5,  0.5)));
		TRUE_(check_f<13>(0.57395971119318700, hype_t<-3,  0>{}.template method< 2>( 0.5,  0.5)));
		TRUE_(check_f< 8>(0.57395971119318700, hype_t<-3,  0>{}.template method< 1>( 0.5,  0.5)));
		TRUE_(check_f< 2>(0.57395971119318700, hype_t<-3,  0>{}.template method< 0>( 0.5,  0.5)));

		TRUE_(check_f<15>(exp(0.5) - one, hype_t<-3>{}.template method<9>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, hype_t<-3>{}.template method<8>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, hype_t<-3>{}.template method<7>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, hype_t<-3>{}.template method<6>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, hype_t<-3>{}.template method<5>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, hype_t<-3>{}.template method<4>(0.5, 1.0)));
		TRUE_(check_f<13>(exp(0.5) - one, hype_t<-3>{}.template method<3>(0.5, 1.0)));
		TRUE_(check_f<11>(exp(0.5) - one, hype_t<-3>{}.template method<2>(0.5, 1.0)));
		TRUE_(check_f< 5>(exp(0.5) - one, hype_t<-3>{}.template method<1>(0.5, 1.0)));
		TRUE_(check_f< 1>(exp(0.5) - one, hype_t<-3>{}.template method<0>(0.5, 1.0)));

		TRUE_(check_f<36>(sinh(0.5), hype_t<-3>{}.template method<9>(0.5, 0.0)));
		TRUE_(check_f<33>(sinh(0.5), hype_t<-3>{}.template method<8>(0.5, 0.0)));
		TRUE_(check_f<27>(sinh(0.5), hype_t<-3>{}.template method<7>(0.5, 0.0)));
		TRUE_(check_f<24>(sinh(0.5), hype_t<-3>{}.template method<6>(0.5, 0.0)));
		TRUE_(check_f<22>(sinh(0.5), hype_t<-3>{}.template method<5>(0.5, 0.0)));
		TRUE_(check_f<18>(sinh(0.5), hype_t<-3>{}.template method<4>(0.5, 0.0)));
		TRUE_(check_f<15>(sinh(0.5), hype_t<-3>{}.template method<3>(0.5, 0.0)));
		TRUE_(check_f<11>(sinh(0.5), hype_t<-3>{}.template method<2>(0.5, 0.0)));
		TRUE_(check_f< 9>(sinh(0.5), hype_t<-3>{}.template method<1>(0.5, 0.0)));
		TRUE_(check_f< 2>(sinh(0.5), hype_t<-3>{}.template method<0>(0.5, 0.0)));

		TRUE_(check_f<31>(sinh(1.0), hype_t<-3>{}.template method<9>(1.0, 0.0)));
		TRUE_(check_f<28>(sinh(1.0), hype_t<-3>{}.template method<8>(1.0, 0.0)));
		TRUE_(check_f<22>(sinh(1.0), hype_t<-3>{}.template method<7>(1.0, 0.0)));
		TRUE_(check_f<18>(sinh(1.0), hype_t<-3>{}.template method<6>(1.0, 0.0)));
		TRUE_(check_f<16>(sinh(1.0), hype_t<-3>{}.template method<5>(1.0, 0.0)));
		TRUE_(check_f<13>(sinh(1.0), hype_t<-3>{}.template method<4>(1.0, 0.0)));
		TRUE_(check_f< 9>(sinh(1.0), hype_t<-3>{}.template method<3>(1.0, 0.0)));
		TRUE_(check_f< 6>(sinh(1.0), hype_t<-3>{}.template method<2>(1.0, 0.0)));
		TRUE_(check_f< 4>(sinh(1.0), hype_t<-3>{}.template method<1>(1.0, 0.0)));
		TRUE_(check_f< 1>(sinh(1.0), hype_t<-3>{}.template method<0>(1.0, 0.0)));

	};
	/***/

	EST_("hype_t<-3>::method<~0, 0> (exp(#) - one)")
	{
		T_alpha o{1};
		for (T_sigma i = 0x100; ~--i;) {
			//\
			o *= hype_t<-3, -0>{}.template method<~0, 0>(_fit::mantissa_f(mt19937_f));
			o *= exp(_fit::mantissa_f(mt19937_f)) - one;
		}
		return o;

	};
	/**/
	EST_("hype_t<-3>::method< 8, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= hype_t<-3, -0>{}.template method< 8, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("hype_t<-3>::method< 4, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= hype_t<-3, -0>{}.template method< 4, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("hype_t<-3>::method< 3, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= hype_t<-3, -0>{}.template method< 3, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("hype_t<-3>::method< 2, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= hype_t<-3, -0>{}.template method< 2, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("hype_t<-3>::method< 1, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= hype_t<-3, -0>{}.template method< 1, 0>(_fit::mantissa_f(mt19937_f));
		}
		return o;

	};
	/***/

}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
