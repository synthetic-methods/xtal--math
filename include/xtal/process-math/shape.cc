#pragma once
#include "./any.cc"
#include "./shape.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("shape")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("shape cardinal")
	{
		TRUE_(check_f<-32>(1.14791942238637400, shape_t<-3, 1>::template function< 4>(0.5, 0.5)));
		TRUE_(check_f<-35>(1.14791942238637400, shape_t<-3, 1>::template function< 3>(0.5, 0.5)));
		TRUE_(check_f<-39>(1.14791942238637400, shape_t<-3, 1>::template function< 2>(0.5, 0.5)));
		TRUE_(check_f<-44>(1.14791942238637400, shape_t<-3, 1>::template function< 1>(0.5, 0.5)));
		TRUE_(check_f<-41>(1.14791942238637400, shape_t<-2, 1>::template function< 4>(0.5, 0.5)));
		TRUE_(check_f<-41>(1.14791942238637400, shape_t<-2, 1>::template function< 3>(0.5, 0.5)));
		TRUE_(check_f<-45>(1.14791942238637400, shape_t<-2, 1>::template function< 2>(0.5, 0.5)));
		TRUE_(check_f<-45>(1.14791942238637400, shape_t<-2, 1>::template function< 1>(0.5, 0.5)));
		TRUE_(check_f<-35>(0.88763432524793617, shape_t< 3, 1>::template function< 4>(0.5, 0.5)));
		TRUE_(check_f<-38>(0.88763432524793617, shape_t< 3, 1>::template function< 3>(0.5, 0.5)));
		TRUE_(check_f<-41>(0.88763432524793617, shape_t< 3, 1>::template function< 2>(0.5, 0.5)));
		TRUE_(check_f<-44>(0.88763432524793617, shape_t< 3, 1>::template function< 1>(0.5, 0.5)));
		TRUE_(check_f<-38>(0.88763432524793617, shape_t< 2, 1>::template function< 4>(0.5, 0.5)));
		TRUE_(check_f<-41>(0.88763432524793617, shape_t< 2, 1>::template function< 3>(0.5, 0.5)));
		TRUE_(check_f<-44>(0.88763432524793617, shape_t< 2, 1>::template function< 2>(0.5, 0.5)));
		TRUE_(check_f<-45>(0.88763432524793617, shape_t< 2, 1>::template function< 1>(0.5, 0.5)));

	}
	TRY_("shape isomorphism")
	{
		TRUE_(check_f<-40>(0.5, shape_t<-3, 0>::template function< 2>(shape_t< 3, 0>::template function< 2>(0.5, 0.5), 0.5)));
		TRUE_(check_f<-40>(0.5, shape_t< 3, 0>::template function< 2>(shape_t<-3, 0>::template function< 2>(0.5, 0.5), 0.5)));

		TRUE_(check_f<- 7>(0.5, shape_t<-2, 0>::template function< 4>(shape_t< 2, 0>::template function< 4>(0.5, 0.5), 0.5)));
		TRUE_(check_f<- 7>(0.5, shape_t< 2, 0>::template function< 4>(shape_t<-2, 0>::template function< 4>(0.5, 0.5), 0.5)));

	}
	/**/
	TRY_("shape evaluation")
	{
		TRUE_(check_f<-2>(0.6180339887498949,  shape_t<-2, -1>::template function<1, 1>(-1.0)));
		TRUE_(check_f<-2>(0.4142135623730949,  shape_t<-2, -1>::template function<1, 1>(-2.0)));
		TRUE_(check_f<-2>(0.3027756377319948,  shape_t<-2, -1>::template function<1, 1>(-3.0)));
		TRUE_(check_f<-2>(0.2360679774997898,  shape_t<-2, -1>::template function<1, 1>(-4.0)));
		TRUE_(check_f<-4>(0.1925824035672532,  shape_t<-2, -1>::template function<1, 1>(-5.0)));
		TRUE_(check_f<-2>(0.1622776601683800,  shape_t<-2, -1>::template function<1, 1>(-6.0)));
		TRUE_(check_f<-2>(0.1400549446402577,  shape_t<-2, -1>::template function<1, 1>(-7.0)));
		TRUE_(check_f<-1>(0.1097722286464435,  shape_t<-2, -1>::template function<1, 1>(-9.0)));

		TRUE_(check_f<-9>(1.6180339887498949,  shape_t<-2, -1>::template function<1, 1>( 1.0)));
		TRUE_(check_f<-9>(2.4142135623730949,  shape_t<-2, -1>::template function<1, 1>( 2.0)));
		TRUE_(check_f<-9>(3.3027756377319948,  shape_t<-2, -1>::template function<1, 1>( 3.0)));
		TRUE_(check_f<-9>(4.2360679774997898,  shape_t<-2, -1>::template function<1, 1>( 4.0)));
		TRUE_(check_f<-9>(5.1925824035672532,  shape_t<-2, -1>::template function<1, 1>( 5.0)));
		TRUE_(check_f<-9>(6.1622776601683800,  shape_t<-2, -1>::template function<1, 1>( 6.0)));
		TRUE_(check_f<-9>(7.1400549446402577,  shape_t<-2, -1>::template function<1, 1>( 7.0)));
		TRUE_(check_f<-9>(9.1097722286464435,  shape_t<-2, -1>::template function<1, 1>( 9.0)));

		TRUE_(check_f<-9>(0.4438171626239681, -shape_t< 3,  0>::template function<~0>(-0.5,-0.5)));
		TRUE_(check_f< 9>(0.4438171626239681, -shape_t< 3,  0>::template function< 9>(-0.5,-0.5)));

		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function<~0>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 9>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 8>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 7>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 6>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 5>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 4>( 0.5, 0.5)));
		TRUE_(check_f< 9>(0.4438171626239681,  shape_t< 3,  0>::template function< 3>( 0.5, 0.5)));
		TRUE_(check_f< 8>(0.4438171626239681,  shape_t< 3,  0>::template function< 2>( 0.5, 0.5)));
		TRUE_(check_f< 7>(0.4438171626239681,  shape_t< 3,  0>::template function< 1>( 0.5, 0.5)));
		TRUE_(check_f< 2>(0.4438171626239681,  shape_t< 3,  0>::template function< 0>( 0.5, 0.5)));
	
		TRUE_(check_f<20>(0.5739597111931870, -shape_t<-3,  0>::template function<~0>(-0.5, -0.5)));
		TRUE_(check_f<20>(0.5739597111931870, -shape_t<-3,  0>::template function< 9>(-0.5, -0.5)));
		TRUE_(check_f<13>(0.5739597111931870, -shape_t<-3,  0>::template function< 2>(-0.5, -0.5)));
		TRUE_(check_f< 8>(0.5739597111931870, -shape_t<-3,  0>::template function< 1>(-0.5, -0.5)));

		TRUE_(check_f<20>(0.57395971119318700, shape_t<-3,  0>::template function<~0>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, shape_t<-3,  0>::template function< 9>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, shape_t<-3,  0>::template function< 8>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, shape_t<-3,  0>::template function< 7>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, shape_t<-3,  0>::template function< 6>( 0.5,  0.5)));
		TRUE_(check_f<20>(0.57395971119318700, shape_t<-3,  0>::template function< 5>( 0.5,  0.5)));
		TRUE_(check_f<19>(0.57395971119318700, shape_t<-3,  0>::template function< 4>( 0.5,  0.5)));
		TRUE_(check_f<17>(0.57395971119318700, shape_t<-3,  0>::template function< 3>( 0.5,  0.5)));
		TRUE_(check_f<13>(0.57395971119318700, shape_t<-3,  0>::template function< 2>( 0.5,  0.5)));
		TRUE_(check_f< 8>(0.57395971119318700, shape_t<-3,  0>::template function< 1>( 0.5,  0.5)));
		TRUE_(check_f< 2>(0.57395971119318700, shape_t<-3,  0>::template function< 0>( 0.5,  0.5)));

		TRUE_(check_f<15>(exp(0.5) - one, shape_t<-3,  0>::template function<9>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, shape_t<-3,  0>::template function<8>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, shape_t<-3,  0>::template function<7>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, shape_t<-3,  0>::template function<6>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, shape_t<-3,  0>::template function<5>(0.5, 1.0)));
		TRUE_(check_f<15>(exp(0.5) - one, shape_t<-3,  0>::template function<4>(0.5, 1.0)));
		TRUE_(check_f<13>(exp(0.5) - one, shape_t<-3,  0>::template function<3>(0.5, 1.0)));
		TRUE_(check_f<11>(exp(0.5) - one, shape_t<-3,  0>::template function<2>(0.5, 1.0)));
		TRUE_(check_f< 5>(exp(0.5) - one, shape_t<-3,  0>::template function<1>(0.5, 1.0)));
		TRUE_(check_f< 1>(exp(0.5) - one, shape_t<-3,  0>::template function<0>(0.5, 1.0)));

		TRUE_(check_f<36>(sinh(0.5), shape_t<-3,  0>::template function<9>(0.5, 0.0)));
		TRUE_(check_f<33>(sinh(0.5), shape_t<-3,  0>::template function<8>(0.5, 0.0)));
		TRUE_(check_f<27>(sinh(0.5), shape_t<-3,  0>::template function<7>(0.5, 0.0)));
		TRUE_(check_f<24>(sinh(0.5), shape_t<-3,  0>::template function<6>(0.5, 0.0)));
		TRUE_(check_f<22>(sinh(0.5), shape_t<-3,  0>::template function<5>(0.5, 0.0)));
		TRUE_(check_f<18>(sinh(0.5), shape_t<-3,  0>::template function<4>(0.5, 0.0)));
		TRUE_(check_f<15>(sinh(0.5), shape_t<-3,  0>::template function<3>(0.5, 0.0)));
		TRUE_(check_f<11>(sinh(0.5), shape_t<-3,  0>::template function<2>(0.5, 0.0)));
		TRUE_(check_f< 9>(sinh(0.5), shape_t<-3,  0>::template function<1>(0.5, 0.0)));
		TRUE_(check_f< 2>(sinh(0.5), shape_t<-3,  0>::template function<0>(0.5, 0.0)));

		TRUE_(check_f<31>(sinh(1.0), shape_t<-3,  0>::template function<9>(1.0, 0.0)));
		TRUE_(check_f<28>(sinh(1.0), shape_t<-3,  0>::template function<8>(1.0, 0.0)));
		TRUE_(check_f<22>(sinh(1.0), shape_t<-3,  0>::template function<7>(1.0, 0.0)));
		TRUE_(check_f<18>(sinh(1.0), shape_t<-3,  0>::template function<6>(1.0, 0.0)));
		TRUE_(check_f<16>(sinh(1.0), shape_t<-3,  0>::template function<5>(1.0, 0.0)));
		TRUE_(check_f<13>(sinh(1.0), shape_t<-3,  0>::template function<4>(1.0, 0.0)));
		TRUE_(check_f< 9>(sinh(1.0), shape_t<-3,  0>::template function<3>(1.0, 0.0)));
		TRUE_(check_f< 6>(sinh(1.0), shape_t<-3,  0>::template function<2>(1.0, 0.0)));
		TRUE_(check_f< 4>(sinh(1.0), shape_t<-3,  0>::template function<1>(1.0, 0.0)));
		TRUE_(check_f< 1>(sinh(1.0), shape_t<-3,  0>::template function<0>(1.0, 0.0)));

	};
	/***/

	EST_("shape_t<-3>::function<~0, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			//\
			o *= shape_t<-3, -0>::template function<~0, 0>(_op::mantissa_f(mt19937_f));
			o *= exp(_op::mantissa_f(mt19937_f)) - one;
		}
		return o;

	};
	EST_("shape_t<-3>::function< 8, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<-3, -0>::template function< 8, 0>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("shape_t<-3>::function< 4, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<-3, -0>::template function< 4, 0>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("shape_t<-3>::function< 3, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<-3, -0>::template function< 3, 0>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("shape_t<-3>::function< 2, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<-3, -0>::template function< 2, 0>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("shape_t<-3>::function< 1, 0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<-3, -0>::template function< 1, 0>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
