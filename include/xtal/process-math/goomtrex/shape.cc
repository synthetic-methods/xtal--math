#pragma once
#include "./any.cc"
#include "./shape.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::goomtrex::_test
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

	TRY_("shape evaluation")
	{
		TRUE_(check_f<30>(exp(0.5) - one, shape_t<3, 0>::template function<7>(0.5, 1.0)));
		TRUE_(check_f<25>(exp(0.5) - one, shape_t<3, 0>::template function<6>(0.5, 1.0)));
		TRUE_(check_f<23>(exp(0.5) - one, shape_t<3, 0>::template function<5>(0.5, 1.0)));
		TRUE_(check_f<19>(exp(0.5) - one, shape_t<3, 0>::template function<4>(0.5, 1.0)));
		TRUE_(check_f<17>(exp(0.5) - one, shape_t<3, 0>::template function<3>(0.5, 1.0)));
		TRUE_(check_f<13>(exp(0.5) - one, shape_t<3, 0>::template function<2>(0.5, 1.0)));
		TRUE_(check_f< 7>(exp(0.5) - one, shape_t<3, 0>::template function<1>(0.5, 1.0)));
		TRUE_(check_f< 1>(exp(0.5) - one, shape_t<3, 0>::template function<0>(0.5, 1.0)));

		TRUE_(check_f<27>(sinh(0.5), shape_t<3, 0>::template function<7>(0.5, 0.0)));
		TRUE_(check_f<24>(sinh(0.5), shape_t<3, 0>::template function<6>(0.5, 0.0)));
		TRUE_(check_f<22>(sinh(0.5), shape_t<3, 0>::template function<5>(0.5, 0.0)));
		TRUE_(check_f<18>(sinh(0.5), shape_t<3, 0>::template function<4>(0.5, 0.0)));
		TRUE_(check_f<15>(sinh(0.5), shape_t<3, 0>::template function<3>(0.5, 0.0)));
		TRUE_(check_f<11>(sinh(0.5), shape_t<3, 0>::template function<2>(0.5, 0.0)));
		TRUE_(check_f< 9>(sinh(0.5), shape_t<3, 0>::template function<1>(0.5, 0.0)));
		TRUE_(check_f< 2>(sinh(0.5), shape_t<3, 0>::template function<0>(0.5, 0.0)));

		TRUE_(check_f<22>(sinh(1.0), shape_t<3, 0>::template function<7>(1.0, 0.0)));
		TRUE_(check_f<18>(sinh(1.0), shape_t<3, 0>::template function<6>(1.0, 0.0)));
		TRUE_(check_f<16>(sinh(1.0), shape_t<3, 0>::template function<5>(1.0, 0.0)));
		TRUE_(check_f<13>(sinh(1.0), shape_t<3, 0>::template function<4>(1.0, 0.0)));
		TRUE_(check_f< 9>(sinh(1.0), shape_t<3, 0>::template function<3>(1.0, 0.0)));
		TRUE_(check_f< 6>(sinh(1.0), shape_t<3, 0>::template function<2>(1.0, 0.0)));
		TRUE_(check_f< 4>(sinh(1.0), shape_t<3, 0>::template function<1>(1.0, 0.0)));
		TRUE_(check_f< 1>(sinh(1.0), shape_t<3, 0>::template function<0>(1.0, 0.0)));

	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function<~0>(0.5, 1.0), 0.64872127070012819));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 9>(0.5, 1.0), 0.64869025019892934));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 8>(0.5, 1.0), 0.64172208700944222));

	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function<~0>(1.0, 1.0), 1.7182818284590451));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 9>(1.0, 1.0), 1.7170793981481478));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 8>(1.0, 1.0), 1.7111788346357946));

	};
	EST_("shape_t<3>::function<  ~0>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			//\
			o *= shape_t<3, -0>::template function<~0>(_op::mantissa_f(mt19937_f));
			o *= exp(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("shape_t<3>::function<  9>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<3, -0>::template function< 9>(_op::mantissa_f(mt19937_f), 1.0);
		}
		return o;

	};
	EST_("shape_t<3>::function<  9>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<3, -0>::template function< 3>(_op::mantissa_f(mt19937_f), 1.0);
		}
		return o;

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
