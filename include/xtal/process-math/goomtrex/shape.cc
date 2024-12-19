#pragma once
#include "./any.cc"
#include "./shape.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::goomtrex::_test
{/////////////////////////////////////////////////////////////////////////////////FIXME
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
	static constexpr T_alpha zero =  0.0;
	static constexpr T_alpha half =  0.5;
	static constexpr T_alpha  one =  1.0;
	static constexpr T_alpha  two =  2.0;
	static constexpr T_alpha  ten = 10.0;

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("shape evaluation")
	{
		echo();
		echo(exp(half));
		echo(shape_t<3, 0>::template function<7>(half, one));
		echo(shape_t<3, 0>::template function<6>(half, one));
		echo(shape_t<3, 0>::template function<5>(half, one));
		echo(shape_t<3, 0>::template function<4>(half, one));
		echo(shape_t<3, 0>::template function<3>(half, one));
		echo(shape_t<3, 0>::template function<2>(half, one));
		echo(shape_t<3, 0>::template function<1>(half, one));
		echo(shape_t<3, 0>::template function<0>(half, one));

		TRUE_(check_f<27>(sinh(half), shape_t<3, 0>::template function<7>(half, zero)));
		TRUE_(check_f<24>(sinh(half), shape_t<3, 0>::template function<6>(half, zero)));
		TRUE_(check_f<22>(sinh(half), shape_t<3, 0>::template function<5>(half, zero)));
		TRUE_(check_f<18>(sinh(half), shape_t<3, 0>::template function<4>(half, zero)));
		TRUE_(check_f<15>(sinh(half), shape_t<3, 0>::template function<3>(half, zero)));
		TRUE_(check_f<11>(sinh(half), shape_t<3, 0>::template function<2>(half, zero)));
		TRUE_(check_f< 9>(sinh(half), shape_t<3, 0>::template function<1>(half, zero)));
		TRUE_(check_f< 2>(sinh(half), shape_t<3, 0>::template function<0>(half, zero)));

		TRUE_(check_f<22>(sinh(one), shape_t<3, 0>::template function<7>(one, zero)));
		TRUE_(check_f<18>(sinh(one), shape_t<3, 0>::template function<6>(one, zero)));
		TRUE_(check_f<16>(sinh(one), shape_t<3, 0>::template function<5>(one, zero)));
		TRUE_(check_f<13>(sinh(one), shape_t<3, 0>::template function<4>(one, zero)));
		TRUE_(check_f< 9>(sinh(one), shape_t<3, 0>::template function<3>(one, zero)));
		TRUE_(check_f< 6>(sinh(one), shape_t<3, 0>::template function<2>(one, zero)));
		TRUE_(check_f< 4>(sinh(one), shape_t<3, 0>::template function<1>(one, zero)));
		TRUE_(check_f< 1>(sinh(one), shape_t<3, 0>::template function<0>(one, zero)));

	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function<~0>(half, one), 0.64872127070012819));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 9>(half, one), 0.64869025019892934));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 8>(half, one), 0.64172208700944222));

	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function<~0>(one, one), 1.7182818284590451));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 9>(one, one), 1.7170793981481478));
	//	TRUE_(check_f<-3>(shape_t<3, -0>::template function< 8>(one, one), 1.7111788346357946));

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
			o *= shape_t<3, -0>::template function< 9>(_op::mantissa_f(mt19937_f), one);
		}
		return o;

	};
	EST_("shape_t<3>::function<  9>")
	{
		T_alpha o{1};

		for (int i{0}; i < 0x100; ++i) {
			o *= shape_t<3, -0>::template function< 3>(_op::mantissa_f(mt19937_f), one);
		}
		return o;

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
