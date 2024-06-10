#pragma once
#include "./any.cc"
#include "./tangent.hh"// testing...






XTAL_ENV_(push)
namespace xtal::process::math::gudermannian::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_t;
	using T_delta = typename _op::delta_t;
	using T_alpha = typename _op::alpha_t;
	using T_aphex = typename _op::aphex_t;
	static constexpr T_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = algebra::d_::circular_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		using _std:: tan;
		using _std::atan;

		TRUE_(check_f<6>(tangent_t< 1>::template function<-1>(0.25), tangent_t< 1>::template function< 0>(0.25)));
		TRUE_(check_f<7>(tangent_t<-1>::template function<-1>(0.25), tangent_t<-1>::template function< 0>(0.25)));
		TRUE_(check_f<7>(tangent_t< 2>::template function<-1>(0.25), tangent_t< 2>::template function< 0>(0.25)));
		TRUE_(check_f<7>(tangent_t<-2>::template function<-1>(0.25), tangent_t<-2>::template function< 0>(0.25)));
		
//		using foo = process::confined_t<dilated<1>, tangent< 2>>;
//		using bar = process::confined_t<dilated<1>, tangent<-2>>;
//
//		echo(foo::function(bar::function(0.50) + bar::function(0.75)));
//		echo(_std::tanh(_std::atanh(0.50) + _std::atanh(0.75)));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
