#pragma once
#include "./any.cc"
#include "./polynomial.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::bernoulli::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;
	static constexpr T_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<11>(polynomial_t<1, 0>::template method_f<-1>(0.250), abs(polynomial_t<1, 0>::template method_f< 7>(0.250))));
		TRUE_(check_f< 8>(polynomial_t<1, 0>::template method_f<-1>(0.250), abs(polynomial_t<1, 0>::template method_f< 5>(0.250))));
		TRUE_(check_f< 5>(polynomial_t<1, 0>::template method_f<-1>(0.250), abs(polynomial_t<1, 0>::template method_f< 3>(0.250))));
		
		TRUE_(check_f<10>(polynomial_t<1, 0>::template method_f<-2>(0.125), abs(polynomial_t<1, 0>::template method_f< 6>(0.125))));
		TRUE_(check_f< 7>(polynomial_t<1, 0>::template method_f<-2>(0.125), abs(polynomial_t<1, 0>::template method_f< 4>(0.125))));
		TRUE_(check_f< 2>(polynomial_t<1, 0>::template method_f<-2>(0.125), abs(polynomial_t<1, 0>::template method_f< 2>(0.125))));
		
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f<-1>( 0.25), polynomial_t<1, 1>::template method_f<-1>( 1.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 7>( 0.25), polynomial_t<1, 1>::template method_f< 7>( 1.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 5>( 0.25), polynomial_t<1, 1>::template method_f< 5>( 1.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 3>( 0.25), polynomial_t<1, 1>::template method_f< 3>( 1.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 1>( 0.25), polynomial_t<1, 1>::template method_f< 1>( 1.25)));
		
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f<-2>( 0.25), polynomial_t<1, 1>::template method_f<-2>(-0.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 6>( 0.25), polynomial_t<1, 1>::template method_f< 6>(-0.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 4>( 0.25), polynomial_t<1, 1>::template method_f< 4>(-0.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 2>( 0.25), polynomial_t<1, 1>::template method_f< 2>(-0.25)));
		TRUE_(check_f<-1>(polynomial_t<1, 0>::template method_f< 0>( 0.25), polynomial_t<1, 1>::template method_f< 0>(-0.25)));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
