#pragma once
#include "./any.cc"
#include "./tang.hh"// testing...

#include "../pade/tangy.hh"



XTAL_ENV_(push)
namespace xtal::process::math::gudermannian::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	static constexpr T_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<6>(tang_t< 1>::template function<-1>(0.25), tang_t< 1>::template function< 0>(0.25)));
		TRUE_(check_f<7>(tang_t<-1>::template function<-1>(0.25), tang_t<-1>::template function< 0>(0.25)));
		TRUE_(check_f<7>(tang_t< 2>::template function<-1>(0.25), tang_t< 2>::template function< 0>(0.25)));
		TRUE_(check_f<7>(tang_t<-2>::template function<-1>(0.25), tang_t<-2>::template function< 0>(0.25)));

		TRUE_(check_f<- 2>(sinh(_op::patio_1*0.33), pade::tangy_t<1>::template function<~0>(tang_t<2>::template function<~0>(0.33))));
		TRUE_(check_f<-46>(sinh(_op::patio_1*0.33), pade::tangy_t<1>::template function< 0>(tang_t<2>::template function< 0>(0.33))));
		TRUE_(check_f<-47>(sinh(_op::patio_1*0.33), pade::tangy_t<1>::template function< 1>(tang_t<2>::template function< 0>(0.33))));
		TRUE_(check_f<-47>(sinh(_op::patio_1*0.33), tang_t<1>::template function<~0>(tang_t<2>::template function< 0>(0.33))*_op::patio_1));

	};
	EST_("Computing `ArcTan`.")
	{
		T_alpha w{0};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f);
			w += tang_t<-1>::template function<-1>(x);
		}
		return w;
	};
	EST_("Computing `ArcTan` via `tang`.")
	{
		T_alpha w{0};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f);
			w += tang_t<-1>::template function< 0>(x);
		}
		return w;
	};
	EST_("Computing `sinh`.")
	{
		T_alpha w{1};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f);
			w *= sinh(x);
		}
		return w;
	};
	EST_("Computing `sinh` via `tangy`/`tang` composition.")
	{
		T_alpha w{1};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f);
			w *= pade::tangy_t<1>::template function< 1>(tang_t<2>::template function< 0>(x));
		}
		return w;
	};
	EST_("Computing `sinh` via `tang`/`tang` composition.")
	{
		T_alpha w{1};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f);
			w *= tang_t<1>::template function<~0>(tang_t<2>::template function< 0>(x))*_op::patio_1;
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
