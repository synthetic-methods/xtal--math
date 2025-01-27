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
	using _fix = bond::fixture<>;

	using T_sigma = typename _fix::sigma_type;
	using T_delta = typename _fix::delta_type;
	using T_alpha = typename _fix::alpha_type;
	using T_aphex = typename _fix::aphex_type;
	static constexpr T_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _fix::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<6>(tang_t< 1>::template static_method<-1>(0.25), tang_t< 1>::template static_method< 0>(0.25)));
		TRUE_(check_f<7>(tang_t<-1>::template static_method<-1>(0.25), tang_t<-1>::template static_method< 0>(0.25)));
		TRUE_(check_f<7>(tang_t< 2>::template static_method<-1>(0.25), tang_t< 2>::template static_method< 0>(0.25)));
		TRUE_(check_f<7>(tang_t<-2>::template static_method<-1>(0.25), tang_t<-2>::template static_method< 0>(0.25)));

		TRUE_(check_f<- 2>(sinh(_fix::patio_1*0.33), pade::tangy_t<1>::template static_method<~0>(tang_t<2>::template static_method<~0>(0.33))));
		TRUE_(check_f<-46>(sinh(_fix::patio_1*0.33), pade::tangy_t<1>::template static_method< 0>(tang_t<2>::template static_method< 0>(0.33))));
		TRUE_(check_f<-47>(sinh(_fix::patio_1*0.33), pade::tangy_t<1>::template static_method< 1>(tang_t<2>::template static_method< 0>(0.33))));
		TRUE_(check_f<-47>(sinh(_fix::patio_1*0.33), tang_t<1>::template static_method<~0>(tang_t<2>::template static_method< 0>(0.33))*_fix::patio_1));

	};
	EST_("Computing `ArcTan`.")
	{
		T_alpha w{0};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _fix::mantissa_f(mt19937_f);
			w += tang_t<-1>::template static_method<-1>(x);
		}
		return w;
	};
	EST_("Computing `ArcTan` via `tang`.")
	{
		T_alpha w{0};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _fix::mantissa_f(mt19937_f);
			w += tang_t<-1>::template static_method< 0>(x);
		}
		return w;
	};
	EST_("Computing `sinh`.")
	{
		T_alpha w{1};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _fix::mantissa_f(mt19937_f);
			w *= sinh(x);
		}
		return w;
	};
	EST_("Computing `sinh` via `tangy`/`tang` composition.")
	{
		T_alpha w{1};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _fix::mantissa_f(mt19937_f);
			w *= pade::tangy_t<1>::template static_method< 1>(tang_t<2>::template static_method< 0>(x));
		}
		return w;
	};
	EST_("Computing `sinh` via `tang`/`tang` composition.")
	{
		T_alpha w{1};
		for (T_sigma i = 0x10; ~--i;) {
			auto x = _fix::mantissa_f(mt19937_f);
			w *= tang_t<1>::template static_method<~0>(tang_t<2>::template static_method< 0>(x))*_fix::patio_1;
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
