#pragma once
#include "./any.cc"
#include "./tangent.hh"




#include "./logarithm.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("logarithm")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;
	static constexpr T_alpha egg =  1.23456789;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	/**/
	TRY_("evaluation")
	{
		TRUE_(check_f<- 8>(logarithm_t< 1, 1>::template method_f< 3>(egg), log(egg)));
		TRUE_(check_f<- 2>(logarithm_t<-1, 1>::template method_f< 3>(egg), exp(egg)));

		TRUE_(check_f<- 1>(logarithm_t< 1   >::template method_f<-1>(egg), log(egg)));
		TRUE_(check_f<- 1>(logarithm_t< 1   >::template method_f< 3>(egg), log(egg)));
		TRUE_(check_f<- 5>(logarithm_t< 1   >::template method_f< 2>(egg), log(egg)));
		TRUE_(check_f<-19>(logarithm_t< 1   >::template method_f< 1>(egg), log(egg)));
		TRUE_(check_f<-40>(logarithm_t< 1   >::template method_f< 0>(egg), log(egg)));

		TRUE_(check_f<- 1>(logarithm_t<-1   >::template method_f<-1>(egg), exp(egg)));
		TRUE_(check_f<- 1>(logarithm_t<-1   >::template method_f< 3>(egg), exp(egg)));
		TRUE_(check_f<-13>(logarithm_t<-1   >::template method_f< 2>(egg), exp(egg)));
		TRUE_(check_f<-29>(logarithm_t<-1   >::template method_f< 1>(egg), exp(egg)));
	//	UNTRUE_(check_f(logarithm_t<-1>::template method_f< 0>(egg), exp(egg)));

		TRUE_(check_f<-48>(logarithm_t< 1, 1>::template method_f<~0>(T_aphex{0.3, 0.8}), logarithm_t< 1, 1>::template method_f< 3>(T_aphex{0.3, 0.8})));
		TRUE_(check_f<-48>(logarithm_t< 1, 1>::template method_f<~0>(T_aphex{0.3, 0.8}), logarithm_t< 1, 1>::template method_f< 0>(T_aphex{0.3, 0.8})));

	}
	/***/
}
TAG_("logarithm trials")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	EST_("logarithm_t<+1,  1; -1>\n   Log@#&\n   (*complex, native*)")
	{
		T_aphex w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			auto y = mt19937_f() + 1.0;
			w += log(T_aphex{x, y});
		//	w += logarithm_t< 1, 1>::template method_f<-1>(T_aphex{x, y});
		}
		return w;
	};
	EST_("logarithm_t<+1,  1;  2>\n   Log@#&\n   (*complex, approx*)")
	{
		T_aphex w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			auto y = mt19937_f() + 1.0;
			w += logarithm_t< 1, 1>::template method_f< 2>(T_aphex{x, y});
		}
		return w;
	};
	EST_("logarithm_t<+1,  0;  2>\n   Log@#&\n   (*complex, approx*)")
	{
		T_aphex w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			auto y = mt19937_f() + 1.0;
			w += logarithm_t< 1, 0>::template method_f< 2>(T_aphex{x, y});
		}
		return w;
	};

	EST_("real std::log")
	{
		T_alpha w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += log(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1, 1>::template method_f<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1>::template method_f<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1>::template method_f<1>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=0>")
	{
		T_alpha w{};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1>::template method_f<0>(x);
		}
		return w;
	};
	EST_("real std::exp")
	{
		T_alpha o{1};
		for (T_sigma i = 100; ~--i;) {
			o *= exp(mt19937_f());
		}
		return o;
	
	};
	EST_("real antilogarithm... <N_lim=~0>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			//\
			w *= logarithm_t<-1>::template method_f<~0>(mt19937_f());
			w *= exp(mt19937_f());
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1, 1>::template method_f<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1>::template method_f<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1>::template method_f<1>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=0>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1>::template method_f<0>(x);
		}
		return w;
	
	};
};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
