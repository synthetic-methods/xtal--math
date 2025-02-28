#pragma once
#include "./any.cc"
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

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(check_f<-31>(logarithm_t< 1, 1>::template method_f<~0>(T_aphex{0.3, 0.8}), logarithm_t< 1, 1>::template method_f< 2>(T_aphex{0.3, 0.8})));
		TRUE_(check_f<-46>(logarithm_t< 1, 1>::template method_f<~0>(T_aphex{0.3, 0.8}), logarithm_t< 1, 1>::template method_f< 0>(T_aphex{0.3, 0.8})));

		TRUE_(check_f<-13>(logarithm_t< 1, 1>::template method_f< 2>(egg), log(egg)));
		TRUE_(check_f<-26>(logarithm_t<-1, 1>::template method_f< 2>(egg), exp(egg)));

		TRUE_(check_f<- 1>(logarithm_t< 1   >::template method_f<-1>(egg), log(egg)));
		TRUE_(check_f<- 5>(logarithm_t< 1   >::template method_f< 3>(egg), log(egg)));
		TRUE_(check_f<-19>(logarithm_t< 1   >::template method_f< 2>(egg), log(egg)));
		TRUE_(check_f<-33>(logarithm_t< 1   >::template method_f< 1>(egg), log(egg)));
		TRUE_(check_f<-41>(logarithm_t< 1   >::template method_f< 0>(egg), log(egg)));

		TRUE_(check_f<- 1>(logarithm_t<-1   >::template method_f<-1>(egg), exp(egg)));
		TRUE_(check_f<-27>(logarithm_t<-1   >::template method_f< 3>(egg), exp(egg)));
		TRUE_(check_f<-35>(logarithm_t<-1   >::template method_f< 2>(egg), exp(egg)));
		TRUE_(check_f<-43>(logarithm_t<-1   >::template method_f< 1>(egg), exp(egg)));
	//	UNTRUE_(check_f(logarithm_t<-1>::template method_f< 0>(egg), exp(egg)));

	//	_std::array<T_alpha, 0x20> food{};
	//	for (int i{}; i < food.size(); ++i) {
	//		auto x = _fit::ratio_f(i - 0x10, 0x10);
	//		x  = taylor::logarithm_t<-1>::template method_f<0>(x);
	//		x -= one;
	//		x *= _fit::ratio_f(1, 3);
	//		food[i] = x;
	//	}
	//	echo_plot_<49>(food);
	}
	/**/
	TRY_("logarithm base-2")
	{
		T_alpha o{};

		o = half;
		TRUE_(check_f<- 1>(root_f< 2>(2.0), o = logarithm_t<-1, 0, 2>::template method_f<-1>(o)));
		TRUE_(check_f<- 1>(root_f<-1>(2.0), o = logarithm_t< 1, 0, 2>::template method_f<-1>(o)));

		o = half;
		TRUE_(check_f<-20>(root_f< 2>(2.0), o = logarithm_t<-1, 0, 2>::template method_f< 3>(o)));
		TRUE_(check_f<-20>(root_f<-1>(2.0), o = logarithm_t< 1, 0, 2>::template method_f< 3>(o)));

		o = half;
		TRUE_(check_f<-28>(root_f< 2>(2.0), o = logarithm_t<-1, 0, 2>::template method_f< 2>(o)));
		TRUE_(check_f<-28>(root_f<-1>(2.0), o = logarithm_t< 1, 0, 2>::template method_f< 2>(o)));

		o = half;
		TRUE_(check_f<-36>(root_f< 2>(2.0), o = logarithm_t<-1, 0, 2>::template method_f< 1>(o)));
		TRUE_(check_f<-36>(root_f<-1>(2.0), o = logarithm_t< 1, 0, 2>::template method_f< 1>(o)));

		o = half;
		TRUE_(check_f<-44>(root_f< 2>(2.0), o = logarithm_t<-1, 0, 2>::template method_f< 0>(o)));
		TRUE_(check_f<- 1>(root_f<-1>(2.0), o = logarithm_t< 1, 0, 2>::template method_f< 0>(o)));

	};
	TRY_("logarithm base-1")
	{
		T_aphex o{0.5};

		TRUE_(check_f<-1>(-1.0, o = logarithm_t<-1, 0, 1>::template method_f<-1>(o)));
		TRUE_(check_f<-1>( 0.5, o = logarithm_t< 1, 0, 1>::template method_f<-1>(o)));

	};


	EST_("complex std::log")
	{
		T_aphex w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			auto y = _fit::mantissa_f(mt19937_f) + 1.0;
			w += log(T_aphex{x, y});
		}
		return w;
	};
	EST_("complex logarithm... <N_lim=2, M_car=1>")
	{
		T_aphex w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			auto y = _fit::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1, 1>::template method_f<2>(T_aphex{x, y});
		}
		return w;
	};

	EST_("real std::log")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			w += log(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1, 1>::template method_f<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1>::template method_f<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1>::template method_f<1>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=0>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1>::template method_f<0>(x);
		}
		return w;
	};
	/***/
	EST_("real std::exp")
	{
		T_alpha o{1};
		for (T_sigma i = 0x100; ~--i;) {
			o *= exp(_fit::mantissa_f(mt19937_f));
		}
		return o;
	
	};
	/**/
	EST_("real antilogarithm... <N_lim=~0>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			w *= logarithm_t<-1>::template method_f<~0>(_fit::mantissa_f(mt19937_f));
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1, 1>::template method_f<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1>::template method_f<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1>::template method_f<1>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=0>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1>::template method_f<0>(x);
		}
		return w;
	
	};
	/***/

}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
