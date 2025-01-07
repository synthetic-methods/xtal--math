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
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	static constexpr T_alpha egg =  1.23456789;
	static constexpr T_alpha ten = 10;

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(check_f<-31>(logarithm_t< 1, 1>::template function<~0>(T_aphex{0.3, 0.8}), logarithm_t< 1, 1>::template function< 2>(T_aphex{0.3, 0.8})));
		TRUE_(check_f<-46>(logarithm_t< 1, 1>::template function<~0>(T_aphex{0.3, 0.8}), logarithm_t< 1, 1>::template function< 0>(T_aphex{0.3, 0.8})));

		TRUE_(check_f<-13>(logarithm_t< 1, 1>::template function<2>(egg), log(egg)));
		TRUE_(check_f<-25>(logarithm_t<-1, 1>::template function<2>(egg), exp(egg)));

		TRUE_(check_f<- 1>(logarithm_t< 1>::template function<-1>(egg), log(egg)));
		TRUE_(check_f<- 5>(logarithm_t< 1>::template function< 3>(egg), log(egg)));
		TRUE_(check_f<-19>(logarithm_t< 1>::template function< 2>(egg), log(egg)));
		TRUE_(check_f<-33>(logarithm_t< 1>::template function< 1>(egg), log(egg)));
		TRUE_(check_f<-40>(logarithm_t< 1>::template function< 0>(egg), log(egg)));

		TRUE_(check_f<- 1>(logarithm_t<-1>::template function<-1>(egg), exp(egg)));
		TRUE_(check_f<-26>(logarithm_t<-1>::template function< 3>(egg), exp(egg)));
		TRUE_(check_f<-34>(logarithm_t<-1>::template function< 2>(egg), exp(egg)));
		TRUE_(check_f<-43>(logarithm_t<-1>::template function< 1>(egg), exp(egg)));
	//	UNTRUE_(check_f(logarithm_t<-1>::template function< 0>(egg), exp(egg)));

	}
	/**/
	EST_("complex std::log")
	{
		T_aphex w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			auto y = _op::mantissa_f(mt19937_f) + 1.0;
			w += log(T_aphex{x, y});
		}
		return w;
	};
	EST_("complex logarithm... <N_lim=2, M_car=1>")
	{
		T_aphex w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			auto y = _op::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1, 1>::template function<2>(T_aphex{x, y});
		}
		return w;
	};

	EST_("real std::log")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			w += log(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1, 1>::template function<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1>::template function<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1>::template function<1>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=0>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + 1.0;
			w += logarithm_t< 1>::template function<0>(x);
		}
		return w;
	};
	/***/
	EST_("real std::exp")
	{
		T_alpha o{1};
		for (T_sigma i = 0x100; ~--i;) {
			o *= exp(_op::mantissa_f(mt19937_f));
		}
		return o;
	
	};
	/**/
	EST_("real antilogarithm... <N_lim=~0>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			w *= logarithm_t<-1>::template function<~0>(_op::mantissa_f(mt19937_f));
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1, 1>::template function<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1>::template function<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1>::template function<1>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=0>")
	{
		T_alpha w{1};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - 1.0;
			w *= logarithm_t<-1>::template function<0>(x);
		}
		return w;
	
	};
	/***/

}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
