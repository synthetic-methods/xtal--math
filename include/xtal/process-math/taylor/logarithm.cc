#pragma once
#include "./any.cc"
#include "./logarithm.hh"// testing...

#include "./monologarithm.hh"
#include "../pade/unity.hh"



XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("logarithm")
{
	using _std::log;
	using _std::exp;

	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	static constexpr T_alpha egg =  1.23456789;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = algebra::d_::circular_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("inversion")
	{
		TRUE_(check_f<-1>(one, logarithm_t< 1,  1>::template function<0>(egg)*logarithm_t< 1, -1>::template function<0>(egg)));
		TRUE_(check_f<-1>(one, logarithm_t<-1,  1>::template function<0>(egg)*logarithm_t<-1, -1>::template function<0>(egg)));

		TRUE_(check_f<19>(egg, logarithm_t< 1,  1>::template function<0>(logarithm_t<-1,  1>::template function<0>(egg))));
		TRUE_(check_f<19>(egg, logarithm_t<-1,  1>::template function<0>(logarithm_t< 1,  1>::template function<0>(egg))));

	}
	TRY_("evaluation")
	{
		TRUE_(check_f<-13>(logarithm_t< 1, 1, 1>::template function<2>(egg), log(egg)));
		TRUE_(check_f<- 4>(logarithm_t<-1, 1, 1>::template function<2>(egg), exp(egg)));

		TRUE_(check_f<- 1>(logarithm_t< 1>::template function<-1>(egg), log(egg)));
		TRUE_(check_f<- 5>(logarithm_t< 1>::template function< 3>(egg), log(egg)));
		TRUE_(check_f<-19>(logarithm_t< 1>::template function< 2>(egg), log(egg)));
		TRUE_(check_f<-33>(logarithm_t< 1>::template function< 1>(egg), log(egg)));
		TRUE_(check_f<-40>(logarithm_t< 1>::template function< 0>(egg), log(egg)));

		TRUE_(check_f<- 1>(logarithm_t<-1>::template function<-1>(egg), exp(egg)));
		TRUE_(check_f<-14>(logarithm_t<-1>::template function< 3>(egg), exp(egg)));
		TRUE_(check_f<-29>(logarithm_t<-1>::template function< 2>(egg), exp(egg)));
		TRUE_(check_f<-45>(logarithm_t<-1>::template function< 1>(egg), exp(egg)));
	//	UNTRUE_(check_f(logarithm_t<-1>::template function< 0>(egg), exp(egg)));

	}

	EST_("real std::log")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += _std::log(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += logarithm_t< 1, 1, 1>::template function<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += logarithm_t< 1>::template function<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += logarithm_t< 1>::template function<1>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=0>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += logarithm_t< 1>::template function<0>(x);
		}
		return w;
	};

	EST_("real std::exp")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - one;
			w *= _std::exp(x);
		}
		return w;
	};
	EST_("real antilogarithm... <N_lim=2, M_car=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - one;
			w *= logarithm_t<-1, 1, 1>::template function<2>(x);
		}
		return w;
	};
	EST_("real antilogarithm... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - one;
			w *= logarithm_t<-1>::template function<2>(x);
		}
		return w;
	};
	EST_("real antilogarithm... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - one;
			w *= logarithm_t<-1>::template function<1>(x);
		}
		return w;
	};
	EST_("real antilogarithm... <N_lim=0>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) - one;
			w *= logarithm_t<-1>::template function<0>(x);
		}
		return w;
	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
