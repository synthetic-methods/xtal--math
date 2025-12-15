#pragma once
#include "./any.cc"





#include "./octarithm.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("octarithm")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;
	static constexpr T_alpha egg = 1.23456789;
	static constexpr T_alpha ten = 10;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("octarithm tuning")
	{
		T_alpha constexpr A4_Hz = 440, C4_Hz = 261.625565300598634677849993523304L;
		T_alpha constexpr A4_st =  69, C4_st =  60;

		T_alpha constexpr oct_01 = A4_Hz*octarithm_f<-2, 1,~0>(-_fit::ratio_f(A4_st, 12));
		TRUE_(check_f<-23>(oct_01*octarithm_f<-2, 1, 2>(_fit::ratio_f(C4_st, 12)), C4_Hz));
		TRUE_(check_f<-23>(oct_01*octarithm_f<-2, 1, 3>(_fit::ratio_f(C4_st, 12)), C4_Hz));

		T_alpha constexpr oct_12 = A4_Hz*octarithm_f<-2, 12,~0>(-A4_st);
		TRUE_(check_f<-23>(oct_12*octarithm_f<-2, 12, 3>(C4_st), C4_Hz));
	}
	TRY_("octarithm base-2 evaluation (real)")
	{
		T_alpha o{};

		o = half;
		TRUE_(check_f<- 1>(root_f< 2>(2.0), o = octarithm_t<-2>::template method_f<-1>(o)));
		TRUE_(check_f<- 1>(root_f<-1>(2.0), o = octarithm_t< 2>::template method_f<-1>(o)));

		o = half;
		TRUE_(check_f<-14>(root_f< 2>(2.0), o = octarithm_t<-2>::template method_f< 3>(o)));
		TRUE_(check_f<-14>(root_f<-1>(2.0), o = octarithm_t< 2>::template method_f< 3>(o)));

		o = half;
		TRUE_(check_f<-22>(root_f< 2>(2.0), o = octarithm_t<-2>::template method_f< 2>(o)));
		TRUE_(check_f<-22>(root_f<-1>(2.0), o = octarithm_t< 2>::template method_f< 2>(o)));

		o = half;
		TRUE_(check_f<-30>(root_f< 2>(2.0), o = octarithm_t<-2>::template method_f< 1>(o)));
		TRUE_(check_f<-30>(root_f<-1>(2.0), o = octarithm_t< 2>::template method_f< 1>(o)));

		o = half;
		TRUE_(check_f<-38>(root_f< 2>(2.0), o = octarithm_t<-2>::template method_f< 0>(o)));
		TRUE_(check_f<-38>(root_f<-1>(2.0), o = octarithm_t< 2>::template method_f< 0>(o)));

	};
	TRY_("octarithm base-2 evaluation (integral)")
	{
		TRUE_(check_f<-32>(octarithm_t<-2, 12>::template method_f<0>( 5), 1.3348398541700344));
		TRUE_(check_f<-32>(octarithm_t<-2, 12>::template method_f<0>(17), 2.6696797083400687));

	};

}
TAG_("octarithm trials")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	EST_("real octarithm... <N_lim=~0>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			//\
			w *= octarithm_t<-1>::template method_f<~0>(mt19937_f());
			w *= exp(0.693*mt19937_f());
		}
		return w;
	
	};
	EST_("real octarithm... <N_lim=3>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= octarithm_t<-2>::template method_f<3>(x);
		}
		return w;
	
	};
	EST_("real octarithm... <N_lim=2>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= octarithm_t<-2>::template method_f<2>(x);
		}
		return w;
	
	};
	EST_("real octarithm... <N_lim=1>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= octarithm_t<-2>::template method_f<1>(x);
		}
		return w;
	
	};
	EST_("real octarithm... <N_lim=0>")
	{
		T_alpha w{1};
		for (T_sigma i = 100; ~--i;) {
			auto x = mt19937_f();
			w *= octarithm_t<-2>::template method_f<0>(x);
		}
		return w;
	
	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
