#pragma once
#include "./any.cc"





#include "./octave.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("octave")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;
	static constexpr T_alpha egg =  1.23456789;
	static constexpr T_alpha ten = 10;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("octave tuning")
	{
		T_alpha constexpr frq = 440;
		T_alpha constexpr key = 69;

		T_alpha constexpr oct_01 = frq*octave_f<-2, 1,~0>(-_fit::ratio_f(key, 12));
		TRUE_(check_f<-1>(oct_01*octave_f<-2, 1, 3>(_fit::ratio_f(60, 12)), 261.625565300598634677849993523304L));

		T_alpha constexpr oct_12 = frq*octave_f<-2, 12,~0>(-key);
		TRUE_(check_f<-1>(oct_12*octave_f<-2, 12, 3>(60.), 261.625565300598634677849993523304L));
	}
	TRY_("octave base-2 evaluation")
	{
		T_alpha o{};

		o = half;
		TRUE_(check_f<- 1>(root_f< 2>(2.0), o = octave_t<-2>::template method_f<-1>(o)));
		TRUE_(check_f<- 1>(root_f<-1>(2.0), o = octave_t< 2>::template method_f<-1>(o)));

		o = half;
		TRUE_(check_f<-20>(root_f< 2>(2.0), o = octave_t<-2>::template method_f< 3>(o)));
		TRUE_(check_f<-20>(root_f<-1>(2.0), o = octave_t< 2>::template method_f< 3>(o)));

		o = half;
		TRUE_(check_f<-28>(root_f< 2>(2.0), o = octave_t<-2>::template method_f< 2>(o)));
		TRUE_(check_f<-28>(root_f<-1>(2.0), o = octave_t< 2>::template method_f< 2>(o)));

		o = half;
		TRUE_(check_f<-36>(root_f< 2>(2.0), o = octave_t<-2>::template method_f< 1>(o)));
		TRUE_(check_f<-36>(root_f<-1>(2.0), o = octave_t< 2>::template method_f< 1>(o)));

		o = half;
		TRUE_(check_f<-44>(root_f< 2>(2.0), o = octave_t<-2>::template method_f< 0>(o)));
		TRUE_(check_f<-44>(root_f<-1>(2.0), o = octave_t< 2>::template method_f< 0>(o)));

	};

}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
