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

	TRY_("octave base-2 evaluation")
	{
		T_alpha o{};

		o = half;
		TRUE_(check_f<- 1>(root_f< 2>(2.0), o = octave_t<-2, 0>::template method_f<-1>(o)));
		TRUE_(check_f<- 1>(root_f<-1>(2.0), o = octave_t< 2, 0>::template method_f<-1>(o)));

		o = half;
		TRUE_(check_f<-20>(root_f< 2>(2.0), o = octave_t<-2, 0>::template method_f< 3>(o)));
		TRUE_(check_f<-20>(root_f<-1>(2.0), o = octave_t< 2, 0>::template method_f< 3>(o)));

		o = half;
		TRUE_(check_f<-28>(root_f< 2>(2.0), o = octave_t<-2, 0>::template method_f< 2>(o)));
		TRUE_(check_f<-28>(root_f<-1>(2.0), o = octave_t< 2, 0>::template method_f< 2>(o)));

		o = half;
		TRUE_(check_f<-36>(root_f< 2>(2.0), o = octave_t<-2, 0>::template method_f< 1>(o)));
		TRUE_(check_f<-36>(root_f<-1>(2.0), o = octave_t< 2, 0>::template method_f< 1>(o)));

		o = half;
		TRUE_(check_f<-44>(root_f< 2>(2.0), o = octave_t<-2, 0>::template method_f< 0>(o)));
		TRUE_(check_f<- 1>(root_f<-1>(2.0), o = octave_t< 2, 0>::template method_f< 0>(o)));

	};
	TRY_("octave base-1 evaluation")
	{
		T_aphex o{0.5};

		TRUE_(check_f<-1>(-1.0, o = octave_t<-1, 0>::template method_f<-1>(o)));
		TRUE_(check_f<-1>( 0.5, o = octave_t< 1, 0>::template method_f<-1>(o)));

	};

}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
